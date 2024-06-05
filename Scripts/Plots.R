# make plots for EDA
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(viridis)
library(genefilter)
library(ComplexHeatmap)
library(data.table)
library(cowplot)
library(InteractiveComplexHeatmap)
library("glmpca")
library(EnhancedVolcano)

# perform lfc shrinkage for better visualisation
KO_df_shrunk <- lfcShrink(dds_KO,contrast = c("sample","shRNA","Control"),type = "ashr")
cTreatment_df_shrunk <- lfcShrink(dds_treatment,contrast =  Dox_control - Ifn_control, type="ashr")
KO_df_shrunk <- lfcShrink(dds_treatment,contrast =  Dox_shRNA - Ifn_shRNA, type="ashr")

# generate pca plots for control vs shRNA
## create transformed values
vsd <- vst(dds_KO,blind=FALSE)
head(assay(vsd),3)
select <- order(rowMeans(counts(dds_KO,normalized=TRUE)),
                decreasing=TRUE)[1:40]

#$PCA Plots
#make pca plot with vst data
pcaData <- plotPCA(vsd,intgroup=c("sample"),returnData=T)
plot(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))
vst_pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with Control vs shRNA")+
  theme_cowplot()
# GLM-PCA plots
tiff("Plots/KO_GPCA_plot_vst.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(vst_pca_plot)
dev.off()

gpca <- glmpca(counts(dds_KO), L=2)
gpca.dat <- gpca$factors
gpca.dat$sample <- dds_KO$sample

vst_pca_plot <- ggplot(gpca.dat, aes(x = dim1, y = dim2, color=sample)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")+
  theme_cowplot()

#make PCA plot with Control vs Treatment
vsd <- vst(dds_treatment,blind=FALSE)
#rld <- rlog(dds_KO,blind=F)
head(assay(vsd),3)
select <- order(rowMeans(counts(dds_treatment,normalized=TRUE)),
                decreasing=TRUE)[1:40]

#$PCA Plots
#make pca plot with vst data
pcaData <- plotPCA(vsd,intgroup=c("sample","condition"),returnData=T)
plot(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))
vst_pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, shape = sample, color=condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA Ifn vs Dox")+
  theme_cowplot()
# GLM-PCA plots
tiff("Plots/Treatment_GPCA_plot_vst.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(vst_pca_plot)
dev.off()

gpca <- glmpca(counts(dds_treatment), L=2)
gpca.dat <- gpca$factors
gpca.dat$sample <- dds_treatment$sample
gpca.dat$condition <- dds_treatment$condition


vst_pca_plot <- ggplot(gpca.dat, aes(x = dim1, y = dim2, color = condition, shape=sample)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")+
  theme_cowplot()




#adjusting normalized counts
#rownames(normalized_counts_KO) <- normalized_counts_KO$GeneID

# Making heatmaps for KO vs control

padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold <- KO_df$padj < padj.cutoff & abs(KO_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
KO_df$threshold <- threshold
KO_sig <- data.frame(subset(KO_df,threshold==TRUE))
KO_sig <- KO_sig[order(KO_sig$padj),]
only_KO <- normalized_counts_KO[rownames(KO_sig),]
only_KO$padj <- KO_sig[match(rownames(only_KO),rownames(KO_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_KO_mat <- as.matrix(only_KO[,1:7])
only_KO_mat_scaled <- t(scale(t(only_KO[,1:7])))
rownames(only_KO_mat_scaled) <- only_KO$GeneName
#only_KO_mat_scaled <- cbind(only_KO_mat_scaled,only_KO[match(rownames(only_KO_mat_scaled),only_KO$GeneName),6])
ht_KO <-ComplexHeatmap::Heatmap(only_KO_mat_scaled,
                                 column_names_side = "bottom",
                                 cluster_rows = T,
                                 cluster_columns = F,
                                 column_names_rot = 45,
                                 column_names_gp = gpar(fontsize=7),
                                 row_names_gp = gpar(fontsize=7),
                                row_names_rot = 0,
                                name = "Z-score",width = 4*unit(25,'mm'),
                                column_title = "Control vs shRNA")
ht_KO
tiff("Plots/Control_vs_shRNA.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_KO)
dev.off()

#making heatmap for Control ifn vs Dox

threshold <- cTreatment_df$padj < padj.cutoff & abs(cTreatment_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
cTreatment_df$threshold <- threshold
cTreatment_sig <- data.frame(subset(cTreatment_df,threshold==TRUE))
cTreatment_sig <- cTreatment_sig[order(cTreatment_sig$padj),]
only_cT <- normalized_counts_tr[rownames(cTreatment_sig),]
only_cT$padj <- cTreatment_sig[match(rownames(only_cT),rownames(cTreatment_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_cT_mat <- as.matrix(only_cT[,1:4])
only_cT_mat_scaled <- t(scale(t(only_cT[,1:4])))
rownames(only_cT_mat_scaled) <- only_cT$GeneName
#only_KO_mat_scaled <- cbind(only_KO_mat_scaled,only_KO[match(rownames(only_KO_mat_scaled),only_KO$GeneName),6])

ht_cT <-ComplexHeatmap::Heatmap(only_cT_mat_scaled,
                                column_names_side = "bottom",
                                cluster_rows = T,
                                cluster_columns = F,
                                column_names_rot = 45,
                                column_names_gp = gpar(fontsize=7),
                                row_names_gp = gpar(fontsize=7),
                                row_names_rot = 0,
                                name = "Z-score",width = 4*unit(25,'mm'),
                                column_title = "Control Ifn vs Dox")
ht_cT
tiff("Plots/Control_Ifn_vs_Dox.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_cT)
dev.off()

#making heatmap for shRNA Ifn vs Dox
threshold <- KOTreatment_df$padj < padj.cutoff & abs(KOTreatment_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
KOTreatment_df$threshold <- threshold
KOTreatment_sig <- data.frame(subset(KOTreatment_df,threshold==TRUE))
KOTreatment_sig <- KOTreatment_sig[order(KOTreatment_sig$padj),]
only_KT <- normalized_counts_tr[rownames(KOTreatment_sig),]
only_KT$padj <- KOTreatment_sig[match(rownames(only_KT),rownames(KOTreatment_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_KT_mat <- as.matrix(only_KT[,5:12])
only_KT_mat_scaled <- t(scale(t(only_KT[,5:12])))
rownames(only_KT_mat_scaled) <- only_KT$GeneName
#only_KO_mat_scaled <- cbind(only_KO_mat_scaled,only_KO[match(rownames(only_KO_mat_scaled),only_KO$GeneName),6])
ht_KT <-ComplexHeatmap::Heatmap(only_KT_mat_scaled,
                                column_names_side = "bottom",
                                cluster_rows = T,
                                cluster_columns = F,
                                column_names_rot = 45,
                                column_names_gp = gpar(fontsize=7),
                                row_names_gp = gpar(fontsize=7),
                                row_names_rot = 0,
                                name = "Z-score",width = 4*unit(25,'mm'),
                                column_title = "shRNA Ifn vs Dox")
ht_KT
tiff("Plots/shRNA_Ifn_vs_Dox.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_KT)
dev.off()
#!--------------------------------------------------------#


#plots to make after enrichment analysis

as.matrix(na.omit(data[Wnt_genes$human_gene_symbol,c(1,2,5,6,3,4,7,8)]))
###
#make heatmap for wnt genes
all_KT <- normalized_counts_tr
all_KT_mat <- as.matrix(all_KT[,5:12])
all_KT_mat_scaled <- t(scale(t(all_KT[,5:12])))
rownames(all_KT_mat_scaled) <- all_KT$GeneName
data <- as.data.frame(only_KT_mat_scaled)

col_fun = colorRamp2(c(-2,-1, 0,1, 2), c("#327eba","#6d7baa", "#967796","#b57283", "#e06663"),space = "XYZ")
col_fun(seq(-3, 3))
ht_KT_bc <-ComplexHeatmap::Heatmap(as.matrix(na.omit(data[Ifn_gamma_genes$human_gene_symbol,c(1,2,5,6,3,4,7,8)])),
                                column_names_side = "bottom",
                                cluster_rows = T,
                                cluster_columns = F,
                                column_names_rot = 45,
                                column_names_gp = gpar(fontsize=7),
                                row_names_gp = gpar(fontsize=6),
                                row_names_rot = 0,
                                name = " ",width = 4*unit(25,'mm'),
                                column_title = "Ifn-gamma pathway genes shRNA Ifn vs Dox",
                                col = col_fun)
  set_enrichplot_color(type = "fill")
ht_KT_bc
tiff("Plots/new_Ifn_pathway_KDIfn_vs_Dox.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_KT_bc)
dev.off()



#making a volcano plot
E1 <- EnhancedVolcano(KO_df,
                      lab = KO_df$GeneName,
                      x='log2FoldChange',
                      y='padj',
                      title = 'shRNA vs control',
                      subtitle = 'Differential Expression',
                      caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      labSize = 3,
                      pointSize = 0.7,
                      colAlpha = 0.8,
                      legendPosition = 'top',
                      legendIconSize = 3,
                      drawConnectors = F,
                      widthConnectors = 0.5)+
  coord_cartesian(xlim = c(-10,10))
tiff("Plots/Control_vs_shRNA_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E1)
dev.off()




#making a volcano plot for control Ifn vs Dox
E2 <- EnhancedVolcano(cTreatment_df,
                      lab = cTreatment_df$GeneName,
                      x='log2FoldChange',
                      y='padj',
                      title = 'Control Ifn vs Dox',
                      subtitle = 'Differential Expression',
                      caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      labSize = 3,
                      pointSize = 0.7,
                      colAlpha = 0.8,
                      legendPosition = 'top',
                      legendIconSize = 3,
                      drawConnectors = F,
                      widthConnectors = 0.5)+
  coord_cartesian(xlim = c(-5,5))
tiff("Plots/Control_Ifn_vs_Dox_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E2)
dev.off()


#making a volcano plot
E3 <- EnhancedVolcano(KOTreatment_df,
                      lab = KOTreatment_df$GeneName,
                      x='log2FoldChange',
                      y='padj',
                      title = 'shRNA Ifn vs Dox',
                      subtitle = 'Differential Expression',
                      caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      labSize = 3,
                      pointSize = 0.7,
                      colAlpha = 0.8,
                      legendPosition = 'top',
                      legendIconSize = 3,
                      drawConnectors = F,
                      widthConnectors = 0.5)+
  coord_cartesian(xlim = c(-10,10))
tiff("Plots/shRNA_Ifn_vs_Dox_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E3)
dev.off()





#making heatmaps with all the columns/samples
all_d5 <- normalized_counts[rownames(d5_sig),]
all_d5$padj <- d5_sig[match(rownames(all_d5),rownames(d5_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
#only_d5_mat <- as.matrix(only_low[,c(2,3,4,5)])
all_d5_mat_scaled <- t(scale(t(all_d5[,c(2:6)])))
rownames(all_d5_mat_scaled) <- all_d5$GeneName

all_d5_mat_scaled <- cbind(all_d5_mat_scaled,all_d5[match(rownames(all_d5_mat_scaled),all_d5$GeneName),8])
ht_d5_all <-ComplexHeatmap::Heatmap(na.omit(all_d5_mat_scaled[c(1:30),c(1,2,4,3,5)]),
                                column_names_side = "bottom",
                                cluster_rows = T,
                                cluster_columns = F,
                                column_names_rot = 45,
                                column_names_gp = gpar(fontsize=7),
                                row_names_gp = gpar(fontsize=7),
                                row_names_rot = 0,
                                name = "Z-score",width = 6*unit(25,'mm'),
                                column_title = "Control vs Day5-sgRNA")
ht_d5_all
tiff("Plots/Day5_vs_control_wt_d3_readable.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_d5_all)
dev.off()



##saving output of sig genes in csv files
write.csv(d3_sig,"Analysis/res_d3_vs_ct_sig.csv")
write.csv(d5_sig,"Analysis/res_d5_vs_ct_sig.csv")
write.csv(d3_vs_d5_sig,"Analysis/res_d5_vs_d3_sig.csv")







#make list of neighboring genes
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("ensembl_gene_id","start_position","end_position","strand","hgnc_symbol","chromosome_name","ucsc","band")
filters <- c("chromosome_name","start","end")
values <- list(chromosome="18",start=" 15169976",end="16169976")
all.genes_up <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.genes_up
values <- list(chromosome="18",start=" 14165346",end="15165346")
all.genes_d <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.genes_d

all.genes <- rbind(all.genes_d,all.genes_up)


# overlap the neighboring genes to find expression changes

##overlap with d3
d3_sig[rownames(d3_sig) %in% all.genes_d$ensembl_gene_id,]
d3_df[rownames(d3_df) %in% all.genes_d$ensembl_gene_id,]