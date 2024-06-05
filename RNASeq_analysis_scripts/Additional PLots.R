#Comparision between DEGs from control vs shRNa and  Ifn vs Dox
library(data.table)
library(ggplot2)
library(viridis)
library(cowplot)
library(ChIPpeakAnno)
library(VennDiagram)
library(RColorBrewer)
library(GenomicRanges)
myCol <- brewer.pal(3, "Pastel2")


KO_genes <- rownames(KO_sig)
Ifn_genes <- rownames(KOTreatment_sig)
Off_genes <- rownames(cTreatment_sig)

which(KO_genes %in% Ifn_genes)
which(Off_genes %in% KO_genes)
#making venn diagrams
venn.diagram(x = list(KO_genes,Ifn_genes,Off_genes),
             category.names = c("Control vs Knockdown","shRNA Ifn vs Dox","Control Ifn vs Dox"),
             filename = "DEG_overlap.png",
             output=T,
             # Output features
             imagetype="png" ,
             height = 1200, 
             width = 1280,  
             resolution = 300,
             compression = "lzw",
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)

#only overlap between exp 1 and 2

#take gene ids from experiment 1
exp1_genes <- rownames(KO_sig)
exp2_genes <- rownames(KOTreatment_sig)


## adding negative control with random genes 
n <- dim(KOTreatment_df)[1]
sampl <- sample.int(n, 500)

df_sampld <- KOTreatment_df[sampl,]
df_sampld$log2FoldChange <- as.numeric(df_sampld$log2FoldChange)
df_sampld[rownames(df_sampld) %nin%,]
logfoldChanges <- c(KOTreatment_df$log2FoldChange,
                    KOTreatment_df[rownames(KO_sig),2],
                    KOTreatment_df[rownames(KOTreatment_sig),2],
                    df_sampld$log2FoldChange)
types <- c(rep("All", dim(KOTreatment_df)[1]),
           rep("OnlyKD",dim(KO_sig)[1]),
           rep("OnlyDox",dim(KOTreatment_sig)[1]),
           rep("NCT",dim(df_sampld)[1]))
comb_df <- data.frame(types, logfoldChanges)
comb_df_summary <- ddply(comb_df, .(types), summarise, median = median(logfoldChanges))

pvalue_list = NULL
pvalue_list[1] <- 1
pvalue_list[4] <- wilcox.test(comb_df[comb_df$types=="All",2], comb_df[comb_df$types=="OnlyKD",2])$p.value

comb_df_summary <- cbind(comb_df_summary, pvalue_list)



ggplot(comb_df,aes(logfoldChanges,colour = types)) +
  stat_ecdf()+
  scale_color_hue(name="Types",labels = c('All','NCT','OnlyDox','OnlyKD')) + 
  theme_cowplot()+
  ggtitle("CDF for all genes")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  geom_label_repel(data=comb_df_summary, aes(x = median, label=median, y= 0.01),
                   size = 3, segment.size=0.25, nudge_x = 0.5, direction = "y",
                   hjust=0) +
  geom_point(data = comb_df_summary,aes(x = median, y= 0.01, shape = types), size=2)
  # annotate("text", x= -1.8, y= 1, label = "21,504 all genes",color = "red")+
  # annotate("text",x=-1.8,y=0.97,label = "(FC=0.005864923)", color = "red") +
  # annotate("text",x=-1.8,y=0.90, label = "324 Gal Genes",color = "#7dad02") + 
  # annotate("text",x=-1.8,y=0.87, label = "(FC=-0.2044665,", color = "#7dad02") + 
  # annotate("text",x=-1.8,y=0.84, label = "P = 6.44e-04)", color = "#7dad02") + 
  # annotate("text",x=-1.8,y=0.78, label = "449 Low Genes", color = "#c57ef8") +
  # annotate("text",x=-1.8,y=0.75,label = "(FC=-0.1975091,", color = "#c57ef8") +
  # annotate("text",x=-1.8,y=0.72,label = "P=9.35e-06)",color = "#c57ef8") +
  # annotate("text",x=-1.8,y=0.67,label= "375 Norm Genes", color = "#00bfc4") +
  # annotate("text",x=-1.8,y=0.64,label = "(FC=-0.2585118,", color = "#00bfc4") + 
  # annotate("text",x=-1.8,y=0.61,label = "P=5.043-06", color = "#00bfc4")

#make  the violin plot for exp1 
#take normalized counts for 
exp1_counts <- KO_sig$log2FoldChange
exp2_counts <-  KOTreatment_df[rownames(KO_sig),2]

violin_df <- data.frame(FC = append(exp1_counts, exp2_counts), Type = c(rep("OnlyKD",549 ), rep("KDIfn", 549)))
FC_violin <- ggplot(violin_df, aes(x = Type, y=FC, fill=Type))+
  geom_violin(trim = F, width=0.4)+
  theme_cowplot()+
  geom_boxplot(width=0.1, color="grey",alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size=11))+
  ggtitle("Log Fold change distribution")+
  geom_jitter(height = 0, width = 0.1, alpha=0.3)+
  scale_fill_viridis_d()+
  ylab("LogFoldChange")+
  xlab("Experiment")

#violin plot for raw counts
lacz_counts <- rowMeans(counts_data[Ifn_gamma_genes$human_ensembl_gene,1:3])
KD_counts <- rowMeans(counts_data[Ifn_gamma_genes$human_ensembl_gene,4:7])
lacz_ifn_counts <- rowMeans(counts_data[Ifn_gamma_genes$human_ensembl_gene,8:9])
lacz_ifn_dox <- rowMeans(counts_data[Ifn_gamma_genes$human_ensembl_gene,10:11])
shRNA_ifn_counts <- rowMeans(counts_data[Ifn_gamma_genes$human_ensembl_gene,c(12,13,16,17)])
shRNA_ifn_dox <- rowMeans(counts_data[Ifn_gamma_genes$human_ensembl_gene,c(14,15,18,19)])
violin_df_2 <- na.omit(data.frame(lacz_counts, KD_counts, lacz_ifn_counts,lacz_ifn_dox, shRNA_ifn_counts, shRNA_ifn_dox))
violin_df_2[,5:6] %>%
  gather(key = "MesureType", value = "val") %>%
  ggplot(aes(x=MesureType, y=log10(val), fill=MesureType))+
  geom_violin(width=0.4, trim = F)+
  theme_cowplot()+
  geom_boxplot(width=0.1,color="grey", alpha=0.5)+
  theme(plot.title = element_text(size = 11))+
  ggtitle("Raw Counts for KD Significant genes")+
  geom_jitter(height = 0, width=0.01, alpha=0.2)+
  scale_fill_viridis_d()+
  ylab("Raw counts")+
  xlab("Type")+
  scale_x_discrete(labels=c("shRNAIfn","shRNADox"))
#violin plot for normalized counts
lacz_counts <- rowMeans(normalized_counts_KO[rownames(KO_sig), 1:3])
KO_counts <- rowMeans(normalized_counts_KO[rownames(KO_sig),4:7])
lacz_ifn_counts <- rowMeans(normalized_counts_tr[rownames(KO_sig),1:2])
lacz_ifn_dox <- rowMeans(normalized_counts_tr[rownames(KO_sig),3:4])
shRNA_ifn_counts <- rowMeans(normalized_counts_tr[rownames(KO_sig),c(5,6,9,10)])
shRNA_ifn_dox <- rowMeans(normalized_counts_tr[rownames(KO_sig),c(7,8,11,12)])
violin_df_3 <- data.frame(lacz_counts, KO_counts, lacz_ifn_counts, lacz_ifn_dox, shRNA_ifn_counts, shRNA_ifn_dox)
violin_df_3 %>%
  gather(key = "MesureType", value = "val") %>%
  ggplot(aes(x=MesureType,y=log10(val), fill=MesureType))+
  geom_violin(width=0.4, trim = F)+
  theme_cowplot()+
  geom_boxplot(width=0.1,color="grey", alpha=0.5)+
  theme(plot.title = element_text(size = 11))+
  ggtitle("Normalized counts for KD significant genes")+
  geom_jitter(height = 0, width=0.01, alpha=0.2)+
  scale_fill_viridis_d()+
  ylab("Norm Counts")+
  xlab("Type")



#neighbor gene analysis

#make list of neighboring genes

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("ensembl_gene_id","start_position","end_position","strand","hgnc_symbol","chromosome_name","ucsc","band")
filters <- c("chromosome_name","start_position","end_position")
values <- list(chromosome_name="2",start_position="5911754",end_position="6918734")
all.genes_up <- biomaRt::getBM(attributes=attributes, filters=filters,  mart=mart,values = values)
getBM(attributes = c('affy_hg_u133_plus_2','ensembl_gene_id'), 
      filters = c('chromosome_name','start','end'),
      values = list(16,1100000,1250000), 
      mart = mart)
values <- list(chromosome="18",start=" 14165346",end="15165346")
all.genes_d <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.genes_d



all.genes <- data.table::fread("close_genes.txt",header = F)

results=getBM(attributes = c("hgnc_symbol","entrezgene", "chromosome_name",
                             "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region"),values = values, 
              mart = mart)
# overlap the neighboring genes to find expression changes

##overlap with KO_df
counts_data[all.genes$V1,]
#neighbor genes heatmap 
all_nbg <- counts_data[all.genes$V1,]
all_nbg <- na.omit(all_nbg)
all_nbg$GeneID <- rownames(all_nbg)
all_nbg$GeneID <- sub("\\.\\d+$", "", all_nbg$GeneID)
all_nbg <- all_nbg[unique(all_nbg$GeneID),]
keep_2f <- rowSums(all_nbg[,1:19]) > 0
all_nbg <- all_nbg[keep_2f,]

# all_nbg$padj <- d5_sig[match(rownames(all_d5),rownames(d5_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
#only_d5_mat <- as.matrix(only_low[,c(2,3,4,5)])
all_nbg_mat_scaled <- t(scale(t(all_nbg[,c(1:19)])))
#rownames(all_d5_mat_scaled) <- all_d5$GeneName

all_d5_mat_scaled <- cbind(all_d5_mat_scaled,all_d5[match(rownames(all_d5_mat_scaled),all_d5$GeneName),8])
ht_nbg_all <-ComplexHeatmap::Heatmap(all_nbg_mat_scaled,
                                     column_names_side = "bottom",
                                     cluster_rows = T,
                                     cluster_columns = F,
                                     column_names_rot = 45,
                                     column_names_gp = gpar(fontsize=7),
                                     row_names_gp = gpar(fontsize=7),
                                     row_names_rot = 0,
                                     name = "Z-score",width = 6*unit(25,'mm'),
                                     column_title = "Raw Neighbor Gene Expression")
ht_nbg_all
tiff("Plots/raw_neighbor_gene_exp.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_nbg_all)
dev.off()


