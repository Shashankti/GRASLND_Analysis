#run gene set enrichemnt and other plots for 

#divide into hot and cold tumorr

#deseq for tcga-laml data
nonZeroCount <- apply(mrna_counts, 1, function(row) {
  return(sum(row > 0))
})
keepInd <- which(nonZeroCount > (length(colnames(mrna_counts)) * .1))
TCGA_df <- mrna_counts[keepInd,]
dds <- DESeq2::DESeqDataSetFromMatrix(TCGA_df,colData = data.frame(sampleID=colnames(TCGA_df)),
                                      design = ~1)
inds <- rownames(TCGA_df)
geoMeansList <- parallel::mclapply(inds,FUN = function(ind){
  row <- TCGA_df[ind,]
  if (all(row == 0)) {
    0
  } else {
    exp(sum(log(row[row != 0]))/length(row))
  }
},mc.cores = 6)
geoMeansList <- unlist(geoMeansList)
dds <- DESeq2::estimateSizeFactors(dds,geoMeans=geoMeansList)
vsd <- DESeq2::vst(dds)
cts <- SummarizedExperiment::assay(vsd)
rownames(cts) <- sub('\\.[0-9]*$', '', rownames(cts))

#make high and low expression of gapdh data
imp <- cts["ENSG00000228203",]
imp <- as.data.frame(imp)
colnames(imp) <- c("GRASLND")
high <- rownames(imp)[which(imp$GRASLND>=median(imp$GRASLND))]
low <- rownames(imp)[which(imp$GRASLND<median(imp$GRASLND))]
dim(high)
high_gapdh <- cts[,high]
low_gapdh <- cts[,low]
#make plot for high vs low
high_low <- imp
high_low$strata <- ifelse(imp$GRASLND>median(imp$GRASLND),"high","low") 

attach(loadNamespace("enrichplot"),name = "enrichplot_all")


#make high vs low plot

expr_plot <- ggplot(high_low, aes(x = strata, y = GRASLND, fill = strata))+
  geom_violin(trim = FALSE, width=0.4)+
  theme_cowplot()+
  geom_boxplot(width = 0.1, color = "grey", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  geom_jitter(height = 0, width = 0.1,alpha=0.7)+
  set_enrichplot_color(type = "color")+
  scale_x_discrete(labels=c("High\n(n=236)","Low\n(n=235)"))+
  ylab("Normalized Expression") +
  xlab("")+
  geom_signif(comparisons = list(c("high","low")),
              map_signif_level = T,
              test = wilcox.test,
              y_position = 11.9, tip_length = 0.1, vjust = 0.2)
# annotate("text",x = 0.58, xend = 1.9, y= 17, label = "Mean Log2FC = 0.1")+
# annotate("text",x = 0.58, xend = 1.9, y= 16.8, label = "pvalue = 6.80e-26")
tiff("../SKCM/high_vs_low.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(expr_plot)
dev.off()

#split data based on cold and hot tumor
#$calculate top 10 and bottom 10% percentage of CD8A
imp <- cts["ENSG00000153563",]
imp <- as.data.frame(imp)
colnames(imp) <- c("CD8A")
high <- imp %>% filter(quantile(CD8A,0.9)<CD8A) %>% rownames()
low <- imp %>% filter(quantile(CD8A,0.1)>CD8A) %>% rownames()
cd8a_counts<- mrna_counts[,c(high,low)]
cd8a_meta <- data.frame("Type" = c(rep("Hot",47),rep("Cold",47)))
rownames(cd8a_meta) <- c(high,low)

dds <- DESeqDataSetFromMatrix(cd8a_counts,colData = cd8a_meta,
                              design = ~Type)
keep_2f <- rowSums(counts(dds)) >= 10
dds<- dds[keep_2f,]
dds

#get normalised counts information
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dds$sample <- relevel(dds$Type,ref="Cold")

# Step 3: Run DESeq ----------------------

dds <- DESeq(dds)
resultnames <- resultsNames(dds)
resultnames

res_hot_vs_cold <- results(dds,contrast = c("Type","Hot","Cold"))
res_hot_vs_cold <- as.data.frame(res_hot_vs_cold)
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(res_hot_vs_cold), mart= mart)

res_hot_vs_cold$GeneName <- gene_IDs[match(rownames(res_hot_vs_cold),gene_IDs[,1]),2]
res_hot_vs_cold <- arrange(res_hot_vs_cold,padj)
#read in list of lncRNAs
lncRNA <- fread("lncRNA.txt",header = F)
lncRNA

hot_col_vol <- EnhancedVolcano::EnhancedVolcano(res_hot_vs_cold[res_hot_vs_cold$GeneName %in% lncRNA$V1,],
                                 lab = res_hot_vs_cold[res_hot_vs_cold$GeneName%in% lncRNA$V1,]$GeneName,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 selectLab = c("GRASLND","CD8A","IFNG","HLA-A"),
                                 title = 'Hot vs Cold SKCM',
                                 subtitle = 'Top and bottom 10% CD8A expression samples(47+47)',
                                 caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                                 pCutoff = 0.05,
                                 FCcutoff = 1,
                                 pointSize = 1.0,
                                 labSize = 3.5,
                                 colAlpha = 1,
                                 legendPosition = 'right',
                                 legendLabSize = 7,
                                 legendIconSize = 1.0,
                                 drawConnectors = TRUE,
                                 widthConnectors = 0.85,
                                 lengthConnectors = unit(0.01,'npc'),
                                 labCol = 'black',
                                 labFace = 'bold',
                                 boxedLabels = TRUE,
                                 colConnectors = 'black')

tiff("../SKCM/hot_vs_cold_differential_expression.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(hot_col_vol)
dev.off()

# correlation plots
#correlation between TCGA SKCM vs genes from LIMIT paper



#import the counts file

library(data.table)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(cowplot)
library(DGEobj.utils)
library(ggpubr)

#get gene lenghts
library(GenomicFeatures)
txdb <- GenomicFeatures::makeTxDbFromGFF("~/CLIP_Seq/Data/Genome/Homo_sapiens.GRCh38.108.gtf",format = "gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))

df <- as.data.frame(exonic.gene.sizes[rownames(mrna_counts),])
rownames(df) <- rownames(mrna_counts)
#tpm counts 
#TCGA LAML raw counts
TCGA_df
#generate TPM normalized counts
nm_tcga <- convertCounts(mrna_counts,
                         unit = "TPM",
                         geneLength = as.matrix(df),
                         log = FALSE
)
nm_tcga <- as.data.frame(nm_tcga)
nm_tcga$GeneName <-  gene_IDs[match(rownames(nm_tcga), gene_IDs[,1]),2]

imp <- nm_tcga[nm_tcga$GeneName %in% c("HLA-A", "GRASLND"),]
imp <- as.data.frame(t(imp))
colnames(imp) <- c("IFNG", "GRASLND")
imp <- head(imp,-1)
imp[,1] <- as.numeric(imp[,1])
imp[,2] <- as.numeric(imp[,2])
#calculate correlation
cor.test(imp$IFNG,imp$GRASLND,method = "spearman",exact = FALSE)



p1 <- ggplot(imp, aes(x=GRASLND,y=IFNG))+
  geom_point(size=2,alpha=0.5)+
  theme_cowplot()+
  ylab("MHC-I Normalized Counts")+
  ggtitle("GRASLND vs. MHC-I")+
  theme(plot.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t(30))))+
  xlab(expression(paste("         GRASLND Normalized counts \n Spearmann p=0.03 coeff=-0.09(n=471)")))
# coord_cartesian(ylim = c(0,30000), clip = "off")
tiff("../SKCM/GRASLND_vs_MHC1.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(p1)
dev.off()





#scatter plot
sp <- ggscatter(imp,x="GAPDH",y="RPL13A",
                add = "reg.line",
                add.params = list(color="blue",fill="lightgray"),
                conf.int = TRUE,
                size=2,alpha=0.5)+
  theme_cowplot()+
  ylab("RPL13A Normalized Counts")+
  ggtitle("GAPDH vs. RPL13A")+
  theme(plot.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t(30))))

sp + stat_cor(method = "spearman")

ggplot(imp, aes(x = ,y = log10(Normalized), fill = Type))+
  geom_violin(trim = FALSE, width=0.4)+
  theme_cowplot()+
  geom_boxplot(width = 0.1, color = "grey", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  ggtitle("GAPDH log-normalized counts")+
  geom_jitter(height = 0, width = 0.1,alpha=0.7)+
  scale_fill_viridis_d()+
  ylab("log10-normalized Counts") +
  xlab("Type")




ggscatter(imp, x="MelanA", y="GRSLND")

theme(axis.text.x=element_blank(),
      axis.ticks.x=c(0,50,100,150,200),
      axis.text.y = element_blank(),
      axis.ticks.y=c(0,5000,15000,25000))
#--------------------------------------------------------------------------------------------_###

#GSEA plots
corrMat <- generate_corr(cts,cores = 6,gene = "ENSG00000228203")
corrMat$Spearman <- round(as.numeric(corrMat$Spearman),digits = 9)
corrMat <- corrMat[order(corrMat$Spearman,decreasing = TRUE),]
Spearman <- corrMat$Spearman
names(Spearman) <- corrMat$GeneID

#run GSEA analysis
##Load hallmark data
#h_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name,human_ensembl_gene)
#run fgsea
gsea <- clusterProfiler::GSEA(Spearman,
                              minGSSize = 15,
                              maxGSSize = 500,
                              pAdjustMethod = "BH",
                              seed = TRUE,
                              eps = 1e-100,
                              nPermSimple=1000,
                              by = "fgsea",
                              TERM2GENE = h_gene_set)
gsea_result <- as.data.frame(gsea)
gsea_result_summary <- gsea_result[,c(1,3,5:8)]
gsea_result_summary$adjPvalue <- ifelse(gsea_result_summary$p.adjust <= 0.05, "significant", "non-significant")
cols <- c("non-significant" ="grey","significant"="red")

#general Plot
gsea_hm <- ggplot(gsea_result_summary, aes(reorder(ID, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")+
  theme_cowplot()
tiff("../SKCM/hallmark_pathway.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_hm)
dev.off()

#dotplot for upregulated
##make sorted dataframe
ggdata_up <- gsea_result_summary[gsea_result_summary$NES>0,]
ggdata_up$ID <- reorder(ggdata_up$ID,ggdata_up$NES,decreasing = F)


gsea_dotplot_up <- ggplot(ggdata_up[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
  geom_point(shape=21)+
  scale_size(range = c(3,8))+
  scale_color_continuous(low='red',high='blue')+
  scale_fill_continuous(low='red',high='blue')+
  xlab('Gene set')+
  ylab('Enrichment Score')+
  labs(title = "Summary dotplot (upregulated)")+
  theme_cowplot()+
  theme(
    legend.position = 'right',
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
    axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
    axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
    axis.title = element_text(size = 11, face = 'bold'),
    axis.title.y = element_text(size = 11, face = 'bold'),
    legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 12, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
  guides(size = guide_legend(order = 1))+
  coord_flip()
tiff("../SKCM/hallmark_gsea_up_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dotplot_up)
dev.off()

#dotplot for downregulated
##make sorted dataframe
ggdata_down <- gsea_result_summary[gsea_result_summary$NES<0,]
ggdata_down$ID <- reorder(ggdata_down$ID,ggdata_down$NES,decreasing = F)


gsea_dotplot_down <- ggplot(ggdata_down[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
  geom_point(shape=21)+
  scale_size(range = c(3,8))+
  scale_color_continuous(low='red',high='blue')+
  scale_fill_continuous(low='red',high='blue')+
  xlab('Gene set')+
  ylab('Enrichment Score')+
  labs(title = "Summary dotplot (downregulated)")+
  theme_cowplot()+
  theme(
    legend.position = 'right',
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
    axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
    axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
    axis.title = element_text(size = 11, face = 'bold'),
    axis.title.y = element_text(size = 11, face = 'bold'),
    legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 12, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
  guides(size = guide_legend(order = 1))+
  coord_flip()
tiff("../SKCM/hallmark_gsea_down_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dotplot_down)
dev.off()



glyco_gsea <- enrichplot::gseaplot2(gsea,geneSetID = gsea@result$ID[35],title = gsea@result$Description[35])
tiff("../SKCM/EPITHELIAL_MESENCHYMAL_TRANSITION.tiff",width = 35,height = 21,units = 'cm',res = 300)
glyco_gsea
dev.off()

#try gene ontlogy databsae instead of hallmark database
##run geneset enrichment analysis for database from GeneOntology
gsea_go <- clusterProfiler::gseGO(Spearman,
                                  OrgDb = "org.Hs.eg.db",
                                  ont = "BP",
                                  eps = 1e-300,
                                  by = "fgsea",
                                  keyType = "ENSEMBL")
gsea_result <- as.data.frame(gsea_go)
gsea_result_summary <- gsea_result[,c(1,3,5:8)]
gsea_result_summary$adjPvalue <- ifelse(gsea_result_summary$p.adjust <= 0.05, "significant", "non-significant")
cols <- c("non-significant" ="grey","significant"="red")

#make plot with GO data
gsea_hm <- ggplot(gsea_result_summary, aes(reorder(ID, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO pathways Enrichment Score from GSEA")+
  theme_cowplot()
tiff("../SKCM/gene_ontology_pathway.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_hm)
dev.off()

enrichplot::dotplot(gsea_go,showCategory=20,split=".sign",font.size=7,
        label_format=60)+facet_grid(.~.sign) + ggtitle("")

goplot(gsea_go)
gsea_go_df <- as.data.frame(gsea_go)
selected_pathways <- head(gsea_go_df[gsea_go_df$enrichmentScore>0,2],n=25)
dotplot <- enrichplot::dotplot(gsea_go,showCategory=selected_pathways,split=,font.size=9,
                               label_format=60)
tiff("../SKCM/gene_ontology_upregulated.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(dotplot)
dev.off()

#individual gsea GO plots
glyco_gsea <- enrichplot::gseaplot2(gsea_go,geneSetID = gsea_go@result$ID[48],title = gsea_go@result$Description[48])
tiff("../SKCM/lymphocyte_mediated_immunity.tiff",width = 35,height = 21,units = 'cm',res = 300)
glyco_gsea
dev.off()

#splitting data into high and low 



