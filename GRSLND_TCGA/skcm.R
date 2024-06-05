#SKCM and GTEX data 
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
library(tidyverse)
library(recount)
library(data.table)

#load data from tcga-skcm
mrna_query <- GDCquery(project = "TCGA-SKCM",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       experimental.strategy = "RNA-Seq")


GDCdownload(mrna_query, method = "api", files.per.chunk = 100,
            directory = "~/TCGA_Data/TCGA-SKCM/")
mrna_df <- GDCprepare(mrna_query, directory = "~/TCGA_Data//TCGA-SKCM/")
# extract what you need from mrna_df$ to add to metadata
mrna_meta <- mrna_df$sample
mrna_meta <- cbind(mrna_meta, mrna_df$definition)
mrna_meta <- cbind(mrna_meta,  mrna_df$barcode)
mrna_df <- assay(mrna_df)

dim(mrna_df)
mrna_df





#running edgeR and Wilcoxons rank sum test
#edgeR TMM normalize
# conditions <- as.data.frame(t(col_data))[2,]
# conditions <- factor(t(conditions))

grslnd <- count_norm["ENSG00000228203",]
grslnd <- as.data.frame(t(grslnd))
grslnd$type <- y$samples[match(rownames(grslnd), rownames(y$samples)),1]
colnames(grslnd)[1] <- "Count"

ggplot(grslnd, aes(x = type, y = log(Count+1,base = 2), fill = type))+
  geom_violin(trim = FALSE, width=0.5,)+
  theme_cowplot()+
  geom_boxplot(width = 0.15, color = "black", alpha=0.7)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  ggtitle("GAPDH Raw counts")+
  geom_jitter(width = 0.05,alpha=0.4,size=0.2)+
  scale_fill_viridis_d()+
  ylab("Expression - log2(TPM+1") +
  xlab("SKCM
       (num(N)=701;num(T)=471)")

# recount gtex data load
load("~/Documents/AML/rse_gene_skin.Rdata")
rse <- scale_counts(rse_gene)

#pheno data
pheno <- read.table("~/Documents/AML/SRP012682.tsv", sep = "\t",
                    header = T, 
                    stringsAsFactors = F)

## Obtain correct order for pheno data
pheno <- pheno[match(rse$run, pheno$run),]
identical(pheno$Run_s, rse$run)
head(cbind(pheno$run,rse$run))

#read in the gtex data and cleanup
skin_sunexp <- data.table::fread("~/Documents/AML/gene_reads_2017-06-05_v8_skin_sun_exposed_lower_leg.gct")
skin_sunexp$Name <- sub('\\.[0-9]*$','',skin_sunexp$Name)
skin_sunexp$id <- NULL
skin_sunexp$Description <- NULL
skin_sunexp <- skin_sunexp %>% remove_rownames() %>% column_to_rownames(var = "Name")

#remove normal samples from tcga skcm data
`%nin%` = Negate(`%in%`)
mrna_meta <- as.data.table(mrna_meta)
mrna_meta_2 <- mrna_meta[V2 %nin% c("Additional Metastatic","Solid Tissue Normal")]
mrna_counts <- mrna_df[,colnames(mrna_df) %in% mrna_meta_2$V3]
rownames(mrna_counts) <- sub('\\.[0-9]*$','',rownames(mrna_counts))


counts <- assays(rse)$counts
rownames(counts) <- sub('\\.[0-9]*$', '', rownames(counts))
rownames(mrna_df) <- sub('\\.[0-9]*$','',rownames(mrna_df))


merged_counts <- na.omit(merge(mrna_counts, skin_sunexp , by="row.names", all = TRUE))
merged_counts <- merged_counts %>% remove_rownames() %>% column_to_rownames(var = "Row.names")


col_data_skcm <- c(rep("Tumor",471), rep("Normal",701))
names(col_data_skcm) <- colnames(merged_counts)

#make pca plot for raw data without normalization
p <- PCAtools::pca(as.matrix(merged_counts),metadata = as.data.frame(col_data_skcm),removeVar = NULL)

#scree plot
scree_plot <- PCAtools::screeplot(p,axisLabSize = 18,titleLabSize = 22)
tiff("~/Documents/SKCM/scree_plot_raw.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(scree_plot)
dev.off()
#biplot Plots 
pcat <- PCAtools::biplot(p,colby = 'col_data_skcm',legendPosition = 'right',
                         lab=NULL)
tiff("~/Documents/SKCM/PCA_tools_raw.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(pcat)
dev.off()

# run batch correction with ruvseq
col_data_skcm <- as.data.frame(col_data_skcm)
colnames(col_data_skcm) <- c("Type")
##removing batch effects with RUVSeq and DESeq2
#run DESeq to determine control genes
dds <- DESeq2::DESeqDataSetFromMatrix (countData = as.data.frame(merged_counts), colData = col_data_skcm , design = ~ Type )
dds$Type <- relevel(dds$Type,ref = "Normal")
#DE analysis before RUVSeq
dds <- DESeq(dds)
normalized_counts <- counts(dds,normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts)
result_skcm <- results(dds, contrast = c("Type", "Tumor", "Normal")) 
result_skcm <- as.data.frame(result_skcm)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(result_skcm), mart= mart)
result_skcm$GeneName <- gene_IDs[match(rownames(result_skcm), gene_IDs[,1]),2]
result_skcm <- as.data.frame(result_skcm)
result_skcm$log2FoldChange <- as.numeric(result_skcm$log2FoldChange)
normalized_counts$GeneName <- gene_IDs[match(rownames(normalized_counts), gene_IDs[,1]),2]
write.table(result_skcm,file = "~/Documents/SKCM/DEGs_SKCM_GTeX.tsv",quote = F,
            sep = "\t",col.names = T,row.names = T)
write.table(normalized_counts,file = "~/Documents/SKCM/norm_counts_SKCM_GTeX.tsv",quote = F,
            sep = "\t",col.names = T,row.names = T)

#make pca plot for raw data with normalization
p <- PCAtools::pca(mat,metadata = col_data_skcm,removeVar = NULL)

#scree plot
scree_plot <- PCAtools::screeplot(p,axisLabSize = 18,titleLabSize = 22)
tiff("~/Documents/SKCM/scree_plot_norm.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(scree_plot)
dev.off()
#biplot Plots 
pcat <- PCAtools::biplot(p,colby = 'Type',legendPosition = 'right',
                         lab=NULL)
tiff("~/Documents/SKCM/PCA_tools_norm.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(pcat)
dev.off()


# conditions <- factor(t(conditions))
conditions <- factor(col_data_skcm$Type)
y <- DGEList(counts =merged_counts, group = conditions)
## Remove rows conssitently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep, keep.lib.sizes = FALSE]
## Perform TMM normalization and convert to CPM (Counts Per Million)
y <- calcNormFactors(y, method = "TMM")
count_norm <- cpm(y)
count_norm <- as.data.frame(count_norm)
rownames(count_norm) <- gsub("\\..*","",rownames(count_norm))
# count_norm$gene_name <- NULL
# count_norm$gene_name <- gene_IDs[match(rownames(count_norm), gene_IDs[,1]),2]
# #count_norm2 <- count_norm %>% remove_rownames %>%
# #  column_to_rownames(var = "gene_name")
# #count_norm2 <- count_norm2[rownames(count_norm2) %in% non_canonical_list$V1,]
# # Run the Wilcoxon rank-sum test for each gene
# count_norm$gene_name <- NULL
pvalues <- sapply(1:nrow(count_norm), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
  p <- wilcox.test(gene~conditions, data)$p.value
  return(p)})
fdr <- p.adjust(pvalues, method = "fdr")


# Calculate the fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
fdrThres <- 0.05
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(count_norm), mart= mart)
outRst$GeneName <- gene_IDs[match(rownames(outRst), gene_IDs[,1]),2]
outRst <- outRst[order(outRst$FDR),]
outRst_2 <- outRst[is.finite(outRst$log2foldChange),]


vol_plot <- EnhancedVolcano::EnhancedVolcano(result_skcm,
                                 lab = result_skcm$GeneName,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 selectLab = c("GRASLND"),
                                 title = 'SKCM',
                                 subtitle = 'Differential Expressions',
                                 caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                                 pCutoff = 0.05,
                                 FCcutoff = 2,
                                 pointSize = 1.0,
                                 labSize = 3.0,
                                 colAlpha = 1,
                                 legendPosition = 'right',
                                 legendLabSize = 6,
                                 legendIconSize = 1.0,
                                 drawConnectors = TRUE,
                                 widthConnectors = 0.75,
                                 labCol = 'black',
                                 labFace = 'bold',
                                 boxedLabels = TRUE,
                                 colConnectors = 'black')
png("~/Documents/SKCM/vol_plot.png",width = 35,height = 21,units = 'cm',res = 300)
plot(vol_plot)
dev.off()


grslnd <- as.data.frame(mat)["ENSG00000228203",]
grslnd <- count_norm["ENSG00000228203",]
grslnd <- as.data.frame(t(grslnd))
grslnd$type <- col_data_skcm[match(rownames(grslnd),rownames(col_data_skcm)),1]
grslnd$type <- y$samples[match(rownames(grslnd), rownames(y$samples)),1]
colnames(grslnd)[1] <- "Count"
grslnd$Count <- as.numeric(grslnd$Count)

viol_plot <- ggplot(grslnd[-1173,], aes(x = type, y = log(Count+1,base = 2), fill = type))+
  geom_violin(trim = FALSE, width=0.5,)+
  theme_cowplot()+
  geom_boxplot(width = 0.15, color = "black", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  ggtitle("GRASLND DESeq VST counts")+
  geom_jitter(width = 0.05,alpha=0.4,size=0.2)+
  scale_fill_manual(values = c("#a0044d","#610051"))+
  ylab("Expression - log2(VST+1)") +
  xlab("SKCM
       (num(N)=701;num(T)=471)")+
  geom_signif(comparisons = list(c("Normal","Tumor")),
              map_signif_level = T,
              test = wilcox.test,
              y_position = 3.8, tip_length = 0.1, vjust = 0.2)+
  annotate("text",x = 0.58, xend = 1.9, y= 4, label = "Mean LFC = 4.56")+
  annotate("text",x = 0.58, xend = 1.9, y= 3.9, label = "pvalue = 3e-153")
  
tiff("~/Documents/SKCM/viol_plot_vst.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(viol_plot)
dev.off()

stat.test <- grslnd %>%
    (Count ~ type) %>%
  add_significance()
stat.test
fold.change <- log(mean(grslnd[grslnd$type=="Tumor",1])/mean(grslnd[grslnd$type=="Normal",1]),base=2)
fold.change <- log(median(grslnd[grslnd$type=="Tumor",1])/median(grslnd[grslnd$type=="Normal",1]),base=2)

#run gene set enrichemnt and other plots for 






#-----------------------------------------------------------------------------#
#only tcga data

# conditions <- factor(as.data.frame(mrna_meta)[,2])
# y <- DGEList(counts =mrna_df, group = conditions)
# ## Remove rows conssitently have zero or very low counts
# keep <- filterByExpr(y)
# y <- y[keep, keep.lib.sizes = FALSE]
# ## Perform TMM normalization and convert to CPM (Counts Per Million)
# y <- calcNormFactors(y, method = "TMM")
# count_norm <- cpm(y)
# count_norm <- as.data.frame(count_norm)
# rownames(count_norm) <- gsub("\\..*","",rownames(count_norm))
# count_norm$gene_name <- gene_IDs[match(rownames(count_norm), gene_IDs[,1]),2]
# count_norm2 <- count_norm %>% remove_rownames %>%
#   column_to_rownames(var = "gene_name")
# count_norm2 <- count_norm2[rownames(count_norm2) %in% non_canonical_list$V1,]
# # Run the Wilcoxon rank-sum test for each gene
# count_norm$gene_name <- NULL
# pvalues <- sapply(1:nrow(count_norm), function(i){
#   data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
#   p <- wilcox.test(gene~conditions, data)$p.value
#   return(p)})
# fdr <- p.adjust(pvalues, method = "fdr")
# 
# pvalues <- sapply(1:nrow(count_norm), function(i){
#   data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
#   p <- wilcox.test(gene~conditions, data)$p.value
#   return(p)})
# fdr <- p.adjust(pvalues, method = "fdr")
# 
# # Calculate the fold-change for each gene
# conditionsLevel <- levels(conditions)
# dataCon1 <- count_norm[,c(which(conditions==conditionsLevel[1]))]
# dataCon2 <- count_norm[,c(which(conditions==conditionsLevel[2]))]
# foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
# # Output results based on the FDR threshold 0.05
# outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
# rownames(outRst) <- rownames(count_norm)
# outRst <- na.omit(outRst)
# fdrThres <- 0.05
# gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
#                   values = rownames(count_norm), mart= mart)
# outRst$GeneName <- gene_IDs[match(rownames(outRst), gene_IDs[,1]),2]
# 
# 
# EnhancedVolcano::EnhancedVolcano(outRst,
#                                  lab = outRst$GeneName,
#                                  x = 'log2foldChange',
#                                  y = 'FDR',
#                                  title = 'SKCM',
#                                  subtitle = 'Differential Expressions',
#                                  caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
#                                  pCutoff = 0.05,
#                                  FCcutoff = 0.5,
#                                  pointSize = 1.0,
#                                  labSize = 3.0,
#                                  colAlpha = 1,
#                                  legendPosition = 'right',
#                                  legendLabSize = 6,
#                                  legendIconSize = 1.0,
#                                  drawConnectors = TRUE,
#                                  widthConnectors = 0.75,
#                                  labCol = 'black',
#                                  labFace = 'bold',
#                                  boxedLabels = TRUE,
#                                  colConnectors = 'black')
