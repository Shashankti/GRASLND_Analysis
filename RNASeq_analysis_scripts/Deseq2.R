# DESeq2 analysis from the ZARP counts dta

library (DESeq2)
library (tidyverse)
library (airway)
library(biomaRt)

# Step 1: preparing count data ----------------
# read in counts data
counts_data <- read.table ('Data/genes_numreads.tsv',header = T)
head(counts_data)
counts_data <- counts_data %>% column_to_rownames("Name")

# read in sample info
colData_2f<-read.table ('Data/samples.tsv')
#view (colData_2f)



#making sure the row names in colData matches the column names in counts_data 

all(colnames(counts_data) %in% rownames(colData_2f))

#check the order
all(colnames(counts_data) == rownames(colData_2f))

#do  lacZ vs shRNA comparision
sample_info_1 <- colData_2f[1:7,]
counts_data[,1:7]

# Step 2: construct a DESeqDataSet object
dds_KO <- DESeqDataSetFromMatrix (countData = round(counts_data[,1:7]), colData = colData_2f[1:7,] , design = ~ sample )

keep_2f <- rowSums(counts(dds_KO)) >= 10
dds_KO <- dds_KO[keep_2f,]
dds_KO

#get normalised counts information
dds_KO <- estimateSizeFactors(dds_KO)
sizeFactors(dds_KO)
dds_KO$sample <- relevel(dds_KO$sample,ref="Control")

# Step 3: Run DESeq ----------------------

dds_KO <- DESeq(dds_KO)
resultnames_KO <- resultsNames(dds_KO)
resultnames_KO

res_KO_vs_Control <- results(dds_KO,contrast = c("sample","shRNA","Control"))
normalized_counts_KO <- counts(dds_KO,normalized=TRUE)
normalized_counts_KO <- as.data.frame(normalized_counts_KO)

#comparision between lacZ Ifn vs lacZ ifnDox
dds_treatment <- DESeqDataSetFromMatrix (countData = round(counts_data[,8:19]), colData = colData_2f[8:19,] , design = ~ sample + condition + sample:condition)

keep_2f <- rowSums(counts(dds_treatment)) >= 10
dds_treatment <- dds_treatment[keep_2f,]
dds_treatment

#get normalised counts information
dds_treatment <- estimateSizeFactors(dds_treatment)
sizeFactors(dds_treatment)
#dds_treatment$sample <- relevel(dds_KO$sample,ref="Control")

# Step 3: Run DESeq ----------------------

dds_treatment <- DESeq(dds_treatment)
resultnames_treatment <- resultsNames(dds_treatment)
resultnames_treatment
mod_mat_2f <- model.matrix(design(dds_treatment), colData(dds_treatment))
mod_mat_2f

Ifn_control <- colMeans(mod_mat_2f[dds_treatment$condition == "Ifn" & dds_treatment$sample == "Control", ])
Dox_control <- colMeans(mod_mat_2f[dds_treatment$condition == "Dox" & dds_treatment$sample == "Control", ])
Ifn_shRNA <- colMeans(mod_mat_2f[dds_treatment$condition == "Ifn" & dds_treatment$sample == "shRNA", ])
Dox_shRNA <- colMeans(mod_mat_2f[dds_treatment$condition == "Dox" & dds_treatment$sample == "shRNA", ])

res_cIfn_vs_cDox <-results (dds_treatment, contrast = Dox_control - Ifn_control) 
res_KOIfn_vs_KODox <-results (dds_treatment, contrast = Dox_shRNA - Ifn_shRNA) 


normalized_counts_tr <- counts(dds_treatment,normalized=TRUE)
normalized_counts_tr <- as.data.frame(normalized_counts_tr)



#merge technical replicates
write.table(normalized_counts,file = "./Analysis/Normalized_counts_merged.txt",sep = "\t",quote = F,col.names = NA)



#summary of results
summary(res_KO_vs_Control)
summary(res_cIfn_vs_cDox)
summary(res_KOIfn_vs_KODox)

#convert to dataframe
KO_df <- as.data.frame(res_KO_vs_Control)
cTreatment_df <- as.data.frame(res_cIfn_vs_cDox)
KOTreatment_df <- as.data.frame(res_KOIfn_vs_KODox)

# LFC shrinkage for visualisation

#annotate geneIDs

## control vs d3
mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart = "ensembl",host = "https://useast.ensembl.org"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(KO_df), mart= mart)

KO_df$GeneName <- gene_IDs[match(rownames(KO_df),gene_IDs[,1]),2]
KO_df <- arrange(KO_df,padj)

gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(normalized_counts_KO), mart= mart)
normalized_counts_KO$GeneName <- gene_IDs[match(rownames(normalized_counts_KO),gene_IDs[,1]),2]



## control ifn vs dox
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(cTreatment_df), mart= mart)

cTreatment_df$GeneName <- gene_IDs[match(rownames(cTreatment_df),gene_IDs[,1]),2]
cTreatment_df <- arrange(cTreatment_df,padj)

## shRNA Ifn vs Dox
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(KOTreatment_df), mart= mart)

KOTreatment_df$GeneName <- gene_IDs[match(rownames(KOTreatment_df),gene_IDs[,1]),2]
KOTreatment_df <- arrange(KOTreatment_df,padj)


gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(normalized_counts_tr), mart= mart)
normalized_counts_tr$GeneName <- gene_IDs[match(rownames(normalized_counts_tr),gene_IDs[,1]),2]


#write the output results
write.table(KO_df,file = "Analysis/Control_vs_KO_res.tsv",sep = "\t",quote = F,col.names = NA)
write.table(cTreatment_df,file = "Analysis/Control_Ifv_vs_Dox_res.tsv",sep = "\t",quote = F,col.names = NA)
write.table(KOTreatment_df,file = "Analysis/KO_Ifn_vs_Dox_res.tsv",sep = "\t",quote = F,col.names = NA)
#norm counts
write.table(normalized_counts_KO,file = "Analysis/Norm_counts_KO.tsv",sep = "\t",quote = F,col.names = NA)
write.table(normalized_counts_tr,file = "Analysis/Nomr_counts_tr.tsv",sep = "\t",quote = F,col.names = NA)


# Add gene names to normalised counts
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = normalized_counts[,1], mart= mart)

#!normalized_counts <- merge(normalized_counts,gene_IDs, by.x="GeneID",by.y="ensembl_gene_id")

#addiditional QC steps
fig_countdist <- as.data.frame(counts_data) %>% 
  pivot_longer(everything(), names_to = "library") %>%
  mutate(value = value+1) %>%
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 50)+
  scale_x_log10()+
  facet_wrap(~ library, nrow =3)
print(fig_countdist)

x_log10 <- log10(sweep(counts_data, 2, 1e6 / colSums(counts_data), '*') + 1)
x_log10_highexp <- x_log10[rowSums(x_log10 > 0) > 3, ]

# only use high expressed genes
x_pca <- prcomp(t(x_log10_highexp), scale = FALSE)


# plot(x_pca$x[, c('PC1', 'PC2')])
# text(x_pca$x[, c('PC1', 'PC2')], labels = rownames(x_pca$x))

fig_pca <- as.data.frame(x_pca$x) %>%
  mutate(libname = rownames(x_pca$x)) %>%
  ggplot(aes(x = PC1, y = PC2, label = libname)) +
  geom_label()
print(fig_pca)

