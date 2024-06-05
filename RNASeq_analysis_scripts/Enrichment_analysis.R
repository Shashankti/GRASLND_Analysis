# GO and GSEA for RNA-Seq data

library(clusterProfiler)
library(msigdbr)
library(goseq)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggridges)
library(pathview)
library(ggplot2)
library(magrittr)
library(path)
set.seed(2022)


#define functions
`%nin%` = Negate(`%in%`)


#only for control vs shRNA
#Select up regulated gene

Enrichment_analysis <- function(deg_table){
  bg_genes <- rownames(deg_table)
  gene_pval <- deg_table$pvalue
  gene_qval <- deg_table$qvalue
  gene_fc <- deg_table$log2FoldChange
  names(gene_pval) <- names(gene_qval) <- names(bg_genes)
  up_genes <- rownames(deg_table[deg_table$threshold==TRUE & deg_table$log2FoldChange>0,])
  down_genes <- rownames(deg_table[deg_table$threshold==TRUE & deg_table$log2FoldChange<0,])
  sig_genes <- rownames(deg_table[deg_table$threshold==TRUE,])
  #make the GO res object
  ego_all <- enrichGO(gene = sig_genes,
                      universe = bg_genes,
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENSEMBL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1.0,
                      qvalueCutoff = 1.0,
                      readable = F)
  ego_up <- enrichGO(gene = up_genes,
                      universe = bg_genes,
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENSEMBL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1.0,
                      qvalueCutoff = 1.0,
                      readable = F)
  ego_down <- enrichGO(gene = down_genes,
                      universe = bg_genes,
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENSEMBL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1.0,
                      qvalueCutoff = 1.0,
                      readable = F)
  #prepare for gsea input
  deg_table_GSEA <- deg_table[deg_table$baseMean>50,]
  deg_table_GSEA <- deg_table_GSEA[order(-deg_table_GSEA$stat),]
  gene_list_gsea <- deg_table_GSEA$stat
  names(gene_list_gsea) <- rownames(deg_table_GSEA)
  gse <- gseGO(gene_list_gsea,
               ont = "BP",
               OrgDb = "org.Hs.eg.db",
               keyType = "ENSEMBL",
               eps = 1e-300)
  #make the output plots
  up_go_bar <- barplot(ego_up, showCategory = 25, font.size = 7)
  down_go_bar <- barplot(ego_down, showCategory = 25, font.size = 7)
  all_bar <- barplot(ego_all, showCategory = 25, font.size = 7)
  up_go_dot <- dotplot(ego_up, showCategory = 25, font.size = 7)
  down_go_dot <- dotplot(ego_down, showCategory = 25, font.size =7)
  all_go_dot <- dotplot(ego_all, showCategory = 25, font.size =7)
  #prepare for cnet plot
  edox_up <- setReadable(ego_up, "org.Hs.eg.db", "ENSEMBL")
  edox_down <- setReadable(ego_down, "org.Hs.eg.db", "ENSEMBL")
  cnet_up <- cnetplot(edox_up, showCategory=15, color.params = list(foldChange = gene_fc[gene_fc>0],category="#9370DB"),
                      cex.params = list(category_label=0.9, gene_label=0.5, gene_node=0.3))
  cnet_down <- cnetplot(edox_down, showCategory=15, color.params = list(foldChange = gene_fc[gene_fc<0],category="#9370DB"),
                        cex.params = list(category_label=0.9, gene_label=0.5, gene_node=0.3))
  heat_up <- heatplot(edox_up, showCategory = 10, foldChange = gene_fc[gene_fc>0],
                      pvalue = gene_qval[names(gene_qval) %in% edox_up@gene])
  heat_down <- heatplot(edox_down, showCategory = 10, foldChange = gene_fc[gene_fc>0],
                      pvalue = gene_qval[names(gene_qval) %in% edox_down@gene])
  # cnet_up <- cnetplot(edox_up, showCategory = 15, foldChange = gene_fc[gene_fc>0], categorySize='p.adjust')
  # cnet_down <- cnetplot(edox_down, showCategory = 15, foldChange = gene_fc[gene_fc>0], categorySize='p.adjust')
  #make gsea plot
  gsea_dot <- dotplot(gse,showCategory=20,split=".sign",font.size=7) +
    facet_grid(.~.sign)
  ego_down <- pairwise_termsim(ego_down)
  ego_up <- pairwise_termsim(ego_up)
  gse <- pairwise_termsim(gse)
  return(list(up_go_bar, down_go_bar,all_bar,
         up_go_dot, down_go_dot, all_go_dot,
         cnet_up, cnet_down, gsea_dot,
         heat_up, heat_down,
         gse, ego_up, ego_down))
}

Plots <- Enrichment_analysis(KO_df)
#save output of GO as excel file
write.table(Plots[[]])
#save the output of plots
tiff("Plots/new_GObar_control_vs_shRna_up.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(Plots[[1]] + ggtitle("Up regulated Control vs shRNA")+theme(plot.title = element_text(size=15)))
dev.off()

#write the results from geneontology
wrti

#make emappplot
tiff("Plots/emap_Control_vs_shRNA_down.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(emapplot(Plots[[14]],cex_line=0.5,cex_label_category=0.7, cex_category=0.7, showCategory=15) +
       ggtitle("Down regulated Control vs shRNA")+theme(plot.title = element_text(size=15)))
dev.off()


up_genes <- rownames(KO_sig[KO_sig$log2FoldChange>0,])
#make the over representation test
GO_res <- enrichGO(gene = up_genes,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   ont = "BP")
as.data.frame(GO_res)
# plot the results
fit1 <- plot(barplot(GO_res,showCategory = 30,font.size = 7,title = "Up-regulated genes control vs shRNA"))
tiff("Plots/GO_bar_up_control_vs_shRNA.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(fit1)
dev.off()

# same process for the downregulated genes
down_genes <- rownames(KO_df_shrunk[KO_df_shrunk$log2FoldChange<0,])

GO_res_d <- enrichGO(gene = down_genes,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENSEMBL",
                     ont = "BP")
as.data.frame(GO_res_d)
fit2 <- plot(barplot(GO_res_d,showCategory = 30,font.size = 7,title = "Significantly down-regulated genes in control vs Day5"))
fit2

# make a dot plot
dot_go <- dotplot(GO_res,showCategory=20,title="Up-regulated GO-enrichment for control vs KO",font.size=7)
dot_go_d <- dotplot(GO_res_d,showCategory=20,title="Down-regulated GO-enrichment for controld vs KO",font.size=9)
tiff("Plots/GO_bar_Day5_vs_control_up.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(fit1)
dev.off()
#GSEA

# start the GSEA with DESEQ2 output
KO_df_GSEA <- KO_df[KO_df$baseMean>50,]
KO_df_GSEA <- KO_df_GSEA[order(-KO_df_GSEA$stat),]
gene_list_up <- KO_df_GSEA$stat
names(gene_list_up) <- rownames(KO_df_GSEA)
gse <- gseGO(gene_list_up,
             ont = "BP",
             OrgDb = "org.Hs.eg.db",
             keyType = "ENSEMBL",
             eps = 1e-300)
gse_df <- as.data.frame(gse)

#make a pubmed trends plot
terms <- gse$Description[1:3]
pmcplot(terms,2010:2020,proportion = F)+theme_cowplot()
#make a dotplot
gsea_dot <- dotplot(gse,showCategory=20,split=".sign",font.size=7) +
  facet_grid(.~.sign) + ggtitle("Control vs shRNA")
tiff("Plots/GSEA_dot_control_vs_shRNA.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dot)
dev.off()




# start the GSEA with DESEQ2 output
KO_df_GSEA <- KOTreatment_df[KOTreatment_df$baseMean>50,]
KO_df_GSEA <- KO_df_GSEA[order(-KO_df_GSEA$stat),]
gene_list_up <- KO_df_GSEA$stat
names(gene_list_up) <- rownames(KO_df_GSEA)
gse <- gseGO(gene_list_up,
             ont = "BP",
             OrgDb = "org.Hs.eg.db",
             keyType = "ENSEMBL",
             eps = 1e-300)
gse_df <- as.data.frame(gse)

#make a dotplot
gsea_dot <- dotplot(gse,showCategory=20,split=".sign",font.size=6) +
  facet_grid(.~.sign) + ggtitle("Ifn vs Dox")
tiff("Plots/GSEA_dot_KO_ifn.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dot)
dev.off()



#save gsea plots
g1 <- gseaplot(gse,geneSetID = 15,by='all',title = gse_df$Description[15])+theme_cowplot()
tiff("Plots/ER to Golgi vesicle mediated transport.tiff",width = 35,height = 21,units = 'cm',res = 300)
print(g1)
dev.off()


#make an enrichment map
x2 <- pairwise_termsim(gse)
emapplot(x2,showCategory = 40)

#get list of genes from inf gamma pathway
hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

deg_table_GSEA <- KOTreatment_df[KOTreatment_df$baseMean>50,]
deg_table_GSEA <- deg_table_GSEA[order(-deg_table_GSEA$stat),]
#subtract DEGs from control vs shRNA
deg_table_GSEA <- deg_table_GSEA[rownames(deg_table_GSEA) %nin% rownames(KO_sig),]
gene_list_gsea <- data.frame(stat = deg_table_GSEA$stat, symbol = deg_table_GSEA$GeneName)
res2 <- gene_list_gsea %>% 
  na.omit() %>%
  distinct() %>%
  group_by(symbol) %>%
  summarize(stat=mean(stat))
res2 <- deframe(res2)
#names(gene_list_gsea) <- deg_table_GSEA$GeneName
gene_list_gsea <- gene_list_gsea %>% 
  na.omit() %>%
  distinct()

#read in the pathways
pathways.hallmark <- fgsea::gmtPathways("~/Downloads/h.all.v2023.1.Hs.symbols.gmt")

fgseaRes <- fgsea::fgsea(pathways = pathways.hallmark, stats = res2, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  write.table("Analysis/fgsea_pathway_enrichment_only_shRNAIfn.tsv", sep = "\t", quote = F, row.names = F)
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

search_kegg_organism("Homo sapiens", by = "scientific_name")
kk <- enrichKEGG(gene = rownames(KO_sig), organism = 'hsa',keyType = "ncbi-geneid")

#HALLMARK_GSEA Plots
#to run with gene names
h_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name,human_gene_symbol)
#to run with gene IDs
h_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name,human_ensembl_gene,human_gene_symbol)

deg_table_GSEA <- KOTreatment_df[KOTreatment_df$baseMean>50,]
deg_table_GSEA <- deg_table_GSEA[order(-deg_table_GSEA$stat),]
gene_list_gsea <- deg_table_GSEA$stat
#to run with gene names
names(gene_list_gsea) <- deg_table_GSEA$GeneName
#to run with gene IDs
names(gene_list_gsea) <- rownames(deg_table_GSEA)

#Gsea for hallmark pathways
gsea_hallmark <- clusterProfiler::GSEA(gene_list_gsea,
                                       minGSSize = 15,
                                       maxGSSize = 500,
                                       pAdjustMethod = "BH",
                                       seed = TRUE,
                                       eps = 1e-100,
                                       nPermSimple=10000,
                                       by = "fgsea",
                                       TERM2GENE = h_gene_set)


gsea_result <- as.data.frame(gsea_hallmark)
#write the results
write.table(gsea_result, file = "Analysis/Control_Ifn_vs_Dox_Hallmark_genes.tsv",
            quote = F,sep = "\t")

gsea_result_summary <- gsea_result[,c(1,3,5:8)]
gsea_result_summary$adjPvalue <- ifelse(gsea_result_summary$p.adjust <= 0.05, "significant", "non-significant")
cols <- c("non-significant" ="grey","significant"="red")
#general Plot
gsea_hallmark_plot <- ggplot(gsea_result_summary, aes(reorder(ID, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")+
  ggtitle("Control Ifn vs Ifn+Dox Hallmark Pathways")+
  theme_cowplot()
tiff("Plots/Hallmark_Gsea_bar_dox_vs_ifn.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_hallmark_plot)
dev.off()

gsea_dot <- enrichplot::dotplot(gsea_hallmark,showCategory=20,split=".sign",font.size=5) +
  facet_grid(.~.sign)+
  ggtitle("Hallmark GSEA Dotplot Control Ifn vs Ifn+Dox")
tiff("Plots/Hallmark_GSEA_dot_dox_vs_ifn.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dot)
dev.off()

## hallmark pathways gsea plots

g1 <- enrichplot::gseaplot2(gse,geneSetID = gse@result$ID[15],title = gse@result$Description[15])
tiff("Plots/Control_vs_Day5/chromatin_remodelling.tiff",width = 35,height = 21,units = 'cm',res = 300)
print(g1)
dev.off()



#list genes from Wnt signalling pathway
Wnt_genes <- data.frame(symbol = hs_hallmark_sets[hs_hallmark_sets$gs_name=="HALLMARK_WNT_BETA_CATENIN_SIGNALING","human_gene_symbol"],
                        ID = hs_hallmark_sets[hs_hallmark_sets$gs_name=="HALLMARK_WNT_BETA_CATENIN_SIGNALING","human_ensembl_gene"])
Ifn_gamma_genes <- data.frame(symbol = hs_hallmark_sets[hs_hallmark_sets$gs_name=="HALLMARK_INTERFERON_GAMMA_RESPONSE","human_gene_symbol"],
                              ID = hs_hallmark_sets[hs_hallmark_sets$gs_name=="HALLMARK_INTERFERON_GAMMA_RESPONSE","human_ensembl_gene"])

T_cell_genes <- Plots[[13]]@result[4,"geneID"]
T_cell_genes <- strsplit(T_cell_genes, split = "/")[[1]]

gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = T_cell_genes, mart= mart)
T_cell_genes_df <- data.frame(symbol=gene_IDs$hgnc_symbol,
                              ID = gene_IDs$ensembl_gene_id)