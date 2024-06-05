#plot for Kim

#comparision of normalised counts for GRASSLND vs MLANA

#import the counts file

library(data.table)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(cowplot)


nm_counts <- fread("~/Documents/GSE134432_bulkRNASeq_33MelanomaBaselines_normalized_counts.tsv", header = F)
imp <- nm_counts[nm_counts$V1 %in% c("RNF144A-AS1", "MLANA")]
imp <- as.data.frame(t(imp))
colnames(imp) <- c("MelanA", "GRASLND")
imp <- imp[-1,]
imp$MELAN.A <- as.numeric(imp$MELAN.A)
imp$GRASSLND <- as.numeric(imp$GRASSLND)
p1 <- ggplot(imp, aes(x=GRASLND,y=MelanA))+
  geom_point(size=2,alpha=0.5)+
  theme_cowplot()+
  ylab("MelanA Normalized Counts")+
  ggtitle("GRASLND vs. MelanA")+
  theme(plot.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t(30))))+
  coord_cartesian(ylim = c(0,30000), clip = "off")+
  xlab(expression(paste("         GRASLND Normalized counts \n Spearmann p=2.8466e-06 coeff=0.72(n=33)")))
pdf("Norm_counts_grnlnd_darker.pdf",width = 35,height = 21,units = 'cm',res = 300)
plot(p1)
dev.off()


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