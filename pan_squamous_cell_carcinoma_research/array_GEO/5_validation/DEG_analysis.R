setwd("D:\\鳞癌分析\\泛鳞癌seq数据_标准参考基因组\\8_Combat_DEG_seq")
#######
###分别合并同一部位的数据集并进行差异表达分析
##合并不同批次数据,整理临床数据
load("D:/鳞癌分析/泛鳞癌array数据/7_DEG_combat/allDiff_5.RData")
genes = intersect(rownames(GSE223804_EXPR),rownames(allDiff_5))

CESC_EXPR_seq = cbind(GSE223804_EXPR[genes,],GSE87410_CESC_EXPR[genes,])
                      
pheno_CESC_seq = rbind(pheno_GSE223804[,c("Sample_title","ID","sample_type","organ","ID_SRR","batch")],
                       pheno_GSE87410_CESC[,c("Sample_title","ID","sample_type","organ","ID_SRR","batch")]
                       )
pheno_CESC_seq$organ = "CESC"
pheno_CESC_seq$sample_type = factor(pheno_CESC_seq$sample_type,levels = c("N","T"))
pheno_CESC_seq$Group = as.numeric(ifelse(pheno_CESC_seq$sample_type=="T","1","0"))
rownames(pheno_CESC_seq) = pheno_CESC_seq$ID_SRR

matchSN = match(colnames(CESC_EXPR_seq), rownames(pheno_CESC_seq))
CESC_EXPR_seq = CESC_EXPR_seq[,matchSN]
pheno_CESC_seq = pheno_CESC_seq[matchSN,]

save(pheno_CESC_seq,CESC_EXPR_seq,file = "CESC_EXPR_seq_TPM_log.RData")
###
library("FactoMineR")
library("factoextra")
ddb.pca <- PCA(t(CESC_EXPR_seq.combat), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_CESC_seq$sample_type, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "batch",
             range=0
)

library(limma)
##批次校正
mod = model.matrix(~sample_type, data=pheno_CESC_seq)
batch = factor(pheno_CESC_seq$batch)
CESC_EXPR_seq.combat = ComBat(CESC_EXPR_seq, batch=batch, mod=mod)

save(file = "CESC_EXPR_seq_combat_combined_2datasets.RData",CESC_EXPR_seq.combat,pheno_CESC_seq)

##After QC plot
boxplot(CESC_EXPR_seq.combat, range = 0, col= as.numeric(pheno_CESC_seq$sample_type), main ="Boxplot", ylab = "Intensity")
legend("topright", legend=c("N", "T"), col = 1:2, pch=19)

# Histogram
i = 1; plot(density((CESC_EXPR_seq.combat[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(CESC_EXPR_seq.combat)[2])
  lines(density((CESC_EXPR_seq.combat[,i]), na.rm=T), col =as.numeric(pheno_CESC_seq$sample_type))
legend("topright", levels(pheno_CESC_seq$sample_type), cex=0.7, text.col = 1:2)



##差异表达分析
data = HNSC_EXPR_seq.combat
group <- as.character(pheno_HNSC_seq$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_HNSC_seq=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_HNSC_seq,file = "allDiff_HNSC_seq.Rdata",quote = F, row.names = T)

###将两种数据类型得到的鳞癌之间的相关性值整理成correlation_both.csv
library(ggplot2)
library(ggpubr)
####两种类型数据最后得到的鳞癌相关性值的相关性
correlation = read.csv("correlation_both.csv",header = T,row.names = 1)
correlation = as.data.frame(t(correlation))
correlation$group = rownames(correlation)

ggplot(correlation, aes(x = Microarray, y = RNA_seq,)) +
  geom_point() +
  geom_smooth(method = "lm", color = "#CD8500",se = T, fill = "#EEC591") +
  stat_cor(method = "pearson", label.x = 0.3, label.y = 0.8) +
  theme(
    axis.text.x=element_text(size=10, hjust=1,color = "#828282"),
    axis.text.y=element_text(size=10, vjust=0.5,color = "#828282"), 
    legend.title = element_text(size=12,color = "grey"), 
    plot.title = element_text(size=12, face="bold",color = "#363636"),
    axis.title.x = element_text(size=12,color = "#363636"),
    axis.title.y = element_text(size=12, vjust=0.5,color = "#363636"),
    plot.margin=unit(c(2,2,16,2),"mm"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "#363636", fill = NA, size = 1)
  ) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",size = 0.55)+
  labs(title = "Microarray vs. RNA_seq",
       x = "Transcriptome correlation (Microarray)",
       y = "Transcriptome correlation (RNA_seq)")+
  coord_fixed(ratio = 1, xlim = c(0.2,0.8), ylim = c(0.2,0.8), expand = TRUE, clip = "off")+
  geom_text(aes(label = group), size = 3, vjust = -0.5, hjust = 0.1, color = "#CD8500") 
  
























