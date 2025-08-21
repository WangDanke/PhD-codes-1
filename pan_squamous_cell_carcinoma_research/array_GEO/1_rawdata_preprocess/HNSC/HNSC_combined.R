###合并HNSC的五个数据集，并且做批次校正、移除离群值、校正未知协变量
AA = intersect(rownames(GSE23036_anno),rownames(GSE201777_anno))
BB = intersect(rownames(GSE23036_anno),rownames(GSE3292_anno))
CC = intersect(AA,BB)
DD = intersect(CC,rownames(GSE6631_anno))

GSE172120_anno = GSE172120_anno[DD,]
GSE201777_anno = GSE201777_anno[DD,]
GSE23036_anno = GSE23036_anno[DD,]
GSE3292_anno = GSE3292_anno[DD,]
GSE6631_anno = GSE6631_anno[DD,]

HNSC_expr = cbind(GSE172120_anno,GSE201777_anno,GSE23036_anno,GSE3292_anno,GSE6631_anno)
HNSC_expr = HNSC_expr[,rownames(pheno_HNSC)]


library("FactoMineR")
library("factoextra")

ddb.pca <- PCA(t(HNSC_expr), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_HNSC$batch, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)


table(datMeta$batch)
##批次校正
mod = model.matrix(~sample_type, data=datMeta)
batch = factor(datMeta$batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat
HNSC_expr = as.data.frame(datExpr)

##match the metadata with the expression data
matchSN = match(colnames(HNSC_expr), rownames(pheno_HNSC))
HNSC_expr = HNSC_expr[,matchSN]
save(file = "datExpr_HNSC_combat_combined_5datasets.RData",HNSC_expr,pheno_HNSC)

##After QC plot
boxplot(HNSC_expr, range = 0, col= as.numeric(pheno_HNSC$sample_type), main ="Boxplot", ylab = "Intensity")
legend("topright", legend=c("N", "T"), col = 1:2, pch=19)

# Histogram
i = 1; plot(density((HNSC_expr[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(HNSC_expr)[2])
  lines(density((HNSC_expr[,i]), na.rm=T), col =as.numeric(pheno_HNSC$sample_type))
legend("topright", levels(pheno_HNSC$sample_type), cex=0.7, text.col = 1:2)

save(HNSC_expr,pheno_HNSC,file = "combind_HNSC_5_dataset_combat_remove_regress_all.RData")







