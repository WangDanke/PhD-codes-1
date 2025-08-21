##combined ESCC 7datasets 
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(stringr)
library(GEOquery)
library(ggplot2)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\ESCC")

pheno_GSE17351$batch = "GSE17351"
pheno_GSE20347$batch = "GSE20347"
pheno_GSE38129$batch = "GSE38129"
pheno_GSE23400_96$batch = "GSE23400"
pheno_GSE77861$batch = "GSE77861"
pheno_GSE161533$batch = "GSE161533"
pheno_GSE100942$batch = "GSE100942"


sel_phe = c("Sample_title","ID","sample_type","organ","batch")

pheno_GSE17351 = pheno_GSE17351[,sel_phe]
pheno_GSE20347 = pheno_GSE20347[,sel_phe]
pheno_GSE38129 = pheno_GSE38129[,sel_phe]
pheno_GSE23400 = pheno_GSE23400_96[,sel_phe]
pheno_GSE100942 = pheno_GSE100942[,sel_phe]
pheno_GSE161533 = pheno_GSE161533[,sel_phe]
pheno_GSE77861 = pheno_GSE77861[,sel_phe]


AA = intersect(rownames(GSE17351_anno),rownames(GSE20347_anno))
BB = intersect(rownames(GSE38129_anno),rownames(GSE23400_96_anno))
CC = intersect(rownames(GSE100942_anno),rownames(GSE161533_anno))


EE = intersect(AA,BB)
FF = intersect(CC,GSE77861_anno)
GG = intersect(EE,FF)


GSE17351_anno = GSE17351_anno[GG,]
GSE20347_anno = GSE20347_anno[GG,]
GSE38129_anno = GSE38129_anno[GG,]
GSE23400_anno = GSE23400_96_anno[GG,]
GSE161533_anno = GSE161533_anno[GG,]
GSE100942_anno = GSE100942_anno[GG,]
GSE77861_anno = GSE77861_anno[GG,]


##合并
ESCC_EXPR = cbind(GSE17351_anno,GSE20347_anno,GSE38129_anno,GSE23400_anno,
                  GSE161533_anno, GSE100942_anno,GSE77861_anno)
pheno_ESCC = rbind(pheno_GSE17351,pheno_GSE20347,pheno_GSE38129,pheno_GSE23400,
                   pheno_GSE100942,pheno_GSE161533,pheno_GSE77861)

##match
matchSN = match(colnames(ESCC_EXPR), rownames(pheno_ESCC))
ESCC_EXPR = ESCC_EXPR[,matchSN]
save(file = "ESCC_EXPR_7_datasets.RData",ESCC_EXPR,pheno_ESCC)


#校正之前看一下分布
###clean the metadata
pheno_ESCC$sample_type = factor(pheno_ESCC$sample_type,levels = c("N","T"))
pheno_ESCC$Group = as.numeric(ifelse(pheno_ESCC$sample_type=="T","1","0"))

boxplot(ESCC_EXPR, range = 0, col= as.numeric(pheno_ESCC$Group), main ="Boxplot", ylab = "Intensity")
#批次校正
mod = model.matrix(~sample_type, data=pheno_ESCC)
batch = factor(pheno_ESCC$batch)
ESCC_EXPR = ComBat(ESCC_EXPR, batch=batch, mod=mod)
ESCC_EXPR = as.data.frame(ESCC_EXPR)
##
boxplot(ESCC_EXPR, range = 0, col= as.numeric(pheno_ESCC$sample_type), main ="Boxplot", ylab = "Intensity")
library("FactoMineR")
library("factoextra")
ddb.pca <- PCA(t(ESCC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_ESCC$sample_type, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)



###clean the metadata

######校正未知协变量
mod = model.matrix(~as.factor(Group),data = pheno_ESCC)
mod0 = model.matrix(~1,data = pheno_ESCC)
n.sv = num.sv(ESCC_EXPR, mod, method="be")
ESCC_EXPR = as.matrix(ESCC_EXPR)

svobj = sva(ESCC_EXPR,mod, mod0,n.sv=n.sv)
cor = cor(pheno_ESCC$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 274)

X = svobj$sv
Y = ESCC_EXPR
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
ESCC_EXPR = ESCC_EXPR-t(to_regress)

##After QC plot
boxplot(ESCC_EXPR, range = 0, col= as.numeric(pheno_ESCC$sample_type), main ="Boxplot", ylab = "Intensity")
legend("topright", legend=c("N", "T"), col = 1:2, pch=19)

# Histogram
i = 1; plot(density((ESCC_EXPR[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(ESCC_EXPR)[2])
  lines(density((ESCC_EXPR[,i]), na.rm=T), col =as.numeric(pheno_ESCC$sample_type))
legend("topright", levels(pheno_ESCC$sample_type), cex=0.7, text.col = 1:2)

ddb.pca <- PCA(t(ESCC_EXPR), graph = F )
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno_ESCC$sample_type, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "Groups",
             range=0
)

save(ESCC_EXPR,pheno_ESCC,file = "ESCC_8_datasets_combat_remove_regress.RData")





