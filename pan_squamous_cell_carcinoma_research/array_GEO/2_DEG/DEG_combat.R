############use the function Combat to remove batch effects
###DEG analysis
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio);library(oligo)

###combined_all_types_SCC,

setwd("D:\\鳞癌分析\\泛鳞癌array数据\\DEG_combat")
AA = intersect(rownames(cSCC_EXPR),rownames(Cer_SCC_EXPR))
BB = intersect(rownames(HNSC_EXPR),rownames(ESCC_EXPR))
CC = intersect(AA,rownames(LUSC_EXPR))
DD = intersect(CC,BB)

cSCC_EXPR = cSCC_EXPR[DD,]
Cer_SCC_EXPR = Cer_SCC_EXPR[DD,]
ESCC_EXPR = ESCC_EXPR[DD,]
HNSC_EXPR = HNSC_EXPR[DD,]
LUSC_EXPR = LUSC_EXPR[DD,]

SCC_EXPR = cbind(cSCC_EXPR,Cer_SCC_EXPR,ESCC_EXPR,HNSC_EXPR,LUSC_EXPR)

sel = c("Sample_title","ID","sample_type","organ","batch","Group")
pheno_Cer_SCC = pheno_Cer_SCC[,sel]
pheno_cSCC = pheno_cSCC[,sel]
pheno_ESCC = pheno_ESCC[,sel]
pheno_HNSC = pheno_HNSC[,sel]
pheno_LUSC = pheno_LUSC[,sel]
pheno_SCC = rbind(pheno_Cer_SCC,pheno_cSCC,pheno_ESCC,pheno_HNSC,pheno_LUSC)

save(pheno_SCC,SCC_EXPR,file = "SCC_combined_697_samples.RData")

#####Cervical
pheno_Cervical = subset(pheno_SCC,pheno_SCC$organ=="ovary")
Cervical_EXPR = as.data.frame(SCC_EXPR[,rownames(pheno_Cervical)])

###use the limma package to find the DEGs
#### Differential expression analysis of individual subtype and paired normal tissues
data = Cervical_EXPR
group <- as.character(pheno_Cervical$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_cervical=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_cervical,file = "allDiff_cervical.Rdata",quote = F, row.names = T)


#####ESCC
pheno_ESCC = subset(pheno_SCC,pheno_SCC$organ=="esophagus")
ESCC_EXPR = as.data.frame(SCC_EXPR[,rownames(pheno_ESCC)])

###use the limma package to find the DEGs
#### Differential expression analysis of individual subtype and paired normal tissues
data = ESCC_EXPR
group <- as.character(pheno_ESCC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_ESCC=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_ESCC,file = "allDiff_ESCC.Rdata",quote = F, row.names = T)

#####cSCC
pheno_cSCC = subset(pheno_SCC,pheno_SCC$organ=="skin")
cSCC_EXPR = as.data.frame(SCC_EXPR[,rownames(pheno_cSCC)])

###use the limma package to find the DEGs
#### Differential expression analysis of individual subtype and paired normal tissues
data = cSCC_EXPR
group <- as.character(pheno_cSCC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_cSCC=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_cSCC,file = "allDiff_cSCC.Rdata",quote = F, row.names = T)

#####HNSC
pheno_HNSC_1 = subset(pheno_SCC,pheno_SCC$organ=="Oral")
pheno_HNSC_2 = subset(pheno_SCC,pheno_SCC$organ=="head and neck")
pheno_HNSC = rbind(pheno_HNSC_1,pheno_HNSC_2)
HNSC_EXPR = as.data.frame(SCC_EXPR[,rownames(pheno_HNSC)])

###use the limma package to find the DEGs
#### Differential expression analysis of individual subtype and paired normal tissues
data = HNSC_EXPR
group <- as.character(pheno_HNSC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_HNSC=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_HNSC,file = "allDiff_HNSC.Rdata",quote = F, row.names = T)


#####LUSC
pheno_LUSC = subset(pheno_SCC,pheno_SCC$organ=="Lung")
LUSC_EXPR = as.data.frame(SCC_EXPR[,rownames(pheno_LUSC)])

###use the limma package to find the DEGs
#### Differential expression analysis of individual subtype and paired normal tissues
data = LUSC_EXPR
group <- as.character(pheno_LUSC$sample_type)
group <- factor(group,levels = c("N","T"),ordered = F)
design <- model.matrix(~group)

fit <- lmFit(data,design)
fit2 <- eBayes(fit)
allDiff_LUSC=topTable(fit2,adjust='fdr',coef=2,number=Inf)

save(allDiff_LUSC,file = "allDiff_LUSC.Rdata",quote = F, row.names = T)



##############相关性分析######
AA = rownames(allDiff_cervical)
allDiff_cSCC = allDiff_cSCC[AA,]
allDiff_ESCC = allDiff_ESCC[AA,]
allDiff_HNSC = allDiff_HNSC[AA,]
allDiff_LUSC = allDiff_LUSC[AA,]

allDiff = cbind(allDiff_cervical$logFC,allDiff_cSCC$logFC,allDiff_ESCC$logFC,
                allDiff_HNSC$logFC,allDiff_LUSC$logFC)

colnames(allDiff) = c("cervical","skin","esophagus","head and neck","lung")

rownames(allDiff) = AA
allDiff = as.data.frame(allDiff)

save(allDiff,file = "allDiff.RData")

library(corrplot)
library(ggplot2)
library(ggpubr)
class(allDiff)

###相关性分析和可视化
##方法一：
tdc = cor(allDiff,method = "pearson")
testRes = cor.mtest(allDiff, method="pearson",conf.level = 0.95)

corrplot(tdc)
addcol <- colorRampPalette(c("red", "white", "blue"))
corrplot(tdc, method = "color", col = addcol(100), 
         tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt",
         p.mat = testRes$p, diag = T, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')
corrplot(tdc, method = "number", type = "upper",col = addcol(100), 
         tl.col = "black", tl.cex = 0.8, tl.pos = "n",order = 'AOE',
         add = T)

##可视化方法二：
library(PerformanceAnalytics)
chart.Correlation(allDiff)

##可视化方法三：
install.packages("GGally")
library(GGally)
ggpairs(allDiff)

library(corrplot)
#绘制一个下三角的热图
corrplot(tdc, type="upper", order="hclust", tl.col="black", tl.srt=45)







