#OSCC_1a_GSE37991 80samples, 40patients' paired tumor and normal samples

rm(list=ls());
options(stringsAsFactors=F)

suppressPackageStartupMessages(T)

library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)


##illumina 芯片数据
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\（待定）GSE37991")
getwd()

#read raw data
Y  = readBGX("GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx.gz")
raw_data <- read.delim("GSE37991_non-normalized.txt")

raw_EXP = raw_data[,seq(0,ncol(raw_data),2)]
rownames(raw_EXP) = raw_data$ID_REF
raw_PVALUE = raw_data[,seq(1,ncol(raw_data),2)]
rownames(raw_PVALUE) = raw_data$ID_REF
raw_PVALUE = raw_PVALUE[,-1]
##filt probes ##保留在一半的样本中，显著检测出来的芯片
IsExpr <- rowSums(raw_PVALUE < 0.05) >= 40
table(IsExpr)
raw_EXP_filt = raw_EXP[IsExpr,]

#background correct
raw_EXP_BC = backgroundCorrect(raw_EXP_filt, method="normexp")

# quantile
raw_EXP_QN<- normalizeBetweenArrays(raw_EXP_BC, method="quantile")


##log2 transform
dataExpr = as.data.frame(log2(raw_EXP_QN)) 

##match the metadata with the expression data
datMeta = pheno_GSE37991
AA = as.data.frame(str_split_fixed(datMeta$Sample_title, " ",2))
datMeta$t_N = AA$V1
matchSN = match(colnames(dataExpr), datMeta$t_N)
dataExpr = dataExpr[,matchSN]
colnames(dataExpr) = datMeta$ID

##annota
iddata<-as.data.frame(Y$probes[,c(14,12)])
colnames(iddata) = c("ID","gene")
dataExpr$ID = rownames(dataExpr)
GSE37991_anno<-merge(x=dataExpr,y=iddata,by='ID',all.x=T,all.y=F)
GSE37991_anno = GSE37991_anno[, -1]

###The expression levels of the same gene names were averaged
GSE37991_anno<- aggregate(GSE37991_anno,by = list(GSE37991_anno$`gene`),FUN = mean)
head(GSE37991_anno)
rownames(GSE37991_anno) <- GSE37991_anno$Group.1
GSE37991_anno = GSE37991_anno[,-c(1,82)]

save(file = "GSE37991.RData",pheno_GSE37991,GSE37991_anno)

##match the metadata with the expression data
matchSN = match(colnames(GSE37991_anno), rownames(pheno_GSE37991))
GSE37991_anno = GSE37991_anno[,matchSN]
save(file = "GSE37991.RData",GSE37991_anno,pheno_GSE37991)


##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE37991_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE37991$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE37991_anno)[outliers]); print(table(outliers))
GSE37991_anno = GSE37991_anno[,!outliers]
pheno_GSE37991 = pheno_GSE37991[!outliers,]

save(file = "GSE37991_remove.RData",pheno_GSE37991,GSE37991_anno)

















