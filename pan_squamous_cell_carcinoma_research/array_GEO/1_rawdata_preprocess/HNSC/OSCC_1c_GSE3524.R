#OSCC_1c_GSE3524 20samples，###Affymetrix
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

##选择工作路径
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\GSE3524")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE3524_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE3524_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE3524_raw)<-stri_sub(sampleNames(GSE3524_raw),1,8)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE3524_raw_rma <- rma(GSE3524_raw)

#找到表达量，定义为最后需要的表达数据集
GSE3524 <- exprs(GSE3524_raw_rma) 


##芯片注释

###anno the gene symbol
GSE3524_anno <- as.data.frame(GSE3524)
GSE3524_anno$ID<-rownames(GSE3524_anno)
gpl<-getGEO("GPL96",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE3524_anno<-merge(x=GSE3524_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE3524_anno <- GSE3524_anno[,-1]
###The expression levels of the same gene names were averaged
GSE3524_anno<- aggregate(GSE3524_anno,by = list(GSE3524_anno$`Gene Symbol`),FUN = mean)
head(GSE3524_anno)

GSE3524_anno <- GSE3524_anno[-c(1),]##blank gene name was dropped
rownames(GSE3524_anno) <- GSE3524_anno$Group.1
GSE3524_anno <- GSE3524_anno[,-c(1,22)]
GSE3524_anno[1:2,1:2]
boxplot(GSE3524_anno)
save(GSE3524_anno,file = "GSE3524_anno.RData")

##match the metadata with the expression data
pheno_GSE3524$sample_type = c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","N","N","N")
matchSN = match(colnames(GSE3524_anno), rownames(pheno_GSE3524))
GSE3524_anno = GSE3524_anno[,matchSN]
save(file = "GSE3524.RData",GSE3524_anno,pheno_GSE3524)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE3524_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE3524$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE3524_anno)[outliers]); print(table(outliers))
GSE3524_anno = GSE3524_anno[,!outliers]
pheno_GSE3524 = pheno_GSE3524[!outliers,]

#clean the datMet
pheno_GSE3524$sample_type = factor(pheno_GSE3524$sample_type,levels = c("N","T"))
pheno_GSE3524$Group = as.numeric(ifelse(pheno_GSE3524$sample_type=="T","1","0"))
pheno_GSE3524$gender = as.factor(pheno_GSE3524$gender)
pheno_GSE3524$age = as.numeric(pheno_GSE3524$age)
pheno_GSE3524$t_stage = as.numeric(substr(pheno_GSE3524$tnm_stage,2,2))



save(file = "GSE3524_remove.RData",GSE3524_anno,pheno_GSE3524)





