#OSCC_1d_GSE23558 32samples，###Agilent

library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)
####Agilent 芯片数据:GSE23558
raw_dir = "D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\GSE23558"
raw_datas = paste0(raw_dir,"/",dir(raw_dir))
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\GSE23558")
GSE23558_raw <- read.maimages(raw_datas,
                              source="agilent", 
                              green.only=TRUE,
                              other.columns = "gIsWellAboveBG")
dim(GSE23558_raw)

##背景校正
GSE23558_1 <- backgroundCorrect(GSE23558_raw,method = "normexp") 
##标准化(quantile normalization)
GSE23558_1 <- normalizeBetweenArrays(GSE23558_1, method="quantile")
class(GSE23558_1)

##基因过滤
##找到对照探针
Control <- GSE23558_1$genes$ControlType==1L
table(Control)

##找到探针名是空白
NoSymbol <- is.na(GSE23558_1$genes$ProbeName)
table(NoSymbol)

##找到低表达量的探针，超过一半数量的样本的表达量大于背景噪音，通过
IsExpr <- rowSums(GSE23558_1$other$gIsWellAboveBG > 0) >= 8
table(IsExpr)

#最后保留的是，非对照探针，且探针名不空白，且超过一半数量的样本的表达量大于背景噪音
GSE23558_1_filt <- GSE23558_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE23558_1_filt)

#最后，得到表达矩阵
GSE23558 = GSE23558_1_filt@.Data[[1]]
boxplot(GSE23558) #看一下基本情况

#给列名命名为样本编号
colnames(GSE23558) = str_extract(colnames(GSE23558),"GSM\\d*")
GSE23558[1:2,1:2]
#给行名命名为探针编号
rownames(GSE23558)= GSE23558_1_filt$genes$ProbeName

save(GSE23558,file = "GSE23558.RData")

##芯片注释

###anno the gene symbol
GSE23558_anno <- as.data.frame(GSE23558)
GSE23558_anno$ID<-rownames(GSE23558_anno)
gpl<-getGEO("GPL10526",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','SYMBOL')])
GSE23558_anno<-merge(x=GSE23558_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE23558_anno <- GSE23558_anno[,-1]
###The expression levels of the same gene names were averaged
GSE23558_anno<- aggregate(GSE23558_anno,by = list(GSE23558_anno$`SYMBOL`),FUN = mean)
head(GSE23558_anno)

GSE23558_anno <- GSE23558_anno[-c(1),]##blank gene name was dropped
rownames(GSE23558_anno) <- GSE23558_anno$Group.1
GSE23558_anno <- GSE23558_anno[,-c(1,98)]
GSE23558_anno[1:2,1:2]
boxplot(GSE23558_anno)
save(GSE23558_anno,file = "GSE23558_anno.RData")

##match the metadata with the expression data
pheno_GSE23558$sample_type = c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","N","N","N","N")
matchSN = match(colnames(GSE23558_anno), rownames(pheno_GSE23558))
GSE23558_anno = GSE23558_anno[,matchSN]
save(file = "GSE23558.RData",GSE23558_anno,pheno_GSE23558)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE23558_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE23558$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE23558_anno)[outliers]); print(table(outliers))
GSE23558_anno = GSE23558_anno[,!outliers]
pheno_GSE23558 = pheno_GSE23558[!outliers,]

#clean the datMet
pheno_GSE23558$sample_type = factor(pheno_GSE23558$sample_type,levels = c("N","T"))
pheno_GSE23558$Group = as.numeric(ifelse(pheno_GSE23558$sample_type=="T","1","0"))
pheno_GSE23558$gender = as.factor(pheno_GSE23558$gender)
pheno_GSE23558$age = as.numeric(pheno_GSE23558$age)
pheno_GSE23558$t_stage = as.numeric(pheno_GSE23558$stage)




save(file = "GSE23558_remove.RData",GSE23558_anno,pheno_GSE23558)
















