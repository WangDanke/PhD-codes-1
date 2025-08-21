##OSCC_1f_GSE87593 ,####Agilent

library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(ggplot2); library(Cairo); library(GEOquery)
library(biomaRt); library(sva)

raw_dir = "D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\GSE87593"
raw_datas = paste0(raw_dir,"/",dir(raw_dir))
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\GSE87593")
GSE87593_raw <- read.maimages(raw_datas,
                              source="agilent", 
                              green.only=TRUE,
                              other.columns = "gIsWellAboveBG")
dim(GSE87593_raw)

##背景校正
GSE87593_1 <- backgroundCorrect(GSE87593_raw,method = "normexp") 
##标准化(quantile normalization)
GSE87593_1 <- normalizeBetweenArrays(GSE87593_1, method="quantile")
class(GSE87593_1)

##基因过滤
##找到对照探针
Control <- GSE87593_1$genes$ControlType==1L
table(Control)

##找到探针名是空白
NoSymbol <- is.na(GSE87593_1$genes$ProbeName)
table(NoSymbol)

##找到低表达量的探针，超过一半数量的样本的表达量大于背景噪音，通过
IsExpr <- rowSums(GSE87593_1$other$gIsWellAboveBG > 0) >= 8
table(IsExpr)

#最后保留的是，非对照探针，且探针名不空白，且超过一半数量的样本的表达量大于背景噪音
GSE87593_1_filt <- GSE87593_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE87593_1_filt)

#最后，得到表达矩阵
GSE87593 = GSE87593_1_filt@.Data[[1]]
boxplot(GSE87593) #看一下基本情况

#给列名命名为样本编号
colnames(GSE87593) = str_extract(colnames(GSE87593),"GSM\\d*")
GSE87593[1:2,1:2]
#给行名命名为探针编号
rownames(GSE87593)= GSE87593_1_filt$genes$ProbeName

save(GSE87593,file = "GSE87593.RData")


##芯片注释
###anno the gene symbol
GSE87593_anno <- as.data.frame(GSE87593)
GSE87593_anno$ID<-rownames(GSE87593_anno)

#gpl<- getGEO(filename = "GPL14550"),找不到，直接在网页上下载了相关的注释文件
iddata = read.csv("GPL14550.csv",header = F)
iddata = iddata[-c(1:17),]
colnames(iddata) = iddata[1,]
iddata = iddata[-1,]
iddata = iddata[,c('ID','GENE_SYMBOL')]

GSE87593_anno<-merge(x=GSE87593_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE87593_anno <- GSE87593_anno[,-1]
###The expression levels of the same gene names were averaged
GSE87593_anno<- aggregate(GSE87593_anno,by = list(GSE87593_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE87593_anno)

GSE87593_anno <- GSE87593_anno[-c(1),]##blank gene name was dropped
rownames(GSE87593_anno) <- GSE87593_anno$Group.1
GSE87593_anno <- GSE87593_anno[,-c(1,18)]
GSE87593_anno[1:2,1:2]
boxplot(GSE87593_anno)
save(GSE87593_anno,file = "GSE87593_anno.RData")

##match the metadata with the expression data
matchSN = match(colnames(GSE87593_anno), rownames(pheno_GSE87593))
GSE87593_anno = GSE87593_anno[,matchSN]
save(file = "GSE87593.RData",GSE87593_anno,pheno_GSE87593)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE87593_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE87593$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE87593_anno)[outliers]); print(table(outliers))
GSE87593_anno = GSE87593_anno[,!outliers]
pheno_GSE87593 = pheno_GSE87593[!outliers,]

save(file = "GSE87593_remove.RData",GSE87593_anno,pheno_GSE87593)



















