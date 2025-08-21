BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("biomaRt")
BiocManager::install("sva")
BiocManager::install("stringi")
BiocManager::install("GEOquery")
BiocManager::install("ggplot2")
BiocManager::install("FactoMineR")
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(GEOquery)
library(ggplot2)
library(FactoMineR)

##HNSC
###Affymetrix 芯片数据：GSE3292,GSE6631，GSE23036，GSE201777
##选择工作路径
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE3292")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE201777_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE201777_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE201777_raw)<-stri_sub(sampleNames(GSE201777_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE201777_raw_rma <- rma(GSE201777_raw)

#找到表达量，定义为最后需要的表达数据集
GSE201777 <- exprs(GSE201777_raw_rma) 


##芯片注释
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE201777")
###anno the gene symbol
GSE201777_anno <- as.data.frame(GSE201777)
GSE201777_anno$ID<-rownames(GSE201777_anno)
gpl571<-getGEO("GPL15207",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE201777_anno<-merge(x=GSE201777_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE201777_anno <- GSE201777_anno[,-1]
###The expression levels of the same gene names were averaged
GSE201777_anno<- aggregate(GSE201777_anno,by = list(GSE201777_anno$`Gene Symbol`),FUN = mean)
head(GSE201777_anno)

GSE201777_anno <- GSE201777_anno[-c(1,2),]##blank gene name was dropped
rownames(GSE201777_anno) <- GSE201777_anno$Group.1
GSE201777_anno <- GSE201777_anno[,-c(1,49)]
GSE201777_anno[1:2,1:2]
boxplot(GSE201777_anno)
save(GSE201777_anno,file = "GSE201777_anno.RData")




######
####Agilent 芯片数据:GSE172120
raw_dir = "D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE172120"
raw_datas = paste0(raw_dir,"/",dir(raw_dir))

GSE172120_raw <- read.maimages(raw_datas,
                   source="agilent", 
                   green.only=TRUE,
                   other.columns = "gIsWellAboveBG")
dim(GSE172120_raw)

##背景校正
GSE172120_1 <- backgroundCorrect(GSE172120_raw,method = "normexp") 
##标准化(quantile normalization)
GSE172120_1 <- normalizeBetweenArrays(GSE172120_1, method="quantile")
class(GSE172120_1)

##基因过滤
##找到对照探针
Control <- GSE172120_1$genes$ControlType==1L
table(Control)

##找到探针名是空白
NoSymbol <- is.na(GSE172120_1$genes$ProbeName)
table(NoSymbol)

##找到低表达量的探针，超过一半数量的样本的表达量大于背景噪音，通过
IsExpr <- rowSums(GSE172120_1$other$gIsWellAboveBG > 0) >= 4
table(IsExpr)

#最后保留的是，非对照探针，且探针名不空白，且超过一半数量的样本的表达量大于背景噪音
GSE172120_1_filt <- GSE172120_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE172120_1_filt)

#最后，得到表达矩阵
GSE172120 = GSE172120_1_filt@.Data[[1]]
boxplot(GSE172120) #看一下基本情况

#给列名命名为样本编号
colnames(GSE172120) = str_extract(colnames(GSE172120),"GSM\\d*")
GSE172120[1:2,1:2]
#给行名命名为探针编号
rownames(GSE172120)= GSE172120_1_filt$genes$ProbeName

save(GSE172120,file = "GSE172120.RData")

##芯片注释
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE172120")
###anno the gene symbol
GSE172120_anno <- as.data.frame(GSE172120)
GSE172120_anno$ID<-rownames(GSE172120_anno)
gpl571<-getGEO("GPL21185",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','GENE_SYMBOL')])
GSE172120_anno<-merge(x=GSE172120_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE172120_anno <- GSE172120_anno[,-1]
###The expression levels of the same gene names were averaged
GSE172120_anno<- aggregate(GSE172120_anno,by = list(GSE172120_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE172120_anno)

GSE172120_anno <- GSE172120_anno[-c(1),]##blank gene name was dropped
rownames(GSE172120_anno) <- GSE172120_anno$Group.1
GSE172120_anno <- GSE172120_anno[,-c(1,10)]
GSE172120_anno[1:2,1:2]
boxplot(GSE172120_anno)
save(GSE172120_anno,file = "GSE172120_anno.RData")
