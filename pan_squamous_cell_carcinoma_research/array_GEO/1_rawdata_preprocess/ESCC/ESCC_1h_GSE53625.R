#####GSE53625
#####358 samples

###load packages
BiocManager::install("biomaRt")
BiocManager::install("FactoMineR")
BiocManager::install("WGCNA")
BiocManager::install("oligo")

library(WGCNA)
library(affy)
library(limma)
library(biomaRt)
library(sva)
library(stringi)
library(GEOquery)
library(ggplot2)
library(FactoMineR)
library(nlme)

setwd("D:\\ESCC_Subtype_research\\GSE53625_RAW")
getwd()
####整理临床信息，GSE53625
dir <- choose.dir(default = "", caption = "")
setwd(dir)
getwd()
##预处理和整理数据
GSE53625<-read.csv("GSE53625_series_matrix.csv",header=T,sep=',')
##整理临床信息数据
GSE53625_pheno <- t(GSE53625[c(38:55,78),])
rownames(GSE53625_pheno) <- GSE53625_pheno[,19]
colnames(GSE53625_pheno) <- GSE53625_pheno[1,]
GSE53625_pheno <- GSE53625_pheno[-1,-19]
GSE53625_pheno <- as.data.frame(GSE53625_pheno)
#####
###整理临床信息
colnames(GSE53625_pheno)[1] <- "Patient_ID"
GSE53625_pheno$Patient_ID <- sub("^patient id:", "", GSE53625_pheno$Patient_ID)
colnames(GSE53625_pheno)[2] <- "Age"
GSE53625_pheno$Age <- sub("^age:", "", GSE53625_pheno$Age)
colnames(GSE53625_pheno)[3] <- "Sex"
GSE53625_pheno$Sex <- sub("^Sex:", "", GSE53625_pheno$Sex)
colnames(GSE53625_pheno)[4] <- "tobacco"
GSE53625_pheno$tobacco <- sub("^tobacco use:", "", GSE53625_pheno$tobacco)
colnames(GSE53625_pheno)[5] <- "alcohol"
GSE53625_pheno$alcohol <- sub("^alcohol use:", "", GSE53625_pheno$alcohol)
colnames(GSE53625_pheno)[6] <- "tumor_loation"
GSE53625_pheno$tumor_loation <- sub("^tumor loation:", "", GSE53625_pheno$tumor_loation)
colnames(GSE53625_pheno)[7] <- "tumor_grade"
GSE53625_pheno$tumor_grade <- sub("^tumor grade:", "", GSE53625_pheno$tumor_grade)
colnames(GSE53625_pheno)[8] <- "t_stage"
GSE53625_pheno$t_stage <- sub("^t stage:", "", GSE53625_pheno$t_stage)
colnames(GSE53625_pheno)[9] <- "n_stage"
GSE53625_pheno$n_stage <- sub("^n stage:", "", GSE53625_pheno$n_stage)

colnames(GSE53625_pheno)[10] <- "tnm_stage"
GSE53625_pheno$tnm_stage <- sub("^tnm stage:", "", GSE53625_pheno$tnm_stage)
colnames(GSE53625_pheno)[11] <- "Arrhythmia"
GSE53625_pheno$Arrhythmia <- sub("^arrhythmia:", "", GSE53625_pheno$Arrhythmia)
colnames(GSE53625_pheno)[12] <- "Pneumonia"
GSE53625_pheno$Pneumonia <- sub("^pneumonia:", "", GSE53625_pheno$Pneumonia)
colnames(GSE53625_pheno)[13] <- "Anastomotic_leak"
GSE53625_pheno$Anastomotic_leak <- sub("^anastomotic leak:", "", GSE53625_pheno$Anastomotic_leak)
colnames(GSE53625_pheno)[14] <- "Adjuvant_therapy"
GSE53625_pheno$Adjuvant_therapy <- sub("^adjuvant therapy:", "", GSE53625_pheno$Adjuvant_therapy)
colnames(GSE53625_pheno)[15] <- "Status"
GSE53625_pheno$Status <- sub("^death at fu:", "", GSE53625_pheno$Status)

colnames(GSE53625_pheno)[16] <- "Time"
GSE53625_pheno$Time <- sub("^survival time", "", GSE53625_pheno$Time)
GSE53625_pheno$Time_day <- as.numeric(substr(GSE53625_pheno$Time,10,14))*30
GSE53625_pheno <- GSE53625_pheno[,-17]

colnames(GSE53625_pheno)[17] <- "Sample_type"
GSE53625_pheno$Sample_type <- substr(GSE53625_pheno$Sample_type,8,14)

##use function choose.dir() to filter working directory
dir <- choose.dir(caption = "D:\\ESCC_Subtype_research\\GSE53625_RAW")##dir <- choose.dir(caption = "D:\\ESCCmRNA-from rawdata\\GSE53625_RAW")
setwd("D:\\ESCC_Subtype_research\\GSE53625_RAW")
getwd()

txt_files <- list.files("D:\\ESCC_Subtype_research\\GSE53625_RAW\\GSE53625_RAW", pattern = "\\.txt.gz", full.names = TRUE)


GSE53625_raw <- read.maimages(txt_files,
                               source="agilent", 
                               green.only=TRUE,
                               other.columns = "gIsWellAboveBG")
dim(GSE53625_raw)

##背景校正
GSE53625_1 <- backgroundCorrect(GSE53625_raw,method = "normexp") 
##标准化(quantile normalization)
GSE53625_1 <- normalizeBetweenArrays(GSE53625_1, method="quantile")
class(GSE53625_1)

##基因过滤
##找到对照探针
Control <- GSE53625_1$genes$ControlType==1L
table(Control)

##找到探针名是空白
NoSymbol <- is.na(GSE53625_1$genes$ProbeName)
table(NoSymbol)

##找到低表达量的探针，超过一半数量的样本的表达量大于背景噪音，通过
IsExpr <- rowSums(GSE53625_1$other$gIsWellAboveBG > 0) >= 179
table(IsExpr)

#最后保留的是，非对照探针，且探针名不空白，且超过一半数量的样本的表达量大于背景噪音
GSE53625_1_filt <- GSE53625_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE53625_1_filt)

#最后，得到表达矩阵
GSE53625 = GSE53625_1_filt@.Data[[1]]
boxplot(GSE53625) #看一下基本情况


#给列名命名为样本编号
colnames(GSE53625) = str_extract(colnames(GSE53625),"GSM\\d*")
GSE53625[1:2,1:2]
#给行名命名为探针编号
rownames(GSE53625)= GSE53625_1_filt$genes$ProbeName

########################################
###芯片注释，有点绕，先使用GPL18109的soft文件对应上芯片ID信息，再使用idmap将芯片ID信息对应上gene symbol
gpl570<-getGEO("GPL18109",destdir=".")
iddata<-as.data.frame(Table(gpl570)[c('ID','SPOT_ID')])

GSE53625_anno = as.data.frame(GSE53625)
GSE53625_anno$SPOT_ID = rownames(GSE53625_anno)


GSE53625_anno<-merge(x=GSE53625_anno,y=iddata,by='SPOT_ID',all.x=T,all.y=F)
GSE53625_anno <- GSE53625_anno[,-1]

#得到探针对应的基因名字
library(AnnoProbe)
library(devtools)
probe2gene=idmap('GPL18109',type = 'pipe')
#展示前10条结果
head(probe2gene)
#probe_id  symbol
#1    80108 DDX11L1
#2    80108  WASH7P
#3     4320 DDX11L1
#4     4320  WASH7P
#5    97414 DDX11L1
#6    97414  WASH7P
colnames(probe2gene) = c("ID","gene_symbol")
GSE53625_anno<-merge(x=GSE53625_anno,y=probe2gene,by='ID',all.x=T,all.y=F)

GSE53625_anno <- GSE53625_anno[,-1]
GSE53625_anno<- aggregate(GSE53625_anno,by = list(GSE53625_anno$`gene_symbol`),FUN = mean)


rownames(GSE53625_anno) <- GSE53625_anno$Group.1

GSE53625_anno <- GSE53625_anno[,-c(1,360)]
head(GSE53625_anno)

####clean the datMet
pheno_GSE53625$sample_type = factor(pheno_GSE53625$sample_type,levels = c("N","T"))
pheno_GSE53625$Group = as.numeric(ifelse(pheno_GSE53625$sample_type=="T","1","0"))



##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE53625_EXP))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE53625$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE53625_EXP)[outliers]); print(table(outliers))
GSE53625_EXP = GSE53625_EXP[,!outliers]
pheno_GSE53625 = pheno_GSE53625[!outliers,]

save(file = "GSE53625_new.RData",GSE53625_anno,pheno_GSE53625)







