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

###HNSC GSE3292 36samples
##选择工作路径
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE3292")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE3292_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE3292_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE3292_raw)<-stri_sub(sampleNames(GSE3292_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE3292_raw_rma <- rma(GSE3292_raw)

#找到表达量，定义为最后需要的表达数据集
GSE3292 <- exprs(GSE3292_raw_rma) 


##芯片注释
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE3292")
###anno the gene symbol
GSE3292_anno <- as.data.frame(GSE3292)
GSE3292_anno$ID<-rownames(GSE3292_anno)
gpl571<-getGEO("GPL15207",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE3292_anno<-merge(x=GSE3292_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE3292_anno <- GSE3292_anno[,-1]
###The expression levels of the same gene names were averaged
GSE3292_anno<- aggregate(GSE3292_anno,by = list(GSE3292_anno$`Gene Symbol`),FUN = mean)
head(GSE3292_anno)

GSE3292_anno <- GSE3292_anno[-c(1,2),]##blank gene name was dropped
rownames(GSE3292_anno) <- GSE3292_anno$Group.1
GSE3292_anno <- GSE3292_anno[,-c(1,49)]
GSE3292_anno[1:2,1:2]

##match the metadata with the expression data
matchSN = match(colnames(GSE3292_anno), rownames(pheno_GSE3292))
GSE3292_anno = GSE3292_anno[,matchSN]
save(file = "GSE3292.RData",GSE3292_anno,pheno_GSE3292)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE3292_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE3292$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE3292_anno)[outliers]); print(table(outliers))
GSE3292_anno = GSE3292_anno[,!outliers]
pheno_GSE3292 = pheno_GSE3292[!outliers,]

save(file = "GSE3292_remove.RData",GSE3292_anno,pheno_GSE3292)




