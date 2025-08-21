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

###HNSC GSE201777 31samples
##选择工作路径
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE201777")
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

##match the metadata with the expression data
matchSN = match(colnames(GSE201777_anno), rownames(pheno_GSE201777))
GSE201777_anno = GSE201777_anno[,matchSN]
save(file = "GSE201777.RData",GSE201777_anno,pheno_GSE201777)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE201777_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE201777$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE201777_anno)[outliers]); print(table(outliers))
GSE201777_anno = GSE201777_anno[,!outliers]
pheno_GSE201777 = pheno_GSE201777[!outliers,]

#clean the datMet
pheno_GSE201777$sample_type = factor(pheno_GSE201777$sample_type,levels = c("N","T"))
pheno_GSE201777$Group = as.numeric(ifelse(pheno_GSE201777$sample_type=="T","1","0"))



####回归除了诊断组（tumor normal）之外所有未知协变量的影响 
mod = model.matrix(~as.factor(Group),data = pheno_GSE201777)
mod0 = model.matrix(~1,data = pheno_GSE201777)
n.sv = num.sv(GSE201777_anno, mod, method="be")
GSE201777_anno = as.matrix(GSE201777_anno)

svobj = sva(GSE201777_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE201777$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 31)

X = svobj$sv
Y = GSE201777_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE201777_anno = GSE201777_anno-t(to_regress)

boxplot(GSE201777_anno)

save(file = "GSE201777.RData",GSE201777_anno,pheno_GSE201777)
