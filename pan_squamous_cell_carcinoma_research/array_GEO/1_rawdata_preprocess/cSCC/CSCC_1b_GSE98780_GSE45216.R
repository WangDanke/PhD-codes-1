###CSCC_1b_GSE98780 #Affymrtrix
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)

setwd("D:\\鳞癌分析\\泛鳞癌array数据\\CSCC\\GSE98780_RAW")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE98780_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE98780_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE98780_raw)<-stri_sub(sampleNames(GSE98780_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE98780_raw_rma <- rma(GSE98780_raw)

#找到表达量，定义为最后需要的表达数据集
GSE98780 <- exprs(GSE98780_raw_rma) 

##将临床信息与表达量信息匹配起来
pheno_GSE98780 = pheno_GSE98780_GPL570[c(19:54),] #去除光化病样本
GSE98780 = GSE98780[,rownames(pheno_GSE98780)]


##芯片注释
###anno the gene symbol
GSE98780_anno <- as.data.frame(GSE98780)
GSE98780_anno$ID<-rownames(GSE98780_anno)
gpl<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE98780_anno<-merge(x=GSE98780_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE98780_anno <- GSE98780_anno[,-1]
###The expression levels of the same gene names were averaged
GSE98780_anno<- aggregate(GSE98780_anno,by = list(GSE98780_anno$`Gene Symbol`),FUN = mean)
head(GSE98780_anno)

GSE98780_anno <- GSE98780_anno[-c(1),]##blank gene name was dropped
rownames(GSE98780_anno) <- GSE98780_anno$Group.1
GSE98780_anno <- GSE98780_anno[,-c(1,38)]
GSE98780_anno[1:2,1:2]
boxplot(GSE98780_anno)

matchSN = match(colnames(GSE98780_anno), rownames(pheno_GSE98780))
GSE98780_anno = GSE98780_anno[,matchSN]
save(file = "GSE98780_match.RData",GSE98780_anno,pheno_GSE98780)

#clean the datMet
pheno_GSE98780$sample_type = factor(pheno_GSE98780$sample_type,levels = c("N","T"))
pheno_GSE98780$Group = as.numeric(ifelse(pheno_GSE98780$sample_type=="T","1","0"))

###处理GSE45216数据集
######
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\CSCC\\GSE45216_RAW")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE45216_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE45216_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE45216_raw)<-stri_sub(sampleNames(GSE45216_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE45216_raw_rma <- rma(GSE45216_raw)

#找到表达量，定义为最后需要的表达数据集
GSE45216 <- exprs(GSE45216_raw_rma) 

##将临床信息与表达量信息匹配起来
pheno_GSE45216 = pheno_GSE45216[c(1:30),] #去除光化病样本
GSE45216 = GSE45216[,rownames(pheno_GSE45216)]


##芯片注释
###anno the gene symbol
GSE45216_anno <- as.data.frame(GSE45216)
GSE45216_anno$ID<-rownames(GSE45216_anno)
gpl<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE45216_anno<-merge(x=GSE45216_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE45216_anno <- GSE45216_anno[,-1]
###The expression levels of the same gene names were averaged
GSE45216_anno<- aggregate(GSE45216_anno,by = list(GSE45216_anno$`Gene Symbol`),FUN = mean)
head(GSE45216_anno)

GSE45216_anno <- GSE45216_anno[-c(1),]##blank gene name was dropped
rownames(GSE45216_anno) <- GSE45216_anno$Group.1
GSE45216_anno <- GSE45216_anno[,-c(1,32)]
GSE45216_anno[1:2,1:2]
boxplot(GSE45216_anno)

matchSN = match(colnames(GSE45216_anno), rownames(pheno_GSE45216))
GSE45216_anno = GSE45216_anno[,matchSN]
save(file = "GSE45216.RData",GSE45216_anno,pheno_GSE45216)

#clean the datMet
pheno_GSE45216$sample_type = factor(pheno_GSE45216$sample_type,levels = c("N","T"))
pheno_GSE45216$Group = as.numeric(ifelse(pheno_GSE45216$sample_type=="T","1","0"))

#####
#####
###合并GSE45216和GSE67610
gene = intersect(rownames(GSE45216_anno),rownames(GSE98780_anno))
sel = c("Sample_title","ID","sample_type","organ","Group")
pheno_GSE98780_com = rbind(pheno_GSE45216[,sel],pheno_GSE98780[,sel])
GSE98780_com_anno = cbind(GSE98780_anno[gene,],GSE45216_anno[gene,])
boxplot(GSE98780_com_anno)


##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE98780_com_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE98780_com$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE98780_com_anno)[outliers]); print(table(outliers))
GSE98780_com_anno = GSE98780_com_anno[,!outliers]
pheno_GSE98780_com = pheno_GSE98780_com[!outliers,]


##回归未知协变量
pheno_GSE98780_com$sample_type = factor(pheno_GSE98780_com$sample_type,levels = c("N","T"))
pheno_GSE98780_com$Group = as.numeric(ifelse(pheno_GSE98780_com$sample_type=="T","1","0"))

mod = model.matrix(~as.factor(Group),data = pheno_GSE98780_com)
mod0 = model.matrix(~1,data = pheno_GSE98780_com)
n.sv = num.sv(GSE98780_com_anno, mod, method="be")
GSE98780_com_anno = as.matrix(GSE98780_com_anno)

svobj = sva(GSE98780_com_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE98780_com$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 62)

X = svobj$sv
Y = GSE98780_com_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE98780_com_anno = GSE98780_com_anno-t(to_regress)
boxplot(GSE98780_com_anno)
save(file = "GSE98780_com.RData",GSE98780_com_anno,pheno_GSE98780_com)









