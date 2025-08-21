#cervical_SCC_1b_GSE7803, Affymetrix
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio);

setwd("D:\\鳞癌分析\\泛鳞癌array数据\\cervix_SCC\\GSE7803_RAW")

dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE7803_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE7803_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE7803_raw)<-stri_sub(sampleNames(GSE7803_raw),1,9)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE7803_raw_rma <- rma(GSE7803_raw)

#找到表达量，定义为最后需要的表达数据集
GSE7803 <- exprs(GSE7803_raw_rma)


###anno the gene symbol
GSE7803_anno <- as.data.frame(GSE7803)
GSE7803_anno$ID<-rownames(GSE7803_anno)
gpl<-getGEO("GPL96",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE7803_anno<-merge(x=GSE7803_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE7803_anno <- GSE7803_anno[,-1]
###The expression levels of the same gene names were averaged
GSE7803_anno<- aggregate(GSE7803_anno,by = list(GSE7803_anno$`Gene Symbol`),FUN = mean)
head(GSE7803_anno)

GSE7803_anno <- GSE7803_anno[-c(1),]##blank gene name was dropped
rownames(GSE7803_anno) <- GSE7803_anno$Group.1
GSE7803_anno <- GSE7803_anno[,-c(1,43)]
GSE7803_anno[1:2,1:2]
boxplot(GSE7803_anno)


##match the metadata with the expression data
GSE7803_anno = GSE7803_anno[,rownames(pheno_GSE7803)]
matchSN = match(colnames(GSE7803_anno), rownames(pheno_GSE7803))
GSE7803_anno = GSE7803_anno[,matchSN]
save(file = "GSE7803.RData",GSE7803_anno,pheno_GSE7803)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE7803_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE7803$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE7803_anno)[outliers]); print(table(outliers))
GSE7803_anno = GSE7803_anno[,!outliers]
pheno_GSE7803 = pheno_GSE7803[!outliers,]

save(file = "GSE7803_remove.RData",pheno_GSE7803,GSE7803_anno)

#clean the datMet
pheno_GSE7803$sample_type = factor(pheno_GSE7803$sample_type,levels = c("N","T"))
pheno_GSE7803$Group = as.numeric(ifelse(pheno_GSE7803$sample_type=="T","1","0"))

###回归未知协变量
mod = model.matrix(~as.factor(Group),data = pheno_GSE7803)
mod0 = model.matrix(~1,data = pheno_GSE7803)
n.sv = num.sv(GSE7803_anno, mod, method="be")
GSE7803_anno = as.matrix(GSE7803_anno)

svobj = sva(GSE7803_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE7803$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 32)

X = svobj$sv
Y = GSE7803_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE7803_anno = GSE7803_anno-t(to_regress)
boxplot(GSE7803_anno)
save(file = "GSE7803.RData",pheno_GSE7803,GSE7803_anno)




