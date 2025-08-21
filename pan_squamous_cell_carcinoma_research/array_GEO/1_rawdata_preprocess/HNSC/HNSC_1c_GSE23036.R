library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)

###HNSC GSE23036 68samples
##选择工作路径
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE23036")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE23036_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE23036_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE23036_raw)<-stri_sub(sampleNames(GSE23036_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE23036_raw_rma <- rma(GSE23036_raw)

#找到表达量，定义为最后需要的表达数据集
GSE23036 <- exprs(GSE23036_raw_rma) 


##芯片注释
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE23036")
###anno the gene symbol
GSE23036_anno <- as.data.frame(GSE23036)
GSE23036_anno$ID<-rownames(GSE23036_anno)
gpl571<-getGEO("GPL15207",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE23036_anno<-merge(x=GSE23036_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE23036_anno <- GSE23036_anno[,-1]
###The expression levels of the same gene names were averaged
GSE23036_anno<- aggregate(GSE23036_anno,by = list(GSE23036_anno$`Gene Symbol`),FUN = mean)
head(GSE23036_anno)

GSE23036_anno <- GSE23036_anno[-c(1,2),]##blank gene name was dropped
rownames(GSE23036_anno) <- GSE23036_anno$Group.1
GSE23036_anno <- GSE23036_anno[,-c(1,49)]
GSE23036_anno[1:2,1:2]

##match the metadata with the expression data
matchSN = match(colnames(GSE23036_anno), rownames(pheno_GSE23036))
GSE23036_anno = GSE23036_anno[,matchSN]
save(file = "GSE23036.RData",GSE23036_anno,pheno_GSE23036)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE23036_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE23036$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE23036_anno)[outliers]); print(table(outliers))
GSE23036_anno = GSE23036_anno[,!outliers]
pheno_GSE23036 = pheno_GSE23036[!outliers,]

#clean the datMet
pheno_GSE23036$sample_type = factor(pheno_GSE23036$sample_type,levels = c("N","T"))
pheno_GSE23036$Group = as.numeric(ifelse(pheno_GSE23036$sample_type=="T","1","0"))



####回归除了诊断组（tumor normal）之外所有未知协变量的影响 
mod = model.matrix(~as.factor(Group),data = pheno_GSE23036)
mod0 = model.matrix(~1,data = pheno_GSE23036)
n.sv = num.sv(GSE23036_anno, mod, method="be")
GSE23036_anno = as.matrix(GSE23036_anno)

svobj = sva(GSE23036_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE23036$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 65)

X = svobj$sv
Y = GSE23036_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE23036_anno = GSE23036_anno-t(to_regress)
boxplot(GSE23036_anno)
save(file = "GSE23036.RData",GSE23036_anno,pheno_GSE23036)




