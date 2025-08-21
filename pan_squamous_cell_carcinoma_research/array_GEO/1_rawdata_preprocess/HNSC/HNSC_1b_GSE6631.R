library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio)

###HNSC GSE6631 44samples
##选择工作路径
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE6631")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE6631_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE6631_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE6631_raw)<-stri_sub(sampleNames(GSE6631_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE6631_raw_rma <- rma(GSE6631_raw)

#找到表达量，定义为最后需要的表达数据集
GSE6631 <- exprs(GSE6631_raw_rma) 


##芯片注释
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\HNSC\\GSE6631")
###anno the gene symbol
GSE6631_anno <- as.data.frame(GSE6631)
GSE6631_anno$ID<-rownames(GSE6631_anno)
gpl571<-getGEO("GPL15207",destdir=".")
iddata<-as.data.frame(Table(gpl571)[c('ID','Gene Symbol')])
GSE6631_anno<-merge(x=GSE6631_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE6631_anno <- GSE6631_anno[,-1]
###The expression levels of the same gene names were averaged
GSE6631_anno<- aggregate(GSE6631_anno,by = list(GSE6631_anno$`Gene Symbol`),FUN = mean)
head(GSE6631_anno)

GSE6631_anno <- GSE6631_anno[-c(1,2),]##blank gene name was dropped
rownames(GSE6631_anno) <- GSE6631_anno$Group.1
GSE6631_anno <- GSE6631_anno[,-c(1,49)]
GSE6631_anno[1:2,1:2]

##match the metadata with the expression data
matchSN = match(colnames(GSE6631_anno), rownames(pheno_GSE6631))
GSE6631_anno = GSE6631_anno[,matchSN]
save(file = "GSE6631.RData",GSE6631_anno,pheno_GSE6631)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE6631_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE6631$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE6631_anno)[outliers]); print(table(outliers))
GSE6631_anno = GSE6631_anno[,!outliers]
pheno_GSE6631 = pheno_GSE6631[!outliers,]

#clean the datMet
pheno_GSE6631$sample_type = factor(pheno_GSE6631$sample_type,levels = c("N","T"))
pheno_GSE6631$Group = as.numeric(ifelse(pheno_GSE6631$sample_type=="T","1","0"))



####回归除了诊断组（tumor normal）之外所有未知协变量的影响 
mod = model.matrix(~as.factor(Group),data = pheno_GSE6631)
mod0 = model.matrix(~1,data = pheno_GSE6631)
n.sv = num.sv(GSE6631_anno, mod, method="be")
GSE6631_anno = as.matrix(GSE6631_anno)

svobj = sva(GSE6631_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE6631$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 42)

X = svobj$sv
Y = GSE6631_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE6631_anno = GSE6631_anno-t(to_regress)

boxplot(GSE6631_anno)
save(file = "GSE6631.RData",GSE6631_anno,pheno_GSE6631)

















