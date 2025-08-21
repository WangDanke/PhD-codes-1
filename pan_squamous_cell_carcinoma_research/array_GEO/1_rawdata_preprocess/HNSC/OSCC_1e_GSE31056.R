##OSCC_1e_GSE31056 
###Affymetrix 芯片数据
##选择工作路径
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\GSE31056")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE31056_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE31056_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE31056_raw)<-stri_sub(sampleNames(GSE31056_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE31056_raw_rma <- rma(GSE31056_raw)

#找到表达量，定义为最后需要的表达数据集
GSE31056 <- exprs(GSE31056_raw_rma) 
save(GSE31056,file = "GSE31056.RData")

##芯片注释

###anno the gene symbol
GSE31056_anno <- as.data.frame(GSE31056)
GSE31056_anno$ID<-rownames(GSE31056_anno)
library(hgu133plus2.db)
ids = toTable(hgu133plus2SYMBOL)##使用getgpl能注释上的基因只有两百个，于是果断放弃，使用hgu133plus2SYMBOL来注释
colnames(ids) = c("ID","Gene Symbol")
GSE31056_anno<-merge(x=GSE31056_anno,y=ids,by='ID',all.x=T,all.y=F)
GSE31056_anno <- GSE31056_anno[,-1]
###The expression levels of the same gene names were averaged
GSE31056_anno<- aggregate(GSE31056_anno,by = list(GSE31056_anno$`Gene Symbol`),FUN = mean)
head(GSE31056_anno)
##blank gene name was dropped
rownames(GSE31056_anno) <- GSE31056_anno$Group.1
GSE31056_anno <- GSE31056_anno[,-c(1,98)]
GSE31056_anno[1:2,1:2]
boxplot(GSE31056_anno)
save(GSE31056_anno,file = "GSE31056_anno.RData")


##match the metadata with the expression data
matchSN = match(colnames(GSE31056_anno), rownames(pheno_GSE31056))
GSE31056_anno = GSE31056_anno[,matchSN]
save(file = "GSE31056.RData",GSE31056_anno,pheno_GSE31056)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE31056_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE31056$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE31056_anno)[outliers]); print(table(outliers))
GSE31056_anno = GSE31056_anno[,!outliers]
pheno_GSE31056 = pheno_GSE31056[!outliers,]


save(file = "GSE31056_remove.RData",GSE31056_anno,pheno_GSE31056)

























