#OSCC_1g_GSE138206, Affymetrix


setwd("D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\GSE138206")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名
GSE138206_raw<- ReadAffy(filenames = cel.files) #读取文件
sampleNames(GSE138206_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE138206_raw)<-stri_sub(sampleNames(GSE138206_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE138206_raw_rma <- rma(GSE138206_raw)

#找到表达量，定义为最后需要的表达数据集
GSE138206 <- exprs(GSE138206_raw_rma) 
save(GSE138206,file = "GSE138206.RData")

##芯片注释

###anno the gene symbol
GSE138206_anno <- as.data.frame(GSE138206)
GSE138206_anno$ID<-rownames(GSE138206_anno)
gpl<-getGEO("GPL570",destdir=".")
iddata<-as.data.frame(Table(gpl)[c('ID','Gene Symbol')])
GSE138206_anno<-merge(x=GSE138206_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE138206_anno <- GSE138206_anno[,-1]
###The expression levels of the same gene names were averaged
GSE138206_anno<- aggregate(GSE138206_anno,by = list(GSE138206_anno$`Gene Symbol`),FUN = mean)
head(GSE138206_anno)

GSE138206_anno <- GSE138206_anno[-c(1),]##blank gene name was dropped
rownames(GSE138206_anno) <- GSE138206_anno$Group.1
GSE138206_anno <- GSE138206_anno[,-c(1,20)]
GSE138206_anno[1:2,1:2]
boxplot(GSE138206_anno)
save(GSE138206_anno,file = "GSE138206_anno.RData")


##match the metadata with the expression data
matchSN = match(colnames(GSE138206_anno), rownames(pheno_GSE138206))
GSE138206_anno = GSE138206_anno[,matchSN]
save(file = "GSE138206.RData",GSE138206_anno,pheno_GSE138206)

##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE138206_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE138206$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE138206_anno)[outliers]); print(table(outliers))
GSE138206_anno = GSE138206_anno[,!outliers]
pheno_GSE138206 = pheno_GSE138206[!outliers,]


save(file = "GSE138206_remove.RData",GSE138206_anno,pheno_GSE138206)



















