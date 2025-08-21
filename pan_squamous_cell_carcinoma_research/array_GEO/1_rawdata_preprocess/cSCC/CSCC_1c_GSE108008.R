###CSCC_1c_GSE108008 #Affymetrix Human Gene 2.0 ST Array
library(stringi);library(pd.hugene.2.0.st);library(oligo)
 

setwd("D:\\鳞癌分析\\泛鳞癌array数据\\CSCC\\GSE108008_RAW")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名



GSE108008_raw<- read.celfiles( cel.files ) #读取文件
sampleNames(GSE108008_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE108008_raw)<-stri_sub(sampleNames(GSE108008_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE108008_raw_rma <- rma(GSE108008_raw)

#找到表达量，定义为最后需要的表达数据集
GSE108008 <- exprs(GSE108008_raw_rma) 

##去除不要的光华病样本
Ak = pheno_GSE108008$sample_type=="actinic keratosis"
pheno_GSE108008 = pheno_GSE108008[!Ak,]

##将临床信息与表达量信息匹配起来
GSE108008 = GSE108008[,rownames(pheno_GSE108008)]

save(file = "GSE108008.RData",GSE108008,pheno_GSE108008)

##芯片注释
###anno the gene symbol
GSE108008_anno <- as.data.frame(GSE108008)
GSE108008_anno$ID<-rownames(GSE108008_anno)

###加载从ThermorFisher 上下载的芯片注释文件
ANNO=fread(file = "HuGene-2_1-st-v1.na36.hg19.transcript.csv",sep = ",")
probe2gene <- ANNO[,c(2,8)]

library(stringr) 
probe2gene$symbol=trimws(str_split(probe2gene$gene_assignment,'//',simplify = T)[,2])
plot(table(table(probe2gene$symbol)),xlim=c(1,50))
head(probe2gene)

idname = probe2gene[,c(1,3)]
colnames(idname) = c("ID","Gene Symbol")

GSE108008_anno<-merge(x=GSE108008_anno,y=idname,by='ID',all.x=T,all.y=F)
GSE108008_anno <- GSE108008_anno[,-1]
###The expression levels of the same gene names were averaged
GSE108008_anno<- aggregate(GSE108008_anno,by = list(GSE108008_anno$`Gene Symbol`),FUN = mean)
head(GSE108008_anno)

GSE108008_anno <- GSE108008_anno[-c(1),]##blank gene name was dropped
rownames(GSE108008_anno) <- GSE108008_anno$Group.1
GSE108008_anno <- GSE108008_anno[,-c(1,22)]
GSE108008_anno[1:2,1:2]
boxplot(GSE108008_anno)

matchSN = match(colnames(GSE108008_anno), rownames(pheno_GSE108008))
GSE108008_anno = GSE108008_anno[,matchSN]
#保存文件
save(file = "GSE108008.RData",GSE108008_anno,pheno_GSE108008)

#clean the datMet
pheno_GSE108008$sample_type = factor(pheno_GSE108008$sample_type,levels = c("N","T"))
pheno_GSE108008$Group = as.numeric(ifelse(pheno_GSE108008$sample_type=="T","1","0"))



##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE108008_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE108008$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomright", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE108008_anno)[outliers]); print(table(outliers))
GSE108008_anno = GSE108008_anno[,!outliers]
pheno_GSE108008 = pheno_GSE108008[!outliers,]

####校正未知协变量
mod = model.matrix(~as.factor(Group),data = pheno_GSE108008)
mod0 = model.matrix(~1,data = pheno_GSE108008)
n.sv = num.sv(GSE108008_anno, mod, method="be")
GSE108008_anno = as.matrix(GSE108008_anno)

svobj = sva(GSE108008_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE108008$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 4)

X = svobj$sv
Y = GSE108008_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE108008_anno = GSE108008_anno-t(to_regress)
boxplot(GSE108008_anno)

save(file = "GSE108008.RData",GSE108008_anno,pheno_GSE108008)




