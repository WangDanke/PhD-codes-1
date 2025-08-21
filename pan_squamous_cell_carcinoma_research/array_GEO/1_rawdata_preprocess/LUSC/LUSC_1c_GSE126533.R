#LUSC_1c_GSE126533, Agilent
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio);


raw_dir = "D:\\鳞癌分析\\泛鳞癌array数据\\LUSC\\GSE126533"
raw_datas = paste0(raw_dir,"/",dir(raw_dir))
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\LUSC\\GSE126533")



library(affy)
# 获取所有TXT文件的文件名
txt_files <- list.files("D:\\鳞癌分析\\泛鳞癌array数据\\LUSC\\GSE126533", pattern = "\\.txt.gz", full.names = TRUE)

setwd("D:\\鳞癌分析\\泛鳞癌array数据\\LUSC\\GSE126533")
GSE126533_raw <- read.maimages(txt_files,
                               source="agilent", 
                               green.only=TRUE,
                               other.columns = "gIsWellAboveBG")
dim(GSE126533_raw)

##背景校正
GSE126533_1 <- backgroundCorrect(GSE126533_raw,method = "normexp") 
##标准化(quantile normalization)
GSE126533_1 <- normalizeBetweenArrays(GSE126533_1, method="quantile")
class(GSE126533_1)

##基因过滤
##找到对照探针
Control <- GSE126533_1$genes$ControlType==1L
table(Control)

##找到探针名是空白
NoSymbol <- is.na(GSE126533_1$genes$ProbeName)
table(NoSymbol)

##找到低表达量的探针，超过一半数量的样本的表达量大于背景噪音，通过
IsExpr <- rowSums(GSE126533_1$other$gIsWellAboveBG > 0) >= 5
table(IsExpr)

#最后保留的是，非对照探针，且探针名不空白，且超过一半数量的样本的表达量大于背景噪音
GSE126533_1_filt <- GSE126533_1[!Control & !NoSymbol & IsExpr, ]
dim(GSE126533_1_filt)

#最后，得到表达矩阵
GSE126533 = GSE126533_1_filt@.Data[[1]]
boxplot(GSE126533) #看一下基本情况


#给列名命名为样本编号
colnames(GSE126533) = str_extract(colnames(GSE126533),"GSM\\d*")
GSE126533[1:2,1:2]
#给行名命名为探针编号
rownames(GSE126533)= GSE126533_1_filt$genes$ProbeName

##芯片注释
###anno the gene symbol
GSE126533_anno <- as.data.frame(GSE126533)
GSE126533_anno$ID<-rownames(GSE126533_anno)

#gpl<- getGEO(filename = "GPL14550"),找不到，直接在网页上下载了相关的注释文件
gpl<-read.csv("GPL22120-25936.csv")
gpl = gpl[-c(1:32),]
colnames(gpl) = gpl[1,]
gpl = gpl[-1,]
iddata = gpl[,c("ID","GENE_SYMBOL")]

GSE126533_anno<-merge(x=GSE126533_anno,y=iddata,by='ID',all.x=T,all.y=F)
GSE126533_anno <- GSE126533_anno[,-1]
###The expression levels of the same gene names were averaged
GSE126533_anno<- aggregate(GSE126533_anno,by = list(GSE126533_anno$`GENE_SYMBOL`),FUN = mean)
head(GSE126533_anno)

GSE126533_anno <- GSE126533_anno[-c(1),]##blank gene name was dropped
rownames(GSE126533_anno) <- GSE126533_anno$Group.1
GSE126533_anno <- GSE126533_anno[,-c(1,12)]
GSE126533_anno[1:2,1:2]
boxplot(GSE126533_anno)

save(GSE126533_anno,file = "GSE126533_anno.RData")

##match the metadata with the expression data
#clean the datMet
pheno_GSE126533$sample_type = factor(pheno_GSE126533$sample_type,levels = c("N","T"))
pheno_GSE126533$Group = as.numeric(ifelse(pheno_GSE126533$sample_type=="T","1","0"))

matchSN = match(colnames(GSE126533_anno), rownames(pheno_GSE126533))
GSE126533_anno = GSE126533_anno[,matchSN]
save(file = "GSE126533.RData",GSE126533_anno,pheno_GSE126533)


##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE126533_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE126533$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE126533_anno)[outliers]); print(table(outliers))
GSE126533_anno = GSE126533_anno[,!outliers]
pheno_GSE126533 = pheno_GSE126533[!outliers,]

###校正所有未知协变量
mod = model.matrix(~as.factor(Group),data = pheno_GSE126533)
mod0 = model.matrix(~1,data = pheno_GSE126533)
n.sv = num.sv(GSE126533_anno, mod, method="be")
GSE126533_anno = as.matrix(GSE126533_anno)

svobj = sva(GSE126533_anno,mod, mod0,n.sv=n.sv)
cor = cor(pheno_GSE126533$Group, svobj$sv, method = "spearman")
p = corPvalueFisher(cor,nSamples = 9)

X = svobj$sv
Y = GSE126533_anno
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X) %*% (as.matrix(beta)))
GSE126533_anno = GSE126533_anno-t(to_regress)
boxplot(GSE126533_anno) 

save(file = "GSE126533.RData",GSE126533_anno,pheno_GSE126533)











