#OSCC_1b_GSE64216 8samples,##illumina 芯片数据

rm(list=ls());
options(stringsAsFactors=F)

suppressPackageStartupMessages(T)


##load rawdata
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\OSCC\\GSE64216")
getwd()
#read raw data
data.lumi = lumiR("GSE64216_non-normalized_data.txt")
datMeta = pheno_GSE64216
matchSN = match(sampleNames(data.lumi), datMeta$ID)
datMeta = datMeta[matchSN,]

# log2 transform
dataAll.lumi<-lumiT(data.lumi, method="log2"); 
dataSCC.lumi<-dataAll.lumi[,datMeta$sample_type!="OSF"]; 
dataNOSCC.lumi<-dataAll.lumi[,datMeta$sample_type=="OSF"]

#Normalize
dataAll_N.lumi<-lumiN(dataAll.lumi, method="quantile"); 
dataSCC_N.lumi<- lumiN(dataSCC.lumi, method="quantile"); 
dataNOSCCL_N.lumi<- lumiN(dataNOSCC.lumi, method="quantile");

#Extract expression data for normal tissue and SCC tumor sample
datExpr = as.data.frame(exprs(dataSCC_N.lumi))
datMeta = datMeta[datMeta$sample_type!="OSF",]
datExpr.prenorm = exprs(dataSCC.lumi)

##Re-annotate Probes
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
identifier <- "illumina_humanht_12_v4"
getinfo <- c("illumina_humanht_12_v4", "ensembl_gene_id","external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position")
geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
idx = match(rownames(datExpr), geneDat[,"Illumina Human HT 12 V4 probe"])
datProbes = cbind(rownames(datExpr), geneDat[idx,])
datExpr$gene_name = datProbes[,"Associated Gene Name"]
###The expression levels of the same gene names were averaged
datExpr<- aggregate(datExpr,by = list(datExpr$`gene_name`),FUN = mean)
rownames(datExpr) = datExpr$Group.1
datExpr = datExpr[,-c(1,6)]

##保存数据
pheno_GSE64216 = datMeta
GSE64216_anno = datExpr
save(file = "GSE64216.RData",GSE64216_anno,pheno_GSE64216)


##移除离群值
sdout <- 2; normadj <- (0.5+0.5*bicor(GSE64216_anno))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
plot(z.ku, col = as.numeric(as.factor(pheno_GSE64216$sample_type)), pch=19, main="Outliers", ylab="Connectivity (z)")
abline(h=-2, lty=2)
legend("bottomleft", legend=c("N", "T"), col = 1:2, pch=19)
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(GSE64216_anno)[outliers]); print(table(outliers))
GSE64216_anno = GSE64216_anno[,!outliers]
pheno_GSE64216 = pheno_GSE64216[!outliers,]


##保存数据
save(file = "GSE64216_remove.RData",GSE64216_anno,pheno_GSE64216)







