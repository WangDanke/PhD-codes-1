#BLCC_1a_GSE186691, [HTA-2_0] Affymetrix Human Transcriptome Array 2.0
library(stringi)
library(stringr)
library(FactoMineR)
library(WGCNA); library(lumi); library(ggplot2); library(nlme);library(affy)
library(ggplot2); library(Cairo); library(GEOquery);library(nlme);
library(biomaRt); library(sva);library(limma);library(illuminaio);library(oligo)


##选择工作路径
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\BLCC\\GSE186691_RAW")
dir = getwd()

##列出以.cel.gz结尾的文件
cel.files <- list.files(path = dir, pattern = ".+\\.cel.gz$", ignore.case = TRUE,
                        full.names = TRUE, recursive = TRUE)

basename(cel.files) #列出文件名

GSE186691_raw<- read.celfiles(filenames = cel.files) #读取文件
sampleNames(GSE186691_raw) #列出样本名
#the length of name maybe 8,9 or 10 #将样本名定义为前8，9,或10位
sampleNames(GSE186691_raw)<-stri_sub(sampleNames(GSE186691_raw),1,10)#the length of name maybe 8,9 or 10 #将样本名定义为前8位

##预处理原始数据，rma(Robust Multi-Array Average expression measure)：Background correcting，Normalizing，Calculating Expression
GSE186691_raw_rma <- rma(GSE186691_raw)

#找到表达量，定义为最后需要的表达数据集
GSE186691 <- exprs(GSE186691_raw_rma) 
datExpr = as.data.frame(GSE186691)






BiocManager::install("hta20transcriptcluster.db")
require(hta20transcriptcluster.db)
mapping <- mapIds(
  hta20transcriptcluster.db,
  keys = rownames(GSE186691_anno),
  column = 'SYMBOL',
  keytype = 'PROBEID')
















##芯片注释
library(biomaRt)

##
rm(list = ls())
options(stringsAsFactors = F)

#加载R包
library(GEOquery)

##读入下载好的文件
gse <- getGEO(filename = "GSE186691_family.soft.gz",destdir = ".")
str(gse)
length(gse)

##提取探针
id_probe <- gse@gpls$GPL17586@dataTable@table
dim(id_probe)
head(id_probe)
id_probe[1:4,1:15]
View(head(id_probe))

probe2gene <- id_probe[,c(2,8)]

library(stringr) 
probe2gene$symbol=trimws(str_split(probe2gene$gene_assignment,'//',simplify = T)[,2])
plot(table(table(probe2gene$symbol)),xlim=c(1,50))
head(probe2gene)

dim(probe2gene)
View(head(probe2gene))
ids2 <- probe2gene[,c(1,3)]
View(head(ids2))
ids2[1:20,1:2]#含有缺失值
table(table(unique(ids2$symbol)))#30907 ,30906个基因，一个空字符
save(ids2,probe2gene,file='gse-probe2gene.Rdata')



library(biomaRt)

# 加载biomaRt包
library(biomaRt)

value <- probe2gene$probeset_id
attr <- c("affy_hta_2_0","hgnc_symbol")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ids <- getBM(attributes = attr,
             filters = "affy_hta_2_0",
             values = value,
             mart = ensembl,
             useCache = F)

dim(ids)#[1] 1041 2
View(head(ids))

save(ids,file = "GPL17586_ids.Rdata")

#去重之后

attributes <- listAttributes(ensembl)
View(attributes) # 查看支持的芯片转换格式

save(ids,ensembl,y,file = "ensembl.Rdata")

plot(table(table(ids$hgnc_symbol)),xlim=c(1,50))

table(table(unique(ids$hgnc_symbol)))#去重之后有29262，丢失了一很多

# 去重复
ids3 <- ids[!duplicated(ids$hgnc_symbol),]





##Re-annotate Probes
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
identifier <- "affy_hg_u133_plus_2"	
efg_affy_hg_u133_plus_2_bool
getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id","external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position")
geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
idx = match(rownames(datExpr), geneDat[,"affy_hg_u133_plus_2"])
datProbes = cbind(rownames(datExpr), geneDat[idx,])
datExpr$gene_name = datProbes[,"Associated Gene Name"]
###The expression levels of the same gene names were averaged
datExpr<- aggregate(datExpr,by = list(datExpr$`gene_name`),FUN = mean)
rownames(datExpr) = datExpr$Group.1
datExpr = datExpr[,-c(1,6)]

##保存数据
GSE186691 = datExpr





