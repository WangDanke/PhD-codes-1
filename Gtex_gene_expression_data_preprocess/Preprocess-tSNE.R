###!/usr/bin/env Rscript
data = read.table("gtex_RSEM_gene_tpm", header = T)
pheno = read.table("GTEX_phenotype.txt", header = T, sep = "\t", quote = "")
rownames(data) = data[, 1]
data_new = data[, -1]
rownames(pheno) <- gsub('-','.', pheno[, 1])
sample = intersect(rownames(pheno),colnames(data_new))
str(pheno)
str(data_new)
pheno_Gtex <- pheno[sample,]
organs = pheno_Gtex[,3]
str(pheno_Gtex)
####PCA analysis
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("FactoMineR")) BiocManager::install("FactoMineR",update = F,ask = F)
library(FactoMineR)
ddb.pca <- PCA(t(data), graph = FALSE)
fviz_pca_ind(ddb.pca,
             geom.ind = "point",     # show points only (but not "text")
             col.ind = pheno$X_primary_site, # color by groups
             addEllipses = TRUE,     # Concentration ellipses
             legend.title = "groups",
             range=0
)

######t-SNE analysis
install.packages("Rtsne") # Rtsne包安装
library(Rtsne)
data.uni <- unique(data_new)##进行t-SNE分析时一定要去除重复,这里指的是不能有重复观测，也就是样本编号不能重复
set.seed(321) # 设置随机数种子
tsne_out = Rtsne(
  t(data_uni),
  dims = 2,
  pca = T,
  max_iter = 1000,
  theta = 0.4,
  perplexity = 20,
  verbose = F
) # 进行t-SNE降维分析
##保存数据
save(tsne_out, file = "tsne_out.RData")

###可视化
BiocManager::install("ggplot2")
library(ggplot2)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result)= c("tSNE1",'tSNE2')
ggplot(tsne_result,aes(tSNE1, tSNE2, color = organs)) +
geom_point() + 
theme_bw()
###ggsave 默认保存上一个图，只能保存ggplot做的图
ggsave(
  filename = "tsne-all.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 11,             # 宽
  height = 11,            # 高
  units = "in",          # 单位
  dpi = 300              # 分辨率DPI
)

######找到特定组织器官的数据做降维分析
pheno_sel = subset(pheno_Gtex, pheno_Gtex[,3]=="Bladder" | pheno_Gtex[,3]=="Esophagus" | pheno_Gtex[,3]=="Cervix Uteri" | pheno_Gtex[,3]=="Lung" | pheno_Gtex[,3]=="Skin")
organs_sel = pheno_sel[,2]##organs_sel = pheno_sel[,3]
data_sel = data.uni[,rownames(pheno_sel)]
str(data_sel)
str(pheno_sel)
tsne_out_sel = Rtsne(
  t(data_sel),
  dims = 2,
  pca = T,
  max_iter = 1000,
  theta = 0.4,
  perplexity = 20,
  verbose = F
)
###可视化
library(ggplot2)
set.seed(321)
tsne_result_sel = as.data.frame(tsne_out_sel$Y)
colnames(tsne_result_sel)= c("tSNE1",'tSNE2')
ggplot(tsne_result_sel,aes(tSNE1, tSNE2, color = organs_sel)) +
geom_point() + theme_bw()
ggsave(
  filename = "tsne_sel.png", # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7,             # 宽
  height = 7,            # 高
  units = "in",          # 单位
  dpi = 300              # 分辨率DPI
)


