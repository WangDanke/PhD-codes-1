setwd("D:\\鳞癌分析\\泛鳞癌array数据\\12_Subtype_analysis\\")
load("D:\\鳞癌分析\\泛鳞癌array数据\\12_Subtype_analysis\\SCC_combined_combat_for_WGCNA.RData")
load("D:\\鳞癌分析\\泛鳞癌array数据\\12_Subtype_analysis\\TCGA_4_SCCs.RData")
###使用TCGA的数据来验证GEO芯片数据得到的亚型，并且，确定不同的亚型是否有不同的预后
TCGA_SCC_EXPR <- TCGA_SCC_EXPR[!duplicated(TCGA_SCC_EXPR$gene_name, fromLast=T), ]
rownames(TCGA_SCC_EXPR) = TCGA_SCC_EXPR$gene_name
TCGA_SCC_EXPR = TCGA_SCC_EXPR[rownames(SCC_EXPR),]
save(TCGA_SCC_EXPR, pheno_SCC_TCGA,file = "TCGA_4_SCCs_8816x1169.RData")

#####使用WGCNA得到的hub gene,top10%或5% in each module,
hub_gene = read.csv("9_module_hub_gene_10%.csv",header = T)
##hub_gene = read.csv("9_module_hub_gene_5%.csv",header = T)使用top5%的样本
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(dplyr)
library(Hmisc)

##########使用与hub gene基因鉴别分型,两种情况（一种包含LUSC，一种不包含LUSC）,同时分别使用top10%和5%的hub gene

library("ConsensusClusterPlus")
######
##聚类鉴别分型_TCGA
TCGA_SCC_EXPR = TCGA_SCC_EXPR[,-c(1,2,3)]
SCC_TCGA_hub = TCGA_SCC_EXPR[hub_gene$gene,]
SCC_TCGA_hub = log2(SCC_TCGA_hub+1)
pheno_SCC_TCGA = subset(pheno_SCC_TCGA,!pheno_SCC_TCGA$project_id=="TCGA-LUSC")##去除LUSC样本
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(pheno_SCC_TCGA)]

##归一化操作
SCC_TCGA_hub = sweep(SCC_TCGA_hub,1, apply(SCC_TCGA_hub,1,median,na.rm=T))
SCC_TCGA_hub = as.matrix(SCC_TCGA_hub)
SCC_TCGA_hub = SCC_TCGA_hub[hub_gene$gene,]
##鉴别分型
title="maxK_10_TCGA_3SSC_hub_10"
Cluster_TCGA = ConsensusClusterPlus(SCC_TCGA_hub,maxK=10,reps=50,pItem=0.8,pFeature=1,
                                    title = title,clusterAlg="hc",distance="pearson",seed=1262118,plot="png")


##验证分型的预后
library(survival)
library(survminer)
library(preprocessCore)

BB = as.data.frame(Cluster_TCGA[[4]][["consensusClass"]])
pheno_SCC_TCGA$Subtype_3SCC_10= BB$`Cluster_TCGA[[4]][["consensusClass"]]`

pheno_SCC_TCGA$time=as.numeric(pheno_SCC_TCGA$days_to_death)+as.numeric(pheno_SCC_TCGA$days_to_last_follow_up)
pheno_SCC_TCGA$time = pheno_SCC_TCGA$time/365
pheno_SCC_TCGA$event=ifelse(pheno_SCC_TCGA$vital_status=='Alive',0,1)

sfit <- survfit(Surv(time, event)~Subtype_3SCC_10, data=pheno_SCC_TCGA)
print(sfit)
colors <- c("#C6B3D3","#ED9F9B","#80BA8A","#9CD1C8")  # 对应的颜色
labels <- c("1","2","3","4")  
ggsurvplot(sfit, 
           conf.int = FALSE, 
           pval = TRUE, 
           palette = colors,    # 设置颜色
           legend.labs = labels # 设置图例标签
)

p

######绘制hub gene 的分组热图
library(dplyr)
annotation_col = pheno_SCC_TCGA[,c("gender",'project_id',"ajcc_clinical_t","Subtype_3SCC_10")]
write.csv(annotation_col,file = "annotation.csv")
##保存以后，在Excel里将分期不明的样本填上NA，并按照亚型排序。
annotation_col = read.csv("annotation.csv",header = T,row.names = 1)
annotation_row =  read.csv("hub_gene_100.csv",header = T,row.names = 1)
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(annotation_col)]
annotation_col<- annotation_col %>%
  mutate(Subtype_3SCC_10 = case_when(
    Subtype_3SCC_10 == 1 ~ "subtype1",
    Subtype_3SCC_10 == 2 ~ "subtype2",
    Subtype_3SCC_10 == 3 ~ "subtype3",
    Subtype_3SCC_10 == 4 ~ "subtype4"
  ))
annotation_colors = list(
  module = c("CD1"= "turquoise","CD2"="blue","CD3"="brown","CD4"="yellow","CD5"="green","CD6"="red",
             "CD7" = "black","CD8"="pink","CD9"="magenta"),
  gender = c("female"="#FFC1C1","male"="#87CEEB"),
  project_id = c("TCGA-CESC"="#778899","TCGA-ESCA" = "#CD5555","TCGA-HNSC"="#CDAA7D"),
  ajcc_clinical_t = c("NA" = "#E8E8E8","T1"="#CFCFCF","T2" = "#B5B5B5","T2" = "#9C9C9C","T3" = "#828282",
                      "T4" = "#696969","T4a" = "#4F4F4F","T4b"="#363636"),
  Subtype_3SCC_10 = c("subtype1"="#C6B3D3","subtype2"="#ED9F9B","subtype3"="#80BA8A","subtype4"="#9CD1C8")
  
)


###

library(pheatmap)
pheatmap(SCC_TCGA_hub ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =T,
         annotation_col = annotation_col,annotation_row = annotation_row,annotation_colors = annotation_colors)

###########
###########
########TCGA SCCs亚型的分组描述
library(compareGroups)
annotation_col$gender = as.factor(annotation_col$gender)
annotation_col$project_id = as.factor(annotation_col$project_id)
annotation_col$ajcc_clinical_t = as.factor(annotation_col$ajcc_clinical_t)
annotation_col$Subtype_3SCC_10 = as.factor(annotation_col$Subtype_3SCC_10)

pheno_SCC_TCGA = pheno_SCC_TCGA[rownames(annotation_col),]
annotation_col$age = pheno_SCC_TCGA$age_at_index
annotation_col$age = as.numeric(annotation_col$age)
save(annotation_col,file = "annotation_col.RData")

annotation_col = annotation_col[,c(5,1,2,3,4)]
res1 <- compareGroups(Subtype_3SCC_10 ~ ., data = annotation_col, ref = 1)
restab = createTable(res1, show.ratio = TRUE)
restab
export2word(restab, file='table1.docx')



####
##聚类鉴别分型_GEO，验证分型
pheno_SCC_T = subset(pheno_SCC,pheno_SCC$sample_type=="T")
pheno_SCC_4 = subset(pheno_SCC_T,!pheno_SCC_T$organ=="lung")
pheno_SCC_3 = subset(pheno_SCC_4,!pheno_SCC_4$organ=="skin")
SCC_EXPR_T_3 = SCC_EXPR[,rownames(pheno_SCC_3)]

SCC_GEO_hub = SCC_EXPR_T_3[hub_gene$gene,]

##归一化操作
SCC_GEO_hub = sweep(SCC_GEO_hub,1, apply(SCC_GEO_hub,1,median,na.rm=T))
SCC_GEO_hub = as.matrix(SCC_GEO_hub)
title="maxK_10_GEO_3SCC_hub_10"
Cluster_GEO = ConsensusClusterPlus(SCC_GEO_hub,maxK=10,reps=50,pItem=0.8,pFeature=1,
                               title = title,clusterAlg="hc",distance="pearson",seed=1262118,plot="png")
###将分型信息与临床信息匹配上
AA = as.data.frame(Cluster_GEO[[4]][["consensusClass"]])
pheno_SCC_3$Subtype = AA$`Cluster_GEO[[4]][["consensusClass"]]`



pheno_SCC_3<-pheno_SCC_3 %>%
  mutate(Subtype = case_when(
    Subtype == 1 ~ "subtype1",
    Subtype == 2 ~ "subtype2",
    Subtype == 3 ~ "subtype3",
    Subtype == 4 ~ "subtype4"
  ))
pheno_SCC_3 <- pheno_SCC_3[order(pheno_SCC_3$Subtype), ]
annotation_col = pheno_SCC_3[,c(4,5,7)]

annotation_col<-annotation_col %>%
  mutate(Subtype = case_when(
    Subtype == "subtype1" ~ "subtype2",
    Subtype == "subtype2" ~ "subtype1",
    Subtype == "subtype3" ~ "subtype3",
    Subtype == "subtype4" ~ "subtype4"
  ))

annotation_col <- annotation_col[order(annotation_col$Subtype), ]


annotation_row =  read.csv("9_module_hub_gene_10%.csv",header = T,row.names = 1)

SCC_GEO_hub = SCC_GEO_hub[,rownames(annotation_col)]

annotation_colors = list(
  module = c("CD1"= "turquoise","CD2"="blue","CD3"="brown","CD4"="yellow","CD5"="green","CD6"="red",
             "CD7" = "black","CD8"="pink","CD9"="magenta"),
  gender = c("female"="#FFC1C1","male"="#87CEEB"),
  project_id = c("TCGA-CESC"="#778899","TCGA-ESCA" = "#CD5555","TCGA-HNSC"="#CDAA7D"),
  ajcc_clinical_t = c("NA" = "#E8E8E8","T1"="#CFCFCF","T2" = "#B5B5B5","T2" = "#9C9C9C","T3" = "#828282",
                      "T4" = "#696969","T4a" = "#4F4F4F","T4b"="#363636"),
  Subtype_3SCC_10 = c("subtype1"="#C6B3D3","subtype2"="#ED9F9B","subtype3"="#80BA8A","subtype4"="#9CD1C8")
  
)

library(pheatmap)
pheatmap(SCC_GEO_hub ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =T,
         annotation_col = annotation_col,annotation_row = annotation_row)


save(pheno_SCC_3,SCC_EXPR_T_3,file = "GEO_4subtypes_without_CSCC_LUSC.RData")
save(pheno_SCC_3,SCC_GEO_hub,file = "GEO_4subtypes_3_SCC.RData")


#####
library(Hmisc)
SCC_samples_cor <- cor(SCC_GEO_hub,SCC_TCGA_hub,method = "spearman")

GEO_subtype2 <- subset(pheno_SCC_3,pheno_SCC_3$subtype=="subtype2")
SCC_cor_subtype2_GEO <- SCC_samples_cor[rownames(GEO_subtype2),]

cor_means <- colMeans(SCC_cor_subtype2_GEO)
SCC_cor_subtype2_GEO = as.data.frame(t(SCC_cor_subtype2_GEO))
SCC_cor_subtype2_GEO$mean = cor_means
SCC_cor_subtype2_GEO = SCC_cor_subtype2_GEO[rownames(annotation_col),]
SCC_cor_subtype2_GEO$group <- annotation_col$Subtype_3SCC_10

library(ggplot2)
library(ggpubr)
my_comparisons <- list(c('subtype1','subtype2'),
                       c('subtype1','subtype3'),
                       c('subtype1','subtype4'),
                       c('subtype2','subtype3'),
                       c('subtype2','subtype4'),
                       c('subtype3','subtype4'))
ggviolin(SCC_cor_subtype2_GEO,x="group",y="mean",fill='group', 
         palette = c("#C6B3D3","#ED9F9B","#80BA8A","#9CD1C8"),
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = my_comparisons)


















SCC_samples_cor <- cor(SCC_GEO_hub,SCC_TCGA_hub,method = "spearman")
TCGA_subtype1 <- subset(annotation_col,annotation_col$Subtype_3SCC_10=="subtype1")
ESCC_cor_subtype1_TCGA <- SCC_samples_cor[,rownames(TCGA_subtype1)]
ESCC_cor_subtype1_TCGA = ESCC_cor_subtype1_TCGA[rownames(subtype_1_GEO),]

ESCC_cor_subtype1_TCGA$mean_cor <- apply(ESCC_cor_subtype1_TCGA,1,mean)
ESCC_cor_subtype1_TCGA$group <- GSE_pheno$SUBTYPE

save(ESCC_cor_subtype1_TCGA,ESCC_cor_subtype2_TCGA,ESCC_cor_subtype3_TCGA,file = "TCGA_Subtypes_cor_GEO_Subtypes.RData")
###做小提琴图
BiocManager::install("vioplot")
library("vioplot")
vioplot(mean_cor~group, data = ESCC_cor_subtype1_TCGA, 
        main = "Correlations of TCGA_Subtype1 with Subtypes", # 设置标题
        col=c("#999999", "#E69F00", "#56B4E9"))

##另一种可以比较组间差异的做小提琴的方法
library(ggpubr)
my_comparisons <- list(c('subtype1','subtype2'),
                       c('subtype2','subtype3'),c('subtype1','subtype3'))
ggviolin(ESCC_cor_subtype3_TCGA,x="group",y="mean_cor",fill='group', 
         palette = c("#999999", "#E69F00", "#56B4E9"),
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = my_comparisons)



row_standardize <- function(df) {
  scale(df)
}


SCC_TCGA_hub_N = scale(SCC_TCGA_hub)
SCC_GEO_hub_N = row_standardize(SCC_GEO_hub)

SCC_TCGA_hub_N = SCC_TCGA_hub_N[,rownames(annotation_col)]


pheno_SCC_3 <- pheno_SCC_3[order(pheno_SCC_3$subtype), ]
SCC_GEO_hub_N = SCC_GEO_hub_N[,rownames(pheno_SCC_3)]


pheatmap(SCC_TCGA_hub_N ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =F,
         annotation_col = annotation_col,annotation_row = annotation_row,annotation_colors = annotation_colors)


pheatmap(SCC_TCGA_hub_N ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =F,
         annotation_col = annotation_col,annotation_row = annotation_row,annotation_colors = annotation_colors)







                     