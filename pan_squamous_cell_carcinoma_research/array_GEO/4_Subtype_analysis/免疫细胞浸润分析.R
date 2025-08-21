######TCGA 三种鳞癌免疫细胞浸润分析，
setwd("D:/鳞癌分析/泛鳞癌array数据/12_Subtype_analysis/")
load("D:/鳞癌分析/泛鳞癌array数据/12_Subtype_analysis/TCGA_SCC_EXPR_8816x1168.RData")
TCGA_SCC_EXPR = TCGA_SCC_EXPR[,-c(1,2,3)]
TCGA_SCC_EXPR = log2(TCGA_SCC_EXPR+1)

pheno_SCC_TCGA = subset(pheno_SCC_TCGA,!pheno_SCC_TCGA$project_id=="TCGA-LUSC")##去除LUSC样本
SCC_TCGA = TCGA_SCC_EXPR[,rownames(pheno_SCC_TCGA)]

library(dplyr)
annotation_col = read.csv("annotation.csv",header = T,row.names = 1)
annotation_col<- annotation_col %>%
  mutate(Subtype_3SCC_10 = case_when(
    Subtype_3SCC_10 == 1 ~ "subtype1",
    Subtype_3SCC_10 == 2 ~ "subtype2",
    Subtype_3SCC_10 == 3 ~ "subtype3",
    Subtype_3SCC_10 == 4 ~ "subtype4"
  ))
save(annotation_col,file = "annotation_col.RData")

SCC_TCGA = SCC_TCGA[,rownames(annotation_col)]

gene_set<- read.csv('D:\\鳞癌分析\\泛鳞癌array数据\\11_ssGSEA\\gene_reference.csv',
                    header = T)
gene_set<-gene_set[, 1:2]
head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
##ssGSEA
library(GSVA)
data_1_gsva<- gsva(as.matrix(SCC_TCGA), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

data_1_gsva1<- t(scale(t(data_1_gsva)))#
data_1_gsva1[data_1_gsva1< -2] <- -2
data_1_gsva1[data_1_gsva1>2] <- 2
anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
anti<- gsub('^ ','',rownames(data_1_gsva1))%in%anti_tumor
pro<- gsub('^ ','',rownames(data_1_gsva1))%in%pro_tumor
non <- !(anti|pro)##
data_1_gsva1<- rbind(data_1_gsva1[anti,],data_1_gsva1[pro,],data_1_gsva1[non,])#
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}# set the normalization function
nor_data_1_gsva1 <- normalization(data_1_gsva1)

##绘制免疫细胞浸润热图
library(pheatmap)
annotation_colors = list(
  gender = c("female"="#FFC1C1","male"="#87CEEB"),
  project_id = c("TCGA-CESC"="#778899","TCGA-ESCA" = "#CD5555","TCGA-HNSC"="#CDAA7D"),
  ajcc_clinical_t = c("NA" = "#E8E8E8","T1"="#CFCFCF","T2" = "#B5B5B5","T2" = "#9C9C9C","T3" = "#828282",
                      "T4" = "#696969","T4a" = "#4F4F4F","T4b"="#363636"),
  Subtype_3SCC_10 = c("subtype1"="#C6B3D3","subtype2"="#ED9F9B","subtype3"="#80BA8A","subtype4"="#9CD1C8")
  
)

pheatmap(nor_data_1_gsva1 ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =T,
         annotation_col = annotation_col,annotation_colors = annotation_colors)

save(nor_data_1_gsva1,file = 'TCGA_3SCC_nor_immune_infiltration_score.Rdata')

##绘制具体的细胞类型在四个亚型中浸润程度的箱线图
data = nor_data_1_gsva1
data <- as.data.frame(data)
data <- as.data.frame(t(data))
data$SUBTYPE<-annotation_col$Subtype_3SCC_10
data$SUBTYPE<-as.factor(data$SUBTYPE)

my_comparisons <- list(c('subtype1','subtype2'),
                       c('subtype1','subtype3'),
                       c('subtype1','subtype4'),
                       c('subtype2','subtype3'),
                       c('subtype2','subtype4'),
                       c('subtype3','subtype4'))
library(ggplot2)
library(ggsignif)
names(data) <- gsub(" ", "_", names(data))

##Plot a boxplot of a single gene
ggplot(data, aes(x=SUBTYPE, y= Neutrophil, color=SUBTYPE)) + 
  geom_boxplot()+
  scale_color_manual(values=c("#C6B3D3","#ED9F9B","#80BA8A","#9CD1C8"))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  geom_signif(comparisons = my_comparisons,step_increase = 0.1,map_signif_level = F,test = t.test)+
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(color="black", fill="transparent"))






