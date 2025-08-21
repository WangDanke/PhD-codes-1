
#创建示例数据框
data <- as.data.frame(t(SCC_TCGA_hub))
# 创建一个空向量，用于存储每一列的标准差
std_dev <- vector()
# 循环计算每一列的标准差
for (col in colnames(data)) {
  std_dev <- c(std_dev, sd(data[[col]]))
}
# 输出每一列的标准差
print(std_dev)
std_dev = as.data.frame(std_dev)
rownames(std_dev)=rownames(SCC_TCGA_hub)
std_dev$gene = rownames(std_dev)
sorted <- as.data.frame(std_dev[order(std_dev$std_dev, decreasing = TRUE), ])
std_top800 = sorted[c(1:100),]$gene
##
rownames(hub_gene)=hub_gene$gene
hub_gene_top_SD_100 = hub_gene[std_top800,]
write.csv(hub_gene_top_SD_100,file = "hub_gene_100.csv")


SCC_TCGA_hub = sweep(SCC_TCGA_hub,1, apply(SCC_TCGA_hub,1,median,na.rm=T))

##保存以后，在Excel里将分期不明的样本填上NA，并按照亚型排序。
annotation_col = read.csv("annotation.csv",header = T,row.names = 1)
annotation_row =  read.csv("hub_gene_100.csv",header = T,row.names = 1)
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(annotation_col)]
SCC_TCGA_hub = SCC_TCGA_hub[rownames(annotation_row),]

 library(dplyr)
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




library(pheatmap)
SCC_TCGA_hub[SCC_TCGA_hub< -5] <- -5
SCC_TCGA_hub[SCC_TCGA_hub>5] <- 5
pheatmap(SCC_TCGA_hub ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("blue","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =T,fontsize_row = 5, 
         annotation_col = annotation_col,annotation_row = annotation_row,annotation_colors = annotation_colors)







###########验证
pheno_SCC_T = subset(pheno_SCC,pheno_SCC$sample_type=="T")
pheno_SCC_4 = subset(pheno_SCC_T,!pheno_SCC_T$organ=="lung")
pheno_SCC_3 = subset(pheno_SCC_4,!pheno_SCC_4$organ=="skin")
SCC_EXPR_T_3 = SCC_EXPR[,rownames(pheno_SCC_3)]


hub_gene = read.csv("9_module_hub_gene_10%.csv",header = T,row.names = 1)
SCC_GEO_hub = SCC_EXPR_T_3[rownames(hub_gene),]

##归一化操作
library(ConsensusClusterPlus)
SCC_GEO_hub = sweep(SCC_GEO_hub,1, apply(SCC_GEO_hub,1,median,na.rm=T))
SCC_GEO_hub = as.matrix(SCC_GEO_hub)
title="maxK_10_GEO_3SCC_hub_10"
Cluster_GEO = ConsensusClusterPlus(SCC_GEO_hub,maxK=10,reps=50,pItem=0.8,pFeature=1,
                                   title = title,clusterAlg="hc",distance="pearson",seed=1262118,plot="png")
###将分型信息与临床信息匹配上
AA = as.data.frame(Cluster_GEO[[4]][["consensusClass"]])
pheno_SCC_3$Subtype = AA$`Cluster_GEO[[4]][["consensusClass"]]`


library(dplyr)
pheno_SCC_3<-pheno_SCC_3 %>%
  mutate(Subtype = case_when(
    Subtype == 1 ~ "subtype2",
    Subtype == 2 ~ "subtype1",
    Subtype == 3 ~ "subtype3",
    Subtype == 4 ~ "subtype4"
  ))

pheno_SCC_3<-pheno_SCC_3 %>%
  mutate(organ = case_when(
    organ == "esophagus" ~ "GEO-ESCC",
    organ == "head and neck" ~ "GEO-HNSC",
    organ == "Oral" ~ "GEO-HNSC",
    organ == "ovary" ~ "GEO-CESC"
  ))


pheno_SCC_3 <- pheno_SCC_3[order(pheno_SCC_3$Subtype), ]
annotation_col = pheno_SCC_3[,c(4,5,7)]
annotation_row =  read.csv("hub_gene_100.csv",header = T,row.names = 1)
##
SCC_GEO_hub = SCC_GEO_hub[,rownames(annotation_col)]
SCC_GEO_hub = SCC_GEO_hub[rownames(annotation_row),]


annotation_colors = list(
  module = c("CD1"= "turquoise","CD2"="blue","CD3"="brown","CD4"="yellow","CD5"="green","CD6"="red",
             "CD7" = "black","CD8"="pink","CD9"="magenta"),
  organ = c("GEO-CESC"="#778899","GEO-ESCC" = "#CD5555","GEO-HNSC"="#CDAA7D"),
  Subtype = c("subtype1"="#C6B3D3","subtype2"="#ED9F9B","subtype3"="#80BA8A","subtype4"="#9CD1C8")
  
)


library(pheatmap)
range(SCC_GEO_hub)
SCC_GEO_hub[SCC_GEO_hub< -2.5] <- -2.5
SCC_GEO_hub[SCC_GEO_hub>2.5] <- 2.5
pheatmap(SCC_GEO_hub ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("blue","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =T,fontsize_row = 5,
         annotation_col = annotation_col,annotation_row = annotation_row,annotation_colors = annotation_colors)



################将hub 验证hub基因的预后
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\12_Subtype_analysis\\")
load("D:/鳞癌分析/泛鳞癌array数据/12_Subtype_analysis/TCGA_SCC_EXPR_8816x1168.RData")
hub_gene = read.csv("9_module_hub_gene_10%.csv",header = T)

TCGA_SCC_EXPR = TCGA_SCC_EXPR[,-c(1,2,3)]
SCC_TCGA_hub = TCGA_SCC_EXPR[hub_gene$gene,]
SCC_TCGA_hub = log2(SCC_TCGA_hub+1)

pheno_SCC_TCGA = subset(pheno_SCC_TCGA,!pheno_SCC_TCGA$project_id=="TCGA-LUSC")##去除LUSC样本
pheno_SCC_TCGA$time = pheno_SCC_TCGA$time/365
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(pheno_SCC_TCGA)]

data_1 = SCC_TCGA_hub
library(survival)  # R自带，无需安装   
library(survminer)
##
MMP12 <- data_1[c('MMP12'),]
pheno_SCC_TCGA <- cbind(pheno_SCC_TCGA,as.data.frame(t(MMP12)))
q1 <- quantile(pheno_SCC_TCGA$MMP12, 0.25)
q3 <- quantile(pheno_SCC_TCGA$MMP12, 0.75)
pheno_SCC_TCGA$Group_MMP12 <- ifelse(
  pheno_SCC_TCGA$MMP12 > q3, "High",
  ifelse(pheno_SCC_TCGA$MMP12< q1, "Low", "Medium")
)
pheno_MMP12= subset(pheno_SCC_TCGA,pheno_SCC_TCGA$Group_MMP12!="Medium")
##Figure 6
sfit <- survfit(Surv(time, event)~Group_MMP12, data=pheno_MMP12)
print(sfit)
ggsurv = ggsurvplot(sfit, data=pheno_MMP12, 
           conf.int = F,                            # 添加置信区间
           pval = TRUE,                                # 添加p值
           fun = "pct",                                # 将y轴转变为百分比的格式     
           legend = c(0.8, 0.85),                      # 改变legend的位置
           legend.title = "MMP12 ",                    # 改变legend的题目
           legend.labs = c("High", "Low"),          # 改变legend的标签
           surv.median.line = "hv") 
ggsurv$plot <- ggsurv$plot + 
  scale_color_manual(values = c("High" = "red", "Low" = "black"))# 设置分组颜色
ggsurv


########绘制这些有预后作用的基因在4个亚型中的表达量
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\12_Subtype_analysis\\")
load("D:/鳞癌分析/泛鳞癌array数据/12_Subtype_analysis/TCGA_SCC_EXPR_8816x1168.RData")
hub_gene = read.csv("9_module_hub_gene_10%.csv",header = T)

TCGA_SCC_EXPR = TCGA_SCC_EXPR[,-c(1,2,3)]
SCC_TCGA_hub = TCGA_SCC_EXPR[hub_gene$gene,]
SCC_TCGA_hub = log2(SCC_TCGA_hub+1)

pheno_SCC_TCGA = subset(pheno_SCC_TCGA,!pheno_SCC_TCGA$project_id=="TCGA-LUSC")##去除LUSC样本
pheno_SCC_TCGA$time = pheno_SCC_TCGA$time/365
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(pheno_SCC_TCGA)]




data = SCC_TCGA_hub[,rownames(annotation_col)]

data <- as.data.frame(data)
data <- as.data.frame(t(data))
data$SUBTYPE<-annotation_col$Subtype_3SCC_10
data$SUBTYPE<-as.factor(data$SUBTYPE)

compaired <- list(c("subtype1", "subtype2"),c("subtype1","subtype3"),c("subtype1","subtype4"),
                  c("subtype2", "subtype3"),c("subtype2","subtype4"),c("subtype3","subtype4"))

##Plot a boxplot of a single gene
ggplot(data, aes(x=SUBTYPE, y=KRT5, color=SUBTYPE)) + 
  geom_boxplot()+
  scale_color_manual(values=c("subtype1"="#C6B3D3","subtype2"="#ED9F9B","subtype3"="#80BA8A","subtype4"="#9CD1C8"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=0.5)+
  geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)+
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(color="black", fill="transparent"))





























