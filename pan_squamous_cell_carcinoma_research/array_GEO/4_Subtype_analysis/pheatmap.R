
annotation_col = pheno_SCC_TCGA[,c("gender",'project_id',"ajcc_clinical_t","Subtype_3SCC_10")]
write.csv(annotation_col,file = "annotation.csv")
##保存以后，在Excel里将分期不明的样本填上NA，并按照亚型排序。
annotation_col = read.csv("annotation.csv",header = T,row.names = 1)
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(annotation_col)]
annotation_col<- annotation_col %>%
  mutate(Subtype_3SCC_10 = case_when(
    Subtype_3SCC_10 == 1 ~ "subtype1",
    Subtype_3SCC_10 == 2 ~ "subtype2",
    Subtype_3SCC_10 == 3 ~ "subtype3",
    Subtype_3SCC_10 == 4 ~ "subtype4"
  ))

library(pheatmap)
pheatmap(SCC_TCGA_hub ,cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =T,
         annotation_col = annotation_col)


library(org.Hs.eg.db)
library(clusterProfiler)
GO_hub_gene_top百分之10 <-enrichGO(gene=hub_gene$gene,OrgDb = "org.Hs.eg.db",ont="ALL",
                     pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff=0.2,keyType="SYMBOL")

write.csv(GO_hub_gene_top百分之10@result,file = "GO_hub_genes_top10.csv")


setwd("D:/鳞癌分析/泛鳞癌array数据/12_Subtype_analysis")
library(readxl)
pheno = read_excel("多组学分型-信息-泛鳞癌-TCGA.xlsx") 
pheno = pheno[,c("Patient Barcodes","HPV_status","Total Mutation Rate (per Mb)","CNV clusters (6)",
                 "MethylMix clusters (5)","miRNA clusters (5)","MDSC clusters (4)")]


pheno_merge<-merge(x=pheno,y=annotation_col,by='Patient Barcodes',all.x=F,all.y=T)

Immune_subtype = read.csv("Immune_subtype.csv",header = T,row.names = 1)
Immune_subtype$"Patient Barcodes" = rownames(Immune_subtype)

pheno_merge<-merge(x=pheno_merge,y=Immune_subtype,by='Patient Barcodes',all.x=T,all.y=F)




chisq.test(subtype_cell_reports_filter$my_subtype,subtype_cell_reports_filter$HPV_status)
subtype_cell_reports_filter$Immune_subtype = A$Immune.subtype
chisq.test(subtype_cell_reports_filter$my_subtype,subtype_cell_reports_filter$Immune_subtype)


setwd("D:/鳞癌分析/泛鳞癌array数据/12_Subtype_analysis/")
##保存以后，在Excel里将分期不明的样本填上NA，并按照亚型排序。
TCGA_SCC_EXPR = TCGA_SCC_EXPR[,-c(1:3)]
SCC_TCGA_hub = log(TCGA_SCC_EXPR+1)
annotation_col = read.csv("annotation.csv",header = T,row.names = 1)
annotation_row =  read.csv("hub_gene_100.csv",header = T,row.names = 1)
pheno_merge = pheno_merge[rownames(annotation_col),]
annotation_col = pheno_merge

SCC_TCGA_hub = SCC_TCGA_hub[,rownames(annotation_col)]
SCC_TCGA_hub = SCC_TCGA_hub[rownames(annotation_row),]

SCC_TCGA_hub = sweep(SCC_TCGA_hub,1, apply(SCC_TCGA_hub,1,median,na.rm=T))

annotation_col = annotation_col[,c(11,10,8,12,2,9,13,4,5,6,7)]
annotation_col$Subtype_3SCC_10 = as.factor(annotation_col$Subtype_3SCC_10)
annotation_col$ajcc_clinical_t = as.factor(annotation_col$ajcc_clinical_t)
annotation_col$gender = as.factor(annotation_col$gender)
annotation_col$HPV_status = as.factor(annotation_col$HPV_status)
annotation_col$project_id = as.factor(annotation_col$project_id)
annotation_col$Immune.subtype = as.factor(annotation_col$Immune.subtype)
annotation_col$`CNV clusters (6)` = as.factor(annotation_col$`CNV clusters (6)`)
annotation_col$`MethylMix clusters (5)` = as.factor(annotation_col$`MethylMix clusters (5)`)
annotation_col$`miRNA clusters (5)` = as.factor(annotation_col$`miRNA clusters (5)`)
annotation_col$`MDSC clusters (4)` = as.factor(annotation_col$`MDSC clusters (4)`)


annotation_colors = list(
  module = c("CD1"= "turquoise","CD2"="blue","CD3"="brown","CD4"="yellow","CD5"="green","CD6"="red",
             "CD7" = "black","CD8"="pink","CD9"="magenta"),
  gender = c("female"="#FFC1C1","male"="#87CEEB"),
  project_id = c("TCGA-CESC"="#778899","TCGA-ESCA" = "#CD5555","TCGA-HNSC"="#CDAA7D"),
  ajcc_clinical_t = c("NA" = "#E8E8E8","T1"="#CFCFCF","T2" = "#B5B5B5","T2" = "#9C9C9C","T3" = "#828282",
                      "T4" = "#696969","T4a" = "#4F4F4F","T4b"="#363636"),
  Subtype_3SCC_10 = c("subtype1"="#C6B3D3","subtype2"="#ED9F9B","subtype3"="#80BA8A","subtype4"="#9CD1C8"),
  HPV_status = c("positive"="#363636","negative" ="#CFCFCF","NA" = "white" )
  
)

library(pheatmap)
SCC_TCGA_hub[SCC_TCGA_hub< -2.5] <- -2.5
SCC_TCGA_hub[SCC_TCGA_hub>2.5] <- 2.5
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(annotation_col)]
pheatmap(SCC_TCGA_hub ,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("blue","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,show_rownames =T,fontsize_row = 5, 
         annotation_col = annotation_col,annotation_row = annotation_row,annotation_colors = annotation_colors)

chisq.test(annotation_col$Subtype_3SCC_10,annotation_col$Immune.subtype)
subtype_cell_reports_filter$Immune_subtype = A$Immune.subtype
chisq.test(subtype_cell_reports_filter$my_subtype,subtype_cell_reports_filter$Immune_subtype)





















