setwd("D:\\桌面")
getwd()

pheno_BLCA = read.table("clinical_BLCA.tsv",header = T,sep = "\t",quote = "")
pheno_CESC = read.table("clinical_CESC.tsv",header = T,sep = "\t",quote = "")
pheno_ESCC = read.table("clinical_ESCC.tsv",header = T,sep = "\t",quote = "")
pheno_HNSC = read.table("clinical_HNSC.tsv",header = T,sep = "\t",quote = "")
pheno_LUSC = read.table("clinical_LUSC.tsv",header = T,sep = "\t",quote = "")

PSI_BLCA = read.table("PSI_download_BLCA.txt",header = T,sep = "",quote = "")
PSI_ESCA = read.table("PSI_download_ESCA.txt",header = T,sep = "",quote = "")
PSI_LUSC = read.table("PSI_download_LUSC.txt",header = T,sep = "",quote = "")
PSI_HNSC = read.table("PSI_download_HNSC.txt",header = T,sep = "",quote = "")
PSI_CESC = read.table("PSI_download_CESC.txt",header = T,sep = "",quote = "")


###下载的临床信息有重复，先去重
install.packages("magrittr")
library(magrittr)
library(dplyr)

phe_BLCA <- pheno_BLCA %>% group_by(case_submitter_id) %>% filter(!duplicated(case_submitter_id))
phe_ESCC <- pheno_ESCC %>% group_by(case_submitter_id) %>% filter(!duplicated(case_submitter_id))
phe_LUSC <- pheno_LUSC %>% group_by(case_submitter_id) %>% filter(!duplicated(case_submitter_id))
phe_HNSC <- pheno_HNSC %>% group_by(case_submitter_id) %>% filter(!duplicated(case_submitter_id))
phe_CESC <- pheno_CESC %>% group_by(case_submitter_id) %>% filter(!duplicated(case_submitter_id))


###
phe_BLCA <- phe_HNSC %>% as.data.frame()
rownames(phe_BLCA) <- phe_BLCA$case_submitter_id

phe_ESCC <- phe_ESCC %>% as.data.frame()
rownames(phe_ESCC) <- phe_ESCC$case_submitter_id

phe_LUSC <- phe_LUSC %>% as.data.frame()
rownames(phe_LUSC) <- phe_LUSC$case_submitter_id

phe_HNSC <- phe_HNSC %>% as.data.frame()
rownames(phe_HNSC) <- phe_HNSC$case_submitter_id

phe_CESC <- phe_CESC %>% as.data.frame()
rownames(phe_CESC) <- phe_CESC$case_submitter_id

save(phe_BLCA,phe_CESC,phe_ESCC,phe_HNSC,phe_LUSC,file = "pheno.RData")
save(PSI_BLCA,PSI_CESC,PSI_ESCA,PSI_HNSC,PSI_LUSC,file = "PSI_TCGA.RData")

###合并不同部位来源的鳞癌样本的临床信息，并挑选需要的变量
pheno = rbind(phe_BLCA,phe_CESC,phe_ESCC,phe_HNSC,phe_LUSC)

phe = pheno[, c("case_submitter_id","project_id","age_at_index","days_to_birth","days_to_death",
                  "gender","race","vital_status","age_at_diagnosis" ,"ajcc_pathologic_n","ajcc_clinical_m",
                  "ajcc_pathologic_t","classification_of_tumor","days_to_best_overall_response","days_to_diagnosis",
                "primary_diagnosis","primary_disease","tumor_regression_grade","treatment_type"
                  )]

###找到匹配的临床信息和PSI数据
BLCA = intersect(colnames(phe),rownames(PSI_BLCA))
ESCC = intersect(rownames(phe),colnames(PSI_ESCA))
LUSC = intersect(rownames(phe),colnames(PSI_LUSC))
HNSC = intersect(rownames(phe),colnames(PSI_HNSC))
CESC = intersect(rownames(phe),colnames(PSI_CESC))

Psi_BLCA = cbind(PSI_BLCA[,1:10],PSI_BLCA[,BLCA]) 
Psi_LUSC = cbind(PSI_LUSC[,1:10],PSI_LUSC[,LUSC])
Psi_HNSC = cbind(PSI_HNSC[,1:10],PSI_HNSC[,HNSC])
Psi_ESCA = cbind(PSI_ESCA[,1:10],PSI_ESCA[,ESCC])
Psi_CESC = cbind(PSI_CESC[,1:10],PSI_CESC[,CESC])

phe_SCC = phe[c(BLCA, LUSC, HNSC, ESCC, CESC),]
#保存数据
save(phe_SCC,pheno, file = "pheno_2.RData")
save(Psi_BLCA,Psi_CESC,Psi_ESCA,Psi_HNSC,Psi_LUSC, file = "PSI-SCC-TCGA.RData")



####将每个基因的每种剪切占比以及位置等，作为每一个单独的变量

rownames(Psi_BLCA) = paste(Psi_BLCA$symbol,Psi_BLCA$splice_type, Psi_BLCA$exons,Psi_BLCA$from_exon, Psi_BLCA$to_exon, sep = "_")
rownames(Psi_CESC) = paste(Psi_CESC$symbol,Psi_CESC$splice_type, 
                           Psi_CESC$exons,Psi_CESC$from_exon, Psi_CESC$to_exon, sep = "_")
rownames(Psi_HNSC) = paste(Psi_HNSC$symbol,Psi_HNSC$splice_type, 
                           Psi_HNSC$exons,Psi_HNSC$from_exon, Psi_HNSC$to_exon, sep = "_")
rownames(Psi_ESCA) = paste(Psi_ESCA$symbol,Psi_ESCA$splice_type, 
                           Psi_ESCA$exons,Psi_ESCA$from_exon, Psi_ESCA$to_exon, sep = "_")

rownames(Psi_LUSC) = paste(Psi_LUSC$symbol,Psi_LUSC$splice_type, 
                           Psi_LUSC$exons,Psi_LUSC$from_exon, Psi_LUSC$to_exon, sep = "_"


###确定共享的唯一变量
c = intersect(rownames(Psi_BLCA),rownames(Psi_CESC),rownames(Psi_ESCA),rownames(Psi_HNSC),rownames(Psi_LUSC))
data = cbind(Psi_BLCA[c,-c(1:10)], Psi_CESC[c,-c(1:10)],Psi_ESCA[c,-c(1:10)],Psi_HNSC[c,-c(1:10)],Psi_LUSC[c,-c(1:10)])



###降维分析

tsne_out_sel = Rtsne(
  t(data_sel),
  dims = 2,
  pca = T,
  max_iter = 1000,
  theta = 0.4,
  perplexity = 20,
  verbose = F
)


