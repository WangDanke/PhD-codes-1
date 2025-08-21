##########prepare the data for permutation test
###
EAC = list(datExpr=EAC_EXPR, datMeta=pheno_EAC)
LUAD = list(datExpr=LUAD_EXPR, datMeta=pheno_LUAD)
LUSC = list(datExpr=LUSC_EXPR, datMeta=pheno_LUSC)
CSCC = list(datExpr=CSCC_EXPR, datMeta=pheno_CSCC)
ESCC = list(datExpr=ESCC_EXPR, datMeta=pheno_ESCC)
HNSC = list(datExpr=HNSC_EXPR, datMeta=pheno_HNSC)
CESC = list(datExpr=CESC_EXPR, datMeta=pheno_CESC)

multiExpr = list(LUSC=LUSC,CSCC=CSCC,ESCC=ESCC,HNSC=HNSC,CESC=CESC,EAC=EAC,LUAD=LUAD,)
save("D:\\鳞癌分析\\泛鳞癌array数据\\8_DEG_Meta\\null_distribution\\Microarray_compiledForPermutationTesting.RData")


##2b_calculateNullDistribution
set.seed(20240704)
options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R"); biocLite("lme4")
library(lme4)
rootdir = "D:\\鳞癌分析\\泛鳞癌array数据\\8_DEG_Meta\\null_distribution"
setwd(rootdir)

load("D:/鳞癌分析/泛鳞癌array数据/8_DEG_Meta/null_distribution/Microarray_compiledForPermutationTesting.RData") ##Compiled expression & metadata from 2a

genes = rownames(multiExpr[[1]]$datExpr)
allmeta = matrix(NA,nrow=length(genes), length(multiExpr))
colnames(allmeta) = c("LUSC","CSCC","ESCC","HNSC","CESC","EAC","LUAD")
allmeta=  as.data.frame(allmeta)

lmer_apply=function(x, datMeta) {
  if(length(unique(datMeta$ID)) < nrow(datMeta)) {
    return(summary(lmer(x ~ Group + batch + (1 | ID),data=datMeta))$coefficients[2,1])
  } else if(length(unique(datMeta$ID)) > 1) {
    return(summary(lm(x ~ Group + batch,data=datMeta))$coefficients[2,1])
  } else {
    return(summary(lm(x ~ Group + batch,data=datMeta))$coefficients[2,1])
  }
}

n_iterations = 1000
for (iteration in 1:n_iterations) {
  print(paste("Iteration", iteration))
  for(i in 1:length(multiExpr)) {
    tt = matrix(NA, nrow=length(genes), ncol=3)
    subj = unique(as.character(multiExpr[[i]]$datMeta$ID))
    subj_group = data.frame(Subject = subj, Group = multiExpr[[i]]$datMeta$Group[match(subj, multiExpr[[i]]$datMeta$ID)])
    subj_group$Group = subj_group$Group[order(runif(nrow(subj_group)))] ##Randomly shuffle group assignment for each subject
    multiExpr[[i]]$datMeta$Group = subj_group$Group[match(multiExpr[[i]]$datMeta$ID,subj_group$Subject)]
    
    allmeta[,i] = apply(multiExpr[[i]]$datExpr,1,lmer_apply, multiExpr[[i]]$datMeta)
  }
  cor_vec = vector(mode="numeric")
  comparisons = t(combn(seq(1,ncol(allmeta)),2))
  
  for(i in 1:nrow(comparisons)) {
    r = cor(allmeta[,comparisons[i,1]], allmeta[,comparisons[i,2]], method = "spearman", use="pairwise.complete.obs")
    cor_vec = c(cor_vec,r)
  }
  # 
  save(cor_vec, file=paste("D:/鳞癌分析/泛鳞癌array数据/8_DEG_Meta/null_distribution/iteration_", iteration, ".RData", sep=""), row.names=F, col.names=F)
}





##########绘图
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\DEG_Meta\\null_distribution")

# 创建一个空白的数据框

file_list <- list.files(pattern = ".RData")

# 初始化一个空数据框
merged_data <- data.frame()

# 循环读取并合并数据
for (file_name in file_list) {
  loaded_data <- load(file_name)
  merged_data = rbind(merged_data,cor_vec)
}

colnames(merged_data) = c("LUSC_CSCC","LUSC_ESCC","LUSC_HNSC","LUSC_CESC","LUSC_EAC","LUSC_LUAD","CSCC_ESCC","CSCC_HNSC","CSCC_CESC","CSCC_EAC","CSCC_LUAD",
"ESCC_HNSC","ESCC_CESC","ESCC_EAC","ESCC_LUAD","HNSC_CESC","HNSC_EAC","HNSC_LUAD","CESC_EAC","CESC_LUAD","EAC_LUAD")
save(merged_data,file = "merge_cor.RData")