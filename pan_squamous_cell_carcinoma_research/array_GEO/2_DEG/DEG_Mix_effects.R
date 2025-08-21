#####batch as a covariate in the final mixed-effect model
##DEG analysis
setwd("D:\\鳞癌分析\\泛鳞癌array数据\\DEG_Meta")

library(nlme)

####cervical
pheno_Cer_SCC$sample_type = factor(pheno_Cer_SCC$sample_type, levels=c("N", "T"))
pheno_Cer_SCC$batch = as.factor(pheno_Cer_SCC$batch)

bd_meta = matrix(NA, nrow=nrow(Cer_SCC_EXPR), ncol=3)
for(i in 1:nrow(Cer_SCC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(Cer_SCC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type, data = pheno_Cer_SCC, random = ~1|batch))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(Cer_SCC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
Cervical_bd_Meta = bd_meta
save(Cervical_bd_Meta,file = "Cervical_bd_Meta.RData")


####lung
pheno_LUSC$sample_type = factor(pheno_LUSC$sample_type, levels=c("N", "T"))
pheno_LUSC$batch = as.factor(pheno_LUSC$batch)

bd_meta = matrix(NA, nrow=nrow(LUSC_EXPR), ncol=3)
for(i in 1:nrow(LUSC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(LUSC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type, data = pheno_LUSC, random = ~1|batch))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(LUSC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
lung_bd_Meta = bd_meta
save(lung_bd_Meta,file = "lung_bd_Meta.RData")


####head_and _neck
pheno_HNSC$sample_type = factor(pheno_HNSC$sample_type, levels=c("N", "T"))
pheno_HNSC$batch = as.factor(pheno_HNSC$batch)

bd_meta = matrix(NA, nrow=nrow(HNSC_EXPR), ncol=3)
for(i in 1:nrow(HNSC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(HNSC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type, data = pheno_HNSC, random = ~1|batch))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(HNSC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
head_and_neck_bd_Meta = bd_meta
save(head_and_neck_bd_Meta,file = "head_and_neck_bd_Meta.RData")


####head_and _neck
pheno_ESCC$sample_type = factor(pheno_ESCC$sample_type, levels=c("N", "T"))
pheno_ESCC$batch = as.factor(pheno_ESCC$batch)

bd_meta = matrix(NA, nrow=nrow(ESCC_EXPR), ncol=3)
for(i in 1:nrow(ESCC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(ESCC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type, data = pheno_ESCC, random = ~1|batch))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(ESCC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
esophagus_bd_Meta = bd_meta
save(esophagus_bd_Meta,file = "esophagus_bd_Meta.RData")

####skin
pheno_cSCC$sample_type = factor(pheno_cSCC$sample_type, levels=c("N", "T"))
pheno_cSCC$batch = as.factor(pheno_cSCC$batch)

bd_meta = matrix(NA, nrow=nrow(cSCC_EXPR), ncol=3)
for(i in 1:nrow(cSCC_EXPR)) {
  if(i%%100==0) print(i)
  expr = as.numeric(cSCC_EXPR[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ sample_type, data = pheno_cSCC, random = ~1|batch))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = rownames(cSCC_EXPR)
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
skin_bd_Meta = bd_meta
save(skin_bd_Meta,file = "skin_bd_Meta.RData")

####合并多种鳞癌的结果
AA = intersect(rownames(Cervical_bd_Meta),rownames(esophagus_bd_Meta))
BB = intersect(rownames(head_and_neck_bd_Meta),rownames(lung_bd_Meta))
CC = intersect(AA,rownames(skin_bd_Meta))
DD = intersect(CC,BB)

DEG_beta = cbind(Cervical_bd_Meta[DD,]$beta,skin_bd_Meta[DD,]$beta,esophagus_bd_Meta[DD,]$beta,
                 head_and_neck_bd_Meta[DD,]$beta,lung_bd_Meta[DD,]$beta)

DEG_beta = as.data.frame(DEG_beta)
rownames(DEG_beta)=DD
colnames(DEG_beta) = c("Cervical","Cutaneous","Esophagus","Head_and_neck","Lung")

#############相关性分析
tdc = cor(DEG_beta,method = "pearson")
testRes = cor.mtest(DEG_beta, method="pearson",conf.level = 0.95)


library(PerformanceAnalytics)
chart.Correlation(DEG_beta)

library(corrplot)
#绘制一个下三角的热图
corrplot(tdc, type="upper", order="hclust", tl.col="black", tl.srt=45)






