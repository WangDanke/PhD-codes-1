library(survival)
library(survminer)
library(preprocessCore)

setwd("D:\\鳞癌分析\\泛鳞癌array数据\\12_Subtype_analysis\\")
load("D:/鳞癌分析/泛鳞癌array数据/12_Subtype_analysis/TCGA_4_SCCs_8816x1169.RData")

hub_gene = read.csv("9_module_hub_gene_10%.csv",header = T)
TCGA_SCC_EXPR = TCGA_SCC_EXPR[,-c(1,2,3)]
SCC_TCGA_hub = TCGA_SCC_EXPR[hub_gene$gene,]
SCC_TCGA_hub = log2(SCC_TCGA_hub+1)
pheno_SCC_TCGA = subset(pheno_SCC_TCGA,!pheno_SCC_TCGA$project_id=="TCGA-LUSC")##去除LUSC样本
pheno_SCC_TCGA$time = pheno_SCC_TCGA$time/365
SCC_TCGA_hub = SCC_TCGA_hub[,rownames(pheno_SCC_TCGA)]

data_1 = SCC_TCGA_hub


###前一半取作为高的组后一半作为低的组
IGF2BP3<- data_1[c('IGF2BP3'),]
pheno_SCC_TCGA <- cbind(pheno_SCC_TCGA,as.data.frame(t(IGF2BP3)))
pheno_SCC_TCGA$Group_IGF2BP3<- ifelse(pheno_SCC_TCGA$IGF2BP3>median(pheno_SCC_TCGA$IGF2BP3),"High","Low")

##Figure 6
sfit <- survfit(Surv(time, event)~Group_IGF2BP3, data=pheno_SCC_TCGA)
print(sfit)
ggsurv = ggsurvplot(sfit, data=pheno_SCC_TCGA, 
                    conf.int = F,                            # 添加置信区间
                    pval = TRUE,                                # 添加p值
                    fun = "pct",                                # 将y轴转变为百分比的格式     
                    legend = c(0.8, 0.85),                      # 改变legend的位置
                    legend.title = "IGF2BP3 ",                    # 改变legend的题目
                    legend.labs = c("High", "Low"),          # 改变legend的标签
                    surv.median.line = "hv") 
ggsurv$plot <- ggsurv$plot + 
  scale_color_manual(values = c("High" = "red", "Low" = "black"))# 设置分组颜色
ggsurv









##分癌种看看
#pheno_HNSC = subset(pheno_SCC_TCGA,pheno_SCC_TCGA$project_id=="TCGA-HNSC")
#data_1 = SCC_TCGA_hub[,rownames(pheno_HNSC)]

##按照前四分之一做高组后四分之一做低组做生存分析
##
COL1A1 <- data_1[c('COL1A1'),]
pheno_HNSC <- cbind(pheno_HNSC,as.data.frame(t(COL1A1)))
q1 <- quantile(pheno_HNSC$COL1A1, 0.25)
q3 <- quantile(pheno_HNSC$COL1A1, 0.75)
pheno_HNSC$Group_COL1A1 <- ifelse(
  pheno_HNSC$COL1A1 > q3, "High",
  ifelse(pheno_HNSC$COL1A1< q1, "Low", "Medium")
)
pheno_COL1A1= subset(pheno_HNSC,pheno_HNSC$Group_COL1A1!="Medium")
##Figure 6
sfit <- survfit(Surv(time, event)~Group_COL1A1, data=pheno_COL1A1)
print(sfit)
ggsurv = ggsurvplot(sfit, data=pheno_COL1A1, 
                    conf.int = F,                            # 添加置信区间
                    pval = TRUE,                                # 添加p值
                    fun = "pct",                                # 将y轴转变为百分比的格式     
                    legend = c(0.8, 0.85),                      # 改变legend的位置
                    legend.title = "COL1A1 ",                    # 改变legend的题目
                    legend.labs = c("High", "Low"),          # 改变legend的标签
                    surv.median.line = "hv") 
ggsurv$plot <- ggsurv$plot + 
  scale_color_manual(values = c("High" = "red", "Low" = "black"))# 设置分组颜色
ggsurv



###使用cox回归批量筛选与生存相关的基因
survival_cancer = cbind(pheno_SCC_TCGA[,c("time","event")],t(SCC_TCGA_hub))
survival_cancer<-survival_cancer[which(survival_cancer$time!=0),]


set.seed(20240829)
data_1 = t(survival_cancer[,-c(1,2)])
index <-  sort(ncol(data_1), ncol(data_1)*0.6)


y <- pheno_train[,c('time','event')]
y$time <- as.double(y$time)
y$event <- as.double(y$event)
y <- as.matrix(survival::Surv(y$time,y$event))

fit <- glmnet(as.data.frame(t(train)),y, family = 'cox', alpha = 1)
plot(fit,xvar = "lambda",label=T)
lasso_fit <- cv.glmnet(survival_cancer[,-c(1,2)],y,family ="cox", type.measure = "deviance")
plot(lasso_fit,xvar = "lambda",label=T)
###Lasso回归筛选具有代表性的变量
coefficient <- coef(lasso_fit,s = lasso_fit$lambda.min)
Active.Index <- which(as.numeric(coefficient)!= 0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
sig_gene_multi_cox 












