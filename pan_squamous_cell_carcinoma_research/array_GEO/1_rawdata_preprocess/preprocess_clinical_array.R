#Oral
setwd("D:\\鳞癌分析\\GEO 泛鳞癌 临床信息\\Oral")
getwd()
library(tidyr)
library(dplyr)
library(stringr)
###
###GSE3524
phe_GSE3524 = read.csv("GSE3524_series_matrix.csv", header = F)
phe_GSE3524 = phe_GSE3524[c(41:57),]
phe_GSE3524 = phe_GSE3524 %>% t()
colnames(phe_GSE3524) = phe_GSE3524[1,]
phe_GSE3524 = phe_GSE3524[-1,]
rownames(phe_GSE3524) = phe_GSE3524[,2]
phe_GSE3524 = as.data.frame(phe_GSE3524)

colnames(phe_GSE3524)
pheno_GSE3524 = phe_GSE3524[,c(1,2,10)]
View(pheno_GSE3524)
pheno_GSE3524$gender = ifelse(pheno_GSE3524$`!Sample_characteristics_ch1.2` == "Sex:F","F","M")
AA = as.data.frame(str_split_fixed(pheno_GSE3524$`!Sample_characteristics_ch1.1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE3524$`!Sample_characteristics_ch1.3`,":",2))
pheno_GSE3524 = cbind(pheno_GSE3524[,c(1,2,3,8)],AA[,2],BB[,2])
colnames(pheno_GSE3524) = c("Sample_title","ID","source","gender","age","tnm_stage")
pheno_GSE3524$organ = "Oral"
save(pheno_GSE3524,file = "pheno_GSE3524.RData")


###GSE23558
phe_GSE23558 = read.csv("GSE23558_series_matrix.csv", header = F)
phe_GSE23558 = phe_GSE23558[c(32:69),]
phe_GSE23558 = phe_GSE23558 %>% t()
colnames(phe_GSE23558) = phe_GSE23558[1,]
phe_GSE23558 = phe_GSE23558[-1,]
rownames(phe_GSE23558) = phe_GSE23558[,2]
phe_GSE23558 = as.data.frame(phe_GSE23558)

colnames(phe_GSE23558)
pheno_GSE23558 = phe_GSE23558[,c(1,2,11,12,13,14)]
View(pheno_GSE23558)

pheno_GSE23558$gender = ifelse(pheno_GSE23558$`!Sample_characteristics_ch1.1` == "gender: Female","F","M")
BB = as.data.frame(str_split_fixed(pheno_GSE23558$`!Sample_characteristics_ch1.3`,":",2))
pheno_GSE23558 = cbind(pheno_GSE23558,BB[,2])
pheno_GSE23558$age = substr(pheno_GSE23558$`!Sample_characteristics_ch1`,5,7)

colnames(pheno_GSE23558) = c("Sample_title","ID","gender","stage","age")
pheno_GSE23558$organ = "Oral"
save(pheno_GSE23558,file = "pheno_GSE23558.RData")


###GSE31056
phe_GSE31056 = read.csv("GSE31056_series_matrix.csv", header = F)
phe_GSE31056 = phe_GSE31056[c(43:62),]
phe_GSE31056 = phe_GSE31056 %>% t()
colnames(phe_GSE31056) = phe_GSE31056[1,]
phe_GSE31056 = phe_GSE31056[-1,]
rownames(phe_GSE31056) = phe_GSE31056[,2]
phe_GSE31056 = as.data.frame(phe_GSE31056)

colnames(phe_GSE31056)
pheno_GSE31056 = phe_GSE31056[,c(1,2,8,10,12,13,14,15,17,18)]
View(pheno_GSE31056)

#移除远端正常样本，只保留癌旁
A = subset(pheno_GSE31056,pheno_GSE31056$`!Sample_source_name_ch1`=="Normal")
B = subset(pheno_GSE31056,pheno_GSE31056$`!Sample_source_name_ch1`=="Tumor")
pheno_GSE31056 = rbind(A,B)

pheno_GSE31056$sample_type = ifelse(pheno_GSE31056$`!Sample_characteristics_ch1.4` == "site: tumor","T","N")
BB = as.data.frame(str_split_fixed(pheno_GSE31056$`!Sample_characteristics_ch1.3`,":",2))
AA = as.data.frame(str_split_fixed(pheno_GSE31056$`!Sample_characteristics_ch1.6`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE31056$`!Sample_characteristics_ch1.5`,":",2))
DD = as.data.frame(str_split_fixed(pheno_GSE31056$`!Sample_characteristics_ch1.2`,":",2))
EE = as.data.frame(str_split_fixed(pheno_GSE31056$`!Sample_characteristics_ch1.1`,":",2))
pheno_GSE31056 = cbind(pheno_GSE31056[,c(1,2,4,11)],AA[,2],BB[,2],CC[,2],DD[,2],EE[,2])
colnames(pheno_GSE31056) = c("Sample_title","ID","patient","sample_type","recurrence_or_lastfollowup_months","grade","recurrence_states","stage","source")
pheno_GSE31056$organ = "Oral"
save(pheno_GSE31056,file = "pheno_GSE31056.RData")


###GSE37991
phe_GSE37991 = read.csv("GSE37991_series_matrix.csv", header = F)
phe_GSE37991 = phe_GSE37991[c(33:62),]
phe_GSE37991 = phe_GSE37991 %>% t()
colnames(phe_GSE37991) = phe_GSE37991[1,]
phe_GSE37991 = phe_GSE37991[-1,]
rownames(phe_GSE37991) = phe_GSE37991[,2]
phe_GSE37991 = as.data.frame(phe_GSE37991)

pheno_GSE37991 = phe_GSE37991[,c(1,2,8)]
View(pheno_GSE37991)
pheno_GSE37991$sample_type = ifelse(pheno_GSE37991$`!Sample_source_name_ch1`=="oral squamous cell carcinoma","T","N")
pheno_GSE37991 = pheno_GSE37991[,c(1,2,4)]
colnames(pheno_GSE37991) = c("Sample_title","ID","sample_type")
getwd()
pheno_GSE37991$organ = "Oral"
save(pheno_GSE37991,file = "pheno_GSE37991.RData")


###GSE41116
###GSE41116
phe_GSE41116 = read.csv("GSE41116_series_matrix.csv", header = F)
phe_GSE41116 = phe_GSE41116[c(32:44),]
phe_GSE41116 = phe_GSE41116 %>% t()
colnames(phe_GSE41116) = phe_GSE41116[1,]
phe_GSE41116 = phe_GSE41116[-1,]
rownames(phe_GSE41116) = phe_GSE41116[,2]
phe_GSE41116 = as.data.frame(phe_GSE41116)
colnames(phe_GSE41116)
pheno_GSE41116 = phe_GSE41116[,c(1,2,10,11)]
pheno_GSE41116$gender = ifelse(pheno_GSE41116$`!Sample_characteristics_ch1`=="gender: M","M","F")
AA = as.data.frame(str_split_fixed(pheno_GSE41116$`!Sample_characteristics_ch1.1`,":",2))
pheno_GSE41116 = cbind(pheno_GSE41116[,c(1,2,5)],AA[,2])
colnames(pheno_GSE41116) = c("Sample_title","ID","gender","source")
getwd()
pheno_GSE41116$organ = "Oral"
pheno_GSE41116$sample_type = "T"
save(pheno_GSE41116,file = "pheno_GSE41116.RData")


###GSE87593
###GSE87593
phe_GSE87593 = read.csv("GSE87593_series_matrix.csv", header = F)
phe_GSE87593 = phe_GSE87593[c(30:42),]
phe_GSE87593 = phe_GSE87593 %>% t()
colnames(phe_GSE87593) = phe_GSE87593[1,]
phe_GSE87593 = phe_GSE87593[-1,]
rownames(phe_GSE87593) = phe_GSE87593[,2]
phe_GSE87593 = as.data.frame(phe_GSE87593)
colnames(phe_GSE87593)
pheno_GSE87593 = phe_GSE87593[,c(1,2,11,12,13)]
pheno_GSE87593$gender = ifelse(pheno_GSE87593$`!Sample_characteristics_ch1`=="gender: male","M","F")
pheno_GSE87593$recur_state = ifelse(pheno_GSE87593$`!Sample_characteristics_ch1.1` =="recurrence status: Recur","Yes","No")
pheno_GSE87593$sample_type = ifelse(pheno_GSE87593$`!Sample_characteristics_ch1.2`=="tissue: Primary tumor","T","N")
pheno_GSE87593 = pheno_GSE87593[,c(1,2,6,7,8)]
colnames(pheno_GSE87593) = c("Sample_title","ID","gender","recur_state","sample_type")
pheno_GSE87593$organ = "Oral"

save(pheno_GSE87593,file = "pheno_GSE87593.RData")


###GSE138206
###GSE138206
phe_GSE138206 = read.csv("GSE138206_series_matrix.csv", header = F)
phe_GSE138206 = phe_GSE138206[c(27:41),]
phe_GSE138206 = phe_GSE138206 %>% t()
colnames(phe_GSE138206) = phe_GSE138206[1,]
phe_GSE138206 = phe_GSE138206[-1,]
rownames(phe_GSE138206) = phe_GSE138206[,2]
phe_GSE138206 = as.data.frame(phe_GSE138206)
colnames(phe_GSE138206)
pheno_GSE138206 = phe_GSE138206[,c(1,2,8,14,15)]

pheno_GSE138206$sample_type = ifelse(pheno_GSE138206$`!Sample_characteristics_ch1`=="tissue disease state: cancer","T","N")
pheno_GSE138206$gender = ifelse(pheno_GSE138206$`!Sample_characteristics_ch1.1` =="gender: male","M","F")

pheno_GSE138206 = pheno_GSE138206[,c(1,2,3,6,7)]

##去除远端正常组织
pheno_GSE138206 = pheno_GSE138206[c(1:12),]
GSE138206_anno = GSE138206_anno[,rownames(pheno_GSE138206)]

colnames(pheno_GSE138206) = c("Sample_title","ID","source","sample_type","gender")
pheno_GSE138206$organ = "Oral"
save(pheno_GSE138206,file = "pheno_GSE138206.RData")


##Skin

setwd("D:\\鳞癌分析\\GEO 泛鳞癌 临床信息\\skin")
getwd()
library(tidyr)
library(dplyr)
library(stringr)

###GSE2503
phe_GSE2503 = read.csv("GSE2503_series_matrix.csv", header = F)
phe_GSE2503 = phe_GSE2503[c(37:38),]
phe_GSE2503 = phe_GSE2503 %>% t()
colnames(phe_GSE2503) = phe_GSE2503[1,]
phe_GSE2503 = phe_GSE2503[-1,]
rownames(phe_GSE2503) = phe_GSE2503[,2]
phe_GSE2503 = as.data.frame(phe_GSE2503)
colnames(phe_GSE2503)
pheno_GSE2503 = phe_GSE2503[,c(1,2)]

pheno_GSE2503$sample_type = c("actinic keratosis","actinic keratosis","actinic keratosis","actinic keratosis","N","N","N","N","N","N","T","T","T","T","T")
pheno_GSE2503$organ = "skin"
colnames(pheno_GSE2503) = c("Sample_title","ID","sample_type","organ")

save(pheno_GSE2503,file = "pheno_GSE2503.RData")



###GSE45216
phe_GSE45216 = read.csv("GSE45216_series_matrix.csv", header = F)
phe_GSE45216 = phe_GSE45216[c(35:48),]
phe_GSE45216 = phe_GSE45216 %>% t()
colnames(phe_GSE45216) = phe_GSE45216[1,]
phe_GSE45216 = phe_GSE45216[-1,]
rownames(phe_GSE45216) = phe_GSE45216[,2]
phe_GSE45216 = as.data.frame(phe_GSE45216)
colnames(phe_GSE45216)
pheno_GSE45216 = phe_GSE45216[,c(1,2,10,11,12,14)]
pheno_GSE45216$sample_type = ifelse(pheno_GSE45216$`!Sample_characteristics_ch1`=="tissue: Cutaneous SCC","T","actinic keratosis")
pheno_GSE45216$gender = ifelse(pheno_GSE45216$`!Sample_characteristics_ch1.3`=="gender: Male","M","F")
AA = as.data.frame(str_split_fixed(pheno_GSE45216$`!Sample_characteristics_ch1.1`,":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE45216$`!Sample_characteristics_ch1.2`,":",2))
pheno_GSE45216 = cbind(pheno_GSE45216[,c(1,2,7,8)],AA[,2],BB[,2])
pheno_GSE45216$sample_type = c("actinic keratosis","actinic keratosis","actinic keratosis","actinic keratosis","N","N","N","N","N","N","T","T","T","T","T")
colnames(pheno_GSE45216) = c("Sample_title","ID","sample_type","gender","differentiated_state","immune_state")
pheno_GSE45216$organ = "skin"
save(pheno_GSE45216,file = "pheno_GSE45216.RData")



###GSE108010_GPL16131
phe_GSE108010_GPL16131 = read.csv("GSE108010-GPL16131_series_matrix.csv", header = F)
phe_GSE108010_GPL16131 = phe_GSE108010_GPL16131[c(32:46),]
phe_GSE108010_GPL16131 = phe_GSE108010_GPL16131 %>% t()
colnames(phe_GSE108010_GPL16131) = phe_GSE108010_GPL16131[1,]
phe_GSE108010_GPL16131 = phe_GSE108010_GPL16131[-1,]
rownames(phe_GSE108010_GPL16131) = phe_GSE108010_GPL16131[,2]
phe_GSE108010_GPL16131 = as.data.frame(phe_GSE108010_GPL16131)
colnames(phe_GSE108010_GPL16131)
pheno_GSE108010_GPL16131 = phe_GSE108010_GPL16131[,c(1,2,11,12,13,14)]
pheno_GSE108010_GPL16131$age = substr(pheno_GSE108010_GPL16131$`!Sample_characteristics_ch1.2`,6,8)
pheno_GSE108010_GPL16131$gender = ifelse(pheno_GSE108010_GPL16131$`!Sample_characteristics_ch1.1`== "gender: male","M","F")
AA = as.data.frame(str_split_fixed(pheno_GSE108010_GPL16131$`!Sample_characteristics_ch1.3`,":",2))

pheno_GSE108010_GPL16131 = cbind(pheno_GSE108010_GPL16131[,c(1,2,7,8)],AA[,2])
pheno_GSE108010_GPL16131$sample_type = c("N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T",
                                         "N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T",
                                         "blood","blood","blood","blood","blood","blood","blood","blood","blood","blood")


pheno_GSE108010_GPL16131$organ = c("skin","skin","skin","skin","skin","skin","skin","skin","skin","skin",
                                   "skin","skin","skin","skin","skin","skin","skin","skin","skin","skin",
                                   "skin","skin","skin","skin","skin","skin","skin","skin","skin","skin",
                                   "blood","blood","blood","blood","blood","blood","blood","blood","blood","blood")

colnames(pheno_GSE108010_GPL16131) = c("Sample_title","ID","age","gender","source","sample_type","organ")
pheno_GSE108010_GPL16131$immune_state = "immunocompetent"
save(pheno_GSE108010_GPL16131,file = "pheno_GSE108010_GPL16131.RData")

###GSE108010_GPL16686
phe_GSE108010_GPL16686 = read.csv("GSE108010-GPL16686_series_matrix.csv", header = F)
phe_GSE108010_GPL16686 = phe_GSE108010_GPL16686[c(32:46),]
phe_GSE108010_GPL16686 = phe_GSE108010_GPL16686 %>% t()
colnames(phe_GSE108010_GPL16686) = phe_GSE108010_GPL16686[1,]
phe_GSE108010_GPL16686 = phe_GSE108010_GPL16686[-1,]
rownames(phe_GSE108010_GPL16686) = phe_GSE108010_GPL16686[,2]
phe_GSE108010_GPL16686 = as.data.frame(phe_GSE108010_GPL16686)
colnames(phe_GSE108010_GPL16686)
pheno_GSE108010_GPL16686 = phe_GSE108010_GPL16686[,c(1,2,8,12,13,14)]

pheno_GSE108010_GPL16686$age = substr(pheno_GSE108010_GPL16686$`!Sample_characteristics_ch1.1`,6,8)
pheno_GSE108010_GPL16686$gender = ifelse(pheno_GSE108010_GPL16686$`!Sample_characteristics_ch1`== "gender: male","M","F")
pheno_GSE108010_GPL16686 = pheno_GSE108010_GPL16686[c(1,2,6,7,8)]
colnames(pheno_GSE108010_GPL16686) = c("Sample_title","ID","source","age","gender")

pheno_GSE108010_GPL16686$sample_type = c("N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T",
                                         "N","actinic keratosis","T","N","actinic keratosis","T","N","actinic keratosis","T","actinic keratosis","N","T","N","actinic keratosis","T")

pheno_GSE108010_GPL16686$organ = "skin"
save(pheno_GSE108010_GPL16686,file = "pheno_GSE108010_GPL16686.RData")


##Lung

setwd("D:\\鳞癌分析\\GEO 泛鳞癌 临床信息\\lung")
getwd()
library(tidyr)
library(dplyr)
library(stringr)
###GSE4573
phe_GSE4573 = read.csv("GSE4573_series_matrix.csv", header = F)
phe_GSE4573 = phe_GSE4573[c(27:36),]
phe_GSE4573 = phe_GSE4573 %>% t()
colnames(phe_GSE4573) = phe_GSE4573[1,]
phe_GSE4573 = phe_GSE4573[-1,]
rownames(phe_GSE4573) = phe_GSE4573[,2]
phe_GSE4573 = as.data.frame(phe_GSE4573)
colnames(phe_GSE4573)

pheno_GSE4573 = phe_GSE4573[,c(1,2)]
pheno_GSE4573$sample_type = "T"
pheno_GSE4573$organ = "lung"
colnames(pheno_GSE4573) = c("Sample_title","ID","sample_type","organ")

save(pheno_GSE4573,file = "pheno_GSE4573.RData")


###GSE5123
phe_GSE5123 = read.csv("GSE5123_series_matrix.csv", header = F)
phe_GSE5123 = phe_GSE5123[c(36:37),]
phe_GSE5123 = phe_GSE5123 %>% t()
colnames(phe_GSE5123) = phe_GSE5123[1,]
phe_GSE5123 = phe_GSE5123[-1,]
rownames(phe_GSE5123) = phe_GSE5123[,2]
phe_GSE5123 = as.data.frame(phe_GSE5123)
colnames(phe_GSE5123) = c("Sample_title","ID")

pheno_GSE5123 = phe_GSE5123
pheno_GSE5123$organ = "lung"
pheno_GSE5123$sample_type = "T"
pheno_GSE5123$recur_state = c("NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO",
                              "NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO","NO",
                              "Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes",
                              "Yes","Yes","Yes")
save(pheno_GSE5123,file = "pheno_GSE5123.RData")

##GSE10245
phe_GSE10245 = read.csv("GSE10245_series_matrix.csv", header = F)
phe_GSE10245 = phe_GSE10245[c(36:45),]
phe_GSE10245 = phe_GSE10245 %>% t()
colnames(phe_GSE10245) = phe_GSE10245[1,]
phe_GSE10245 = phe_GSE10245[-1,]
rownames(phe_GSE10245) = phe_GSE10245[,2]
phe_GSE10245 = as.data.frame(phe_GSE10245)
pheno_GSE10245 = phe_GSE10245[,c(1,2,10)]

pheno_GSE10245$tumor_type = ifelse(pheno_GSE10245$`!Sample_characteristics_ch1`=="disease state: adenocarcinoma","AC","SCC")
pheno_GSE10245 = pheno_GSE10245[,c(1,2,4)]
colnames(pheno_GSE10245) = c("Sample_title","ID","tumor_type")
pheno_GSE10245$organ = "lung"
pheno_GSE10245$sample_type = "T"

save(pheno_GSE10245,file = "pheno_GSE10245.RData")



##GSE43580
phe_GSE43580 = read.csv("GSE43580_series_matrix.csv", header = F)
phe_GSE43580 = phe_GSE43580[c(35:56),]
phe_GSE43580 = phe_GSE43580 %>% t()
write.csv(phe_GSE43580,file = "phe_GSE43580.csv")#因为数据缺失，单元格错位了，转成excel以后手动删除
phe_GSE43580 = read.csv("phe_GSE43580.csv", header = F)
phe_GSE43580 = phe_GSE43580[,-1]
colnames(phe_GSE43580) = phe_GSE43580[2,]
phe_GSE43580 = phe_GSE43580[-c(1,2),]
rownames(phe_GSE43580) = phe_GSE43580[,2]
phe_GSE43580 = as.data.frame(phe_GSE43580)

colnames(phe_GSE43580)
pheno_GSE43580 = phe_GSE43580[,c(1,2,10,11,15)]

AA = as.data.frame(str_split_fixed(pheno_GSE43580$`!Sample_title`,"_",6))
BB = as.data.frame(str_split_fixed(pheno_GSE43580$`!Sample_characteristics_ch1.2`,":",2))
CC = as.data.frame(str_split_fixed(pheno_GSE43580$`!Sample_characteristics_ch1.1`,":",2))
pheno_GSE43580$gender = ifelse(pheno_GSE43580$`!Sample_characteristics_ch1`=="gender: male","M","F")
pheno_GSE43580 = cbind(pheno_GSE43580[,c(1,2,7)],AA[,c(1,3,5)],BB[,2],CC[,2])
colnames(pheno_GSE43580) = c("Sample_title","ID","gender","organ","tumor_type","stage","smoke_state","age")
pheno_GSE43580$sample_type = "T"

save(pheno_GSE43580,file = "pheno_GSE43580.RData")



##GSE126533
phe_GSE126533 = read.csv("GSE126533_series_matrix.csv", header = F)
phe_GSE126533 = phe_GSE126533[c(29:41),]
phe_GSE126533 = phe_GSE126533 %>% t()
colnames(phe_GSE126533) = phe_GSE126533[1,]
phe_GSE126533 = phe_GSE126533[-1,]
rownames(phe_GSE126533) = phe_GSE126533[,2]
phe_GSE126533 = as.data.frame(phe_GSE126533)
pheno_GSE126533$age = substr(pheno_GSE126533$`!Sample_characteristics_ch1.3`,6,7)
pheno_GSE126533$gender = substr(pheno_GSE126533$`!Sample_characteristics_ch1.2`,9,9)
pheno_GSE126533$sample_type = ifelse(pheno_GSE126533$`!Sample_characteristics_ch1.1`=="tissue: normal lung tissue","N","T")

pheno_GSE126533$tumor_type = "SCC"
pheno_GSE126533=pheno_GSE126533[,c(1,2,3,7,8,9,10)] 

colnames(pheno_GSE126533) = c("Sample_title","ID","patient","age","gender","sample_type","tumor_type")
pheno_GSE126533$organ = "lung"
save(pheno_GSE126533,file = "pheno_GSE126533.RData")


##GSE137291
phe_GSE137291 = read.csv("GSE137291_series_matrix.csv", header = F)
phe_GSE137291 = phe_GSE137291[c(42:43),]
phe_GSE137291 = phe_GSE137291 %>% t()
colnames(phe_GSE137291) = phe_GSE137291[1,]
phe_GSE137291 = phe_GSE137291[-1,]
rownames(phe_GSE137291) = phe_GSE137291[,2]
pheno_GSE137291 = as.data.frame(phe_GSE137291)
colnames(pheno_GSE137291) = c("Sample_title","ID")
pheno_GSE137291$tumor_type = "SCC"
pheno_GSE137291$organ = "lung" 
pheno_GSE137291$sample_type = "T"

save(pheno_GSE137291,file = "pheno_GSE137291.RData")



##Ovary

setwd("D:\\鳞癌分析\\GEO 泛鳞癌 临床信息\\ovary")
getwd()
library(tidyr)
library(dplyr)
library(stringr)

###GSE3578
phe_GSE3578 = read.csv("GSE3578_series_matrix.csv", header = F)
phe_GSE3578 = phe_GSE3578[c(29:30),]
phe_GSE3578 = phe_GSE3578 %>% t()
colnames(phe_GSE3578) = phe_GSE3578[1,]
phe_GSE3578 = phe_GSE3578[-1,]
rownames(phe_GSE3578) = phe_GSE3578[,2]
pheno_GSE3578 = as.data.frame(phe_GSE3578)
AA = as.data.frame(str_split_fixed(pheno_GSE3578$`!Sample_title`,"_",3))
pheno_GSE3578 = cbind(pheno_GSE3578,AA)
pheno_GSE3578 =subset(pheno_GSE3578,pheno_GSE3578$V2=="p" & pheno_GSE3578$V3=="1")
pheno_GSE3578 = pheno_GSE3578[,c(1,2)]
colnames(pheno_GSE3578) = c("Sample_title","ID")
pheno_GSE3578$sample_type = "T"
pheno_GSE3578$organ = "ovary"
save(pheno_GSE3578,file = "pheno_GSE3578.RData")


###GSE9750
phe_GSE9750 = read.csv("GSE9750_series_matrix.csv", header = F)
phe_GSE9750 = phe_GSE9750[c(35:44),]
phe_GSE9750 = phe_GSE9750 %>% t()
colnames(phe_GSE9750) = phe_GSE9750[1,]
phe_GSE9750 = phe_GSE9750[-1,]
rownames(phe_GSE9750) = phe_GSE9750[,2]
pheno_GSE9750 = as.data.frame(phe_GSE9750)
pheno_GSE9750 = pheno_GSE9750[-c(1:9),]   ##剔除不需要细胞系数据
pheno_GSE9750 = pheno_GSE9750[,c(1,2,8,10)]
write.csv(pheno_GSE9750, file = "pheno_GSE9750.csv")
pheno_GSE9750 = read.csv("pheno_GSE9750.csv", header = F)
pheno_GSE9750 = pheno_GSE9750[,c(2,3,6,7)]
pheno_GSE9750 = pheno_GSE9750[-1,]
rownames(pheno_GSE9750) = pheno_GSE9750$V3

colnames(pheno_GSE9750) = c("Sample_title","ID","sample_type","age")

pheno_GSE9750$organ = "ovary"
save(pheno_GSE9750,file = "pheno_GSE9750.RData")


###GSE166466
phe_GSE166466 = read.csv("GSE166466_series_matrix.csv", header = F)
phe_GSE166466 = phe_GSE166466[c(26:33),]
phe_GSE166466 = phe_GSE166466 %>% t()
colnames(phe_GSE166466) = phe_GSE166466[1,]
phe_GSE166466 = phe_GSE166466[-1,]
rownames(phe_GSE166466) = phe_GSE166466[,2]
pheno_GSE166466 = phe_GSE166466[,c(1,2)]
pheno_GSE166466 = pheno_GSE166466[-c(8:13),]
pheno_GSE166466 = as.data.frame(pheno_GSE166466)

colnames(pheno_GSE166466) = c("Sample_title","ID")
pheno_GSE166466$sample_type = c("N","N","N","N","N","N","N","T","T","T","T","T","T","T")
pheno_GSE166466$organ = "ovary"
save(pheno_GSE166466,file = "pheno_GSE166466.RData")


###bladder

setwd("D:\\鳞癌分析\\GEO 泛鳞癌 临床信息\\bladder")
getwd()
library(tidyr)
library(dplyr)
library(stringr)

###GSE186691
phe_GSE186691 = read.csv("GSE186691_series_matrix.csv", header = F)
phe_GSE186691 = phe_GSE186691[c(26:36),]
phe_GSE186691 = phe_GSE186691 %>% t()
colnames(phe_GSE186691) = phe_GSE186691[1,]
phe_GSE186691 = phe_GSE186691[-1,]
rownames(phe_GSE186691) = phe_GSE186691[,2]
pheno_GSE186691 = as.data.frame(phe_GSE186691)
pheno_GSE186691 = pheno_GSE186691[,c(1,2,10)]
pheno_GSE186691 = pheno_GSE186691[c(19:34),] ##只要鳞状细胞癌样本
pheno_GSE186691 = pheno_GSE186691[,c(1,2)]
colnames(pheno_GSE186691) = c("Sample_title","ID")
pheno_GSE186691$sample_type = "T"
pheno_GSE186691$tumor_type = "SCC"
pheno_GSE186691$organ = "bladder"

save(pheno_GSE186691,file = "pheno_GSE186691.RData")


###penis

setwd("D:\\鳞癌分析\\GEO 泛鳞癌 临床信息\\penis")
getwd()
library(tidyr)
library(dplyr)
library(stringr)

###GSE57955
phe_GSE57955 = read.csv("GSE57955_series_matrix.csv", header = F)
phe_GSE57955 = phe_GSE57955[c(27:37),]
phe_GSE57955 = phe_GSE57955 %>% t()
colnames(phe_GSE57955) = phe_GSE57955[1,]
phe_GSE57955 = phe_GSE57955[-1,]
rownames(phe_GSE57955) = phe_GSE57955[,2]
pheno_GSE57955 = as.data.frame(phe_GSE57955)
pheno_GSE57955 = pheno_GSE57955[,c(1,2)]

colnames(pheno_GSE57955) = c("Sample_title","ID")
pheno_GSE57955$sample_type = "T"
pheno_GSE57955$tumor_type = "SCC"
pheno_GSE57955$organ = "penis"

save(pheno_GSE57955,file = "pheno_GSE57955.RData")



###External auditory

setwd("D:\\鳞癌分析\\GEO 泛鳞癌 临床信息\\EACC")
getwd()
library(tidyr)
library(dplyr)
library(stringr)

###GSE98912
phe_GSE98912 = read.csv("GSE98912_series_matrix.csv", header = F)
phe_GSE98912 = phe_GSE98912[c(30:46),]
phe_GSE98912 = phe_GSE98912 %>% t()
colnames(phe_GSE98912) = phe_GSE98912[1,]
phe_GSE98912 = phe_GSE98912[-1,]
rownames(phe_GSE98912) = phe_GSE98912[,2]
pheno_GSE98912 = as.data.frame(phe_GSE98912)
pheno_GSE98912 = pheno_GSE98912[,c(1,2,10,12,13,14,15)]
pheno_GSE98912$gender = ifelse(pheno_GSE98912$`!Sample_characteristics_ch1`=="Sex: female","F","M")
pheno_GSE98912$sample_type = ifelse(pheno_GSE98912$`!Sample_characteristics_ch1.4`=="tissue: External auditory canal squamous cell carcinoma tissue collected from patients during surgery","T","N")
pheno_GSE98912$age = substr(pheno_GSE98912$`!Sample_characteristics_ch1.1`,6,7)
pheno_GSE98912$t_stage = substr(pheno_GSE98912$`!Sample_characteristics_ch1.2`,21,22)
AA = as.data.frame(str_split_fixed(pheno_GSE98912$`!Sample_characteristics_ch1.3`,":",2))
pheno_GSE98912 = pheno_GSE98912[,c(1,2,8,9,10,11)]
pheno_GSE98912 = cbind(pheno_GSE98912,AA[,2])
colnames(pheno_GSE98912) = c("Sample_title","ID","gender","sample_type","age","t_stage","grade")

save(pheno_GSE98912,file = "pheno_GSE98912.RData")




###head and neck

setwd("D:\\鳞癌分析\\GEO 泛鳞癌 临床信息\\head and neck")
getwd()
library(tidyr)
library(dplyr)
library(stringr)

###GSE3292
phe_GSE3292 = read.csv("GSE3292_series_matrix.csv", header = F)
phe_GSE3292 = phe_GSE3292[c(47:48),]
phe_GSE3292 = phe_GSE3292 %>% t()
colnames(phe_GSE3292) = phe_GSE3292[1,]
phe_GSE3292 = phe_GSE3292[-1,]
rownames(phe_GSE3292) = phe_GSE3292[,2]
pheno_GSE3292 = as.data.frame(phe_GSE3292)
colnames(pheno_GSE3292) = c("Sample_title","ID")
pheno_GSE3292$sample_type = "T"
pheno_GSE3292$organ = "head and neck"
save(pheno_GSE3292,file = "pheno_GSE3292.RData")


###GSE6631
phe_GSE6631 = read.csv("GSE6631_series_matrix.csv", header = F)
phe_GSE6631 = phe_GSE6631[c(43:44),]
phe_GSE6631 = phe_GSE6631 %>% t()
colnames(phe_GSE6631) = phe_GSE6631[1,]
phe_GSE6631 = phe_GSE6631[-1,]
rownames(phe_GSE6631) = phe_GSE6631[,2]
pheno_GSE6631 = as.data.frame(phe_GSE6631)
colnames(pheno_GSE6631) = c("Sample_title","ID")
pheno_GSE6631$sample_type = c("N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N",
                              "T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T")
pheno_GSE6631$organ = "head and neck"
save(pheno_GSE6631,file = "pheno_GSE6631.RData")


###GSE23036
phe_GSE23036 = read.csv("GSE23036_series_matrix.csv", header = F)
phe_GSE23036 = phe_GSE23036[c(41:55),]
phe_GSE23036 = phe_GSE23036 %>% t()
colnames(phe_GSE23036) = phe_GSE23036[1,]
phe_GSE23036 = phe_GSE23036[-1,]
rownames(phe_GSE23036) = phe_GSE23036[,2]
pheno_GSE23036 = as.data.frame(phe_GSE23036)

pheno_GSE23036$sample_type = c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T",
                               "T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T",
                               "T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","N","N","N","N"
                               )
AA = as.data.frame(str_split_fixed(pheno_GSE23036$`!Sample_characteristics_ch1`,":",2))
pheno_GSE23036 = cbind(pheno_GSE23036[,c(1,2,16)],AA[,2])
pheno_GSE23036$gender = ifelse(pheno_GSE23036$`AA[, 2]`==" normal mucosa","NA",pheno_GSE23036$`AA[, 2]`)
pheno_GSE23036 = pheno_GSE23036[,c(1,2,3,5)]
colnames(pheno_GSE23036) = c("Sample_title","ID","sample_type","gender")
pheno_GSE23036$organ = "head and neck"
pheno_GSE23036$tumor_type = "SCC"
save(pheno_GSE23036,file = "pheno_GSE23036.RData")


###GSE58911
phe_GSE58911 = read.csv("GSE58911_series_matrix.csv", header = F)
phe_GSE58911 = phe_GSE58911[c(30:45),]
phe_GSE58911 = phe_GSE58911 %>% t()
colnames(phe_GSE58911) = phe_GSE58911[1,]
phe_GSE58911 = phe_GSE58911[-1,]
rownames(phe_GSE58911) = phe_GSE58911[,2]
pheno_GSE58911 = as.data.frame(phe_GSE58911)
pheno_GSE58911$sample_type = ifelse(pheno_GSE58911$`!Sample_source_name_ch1`=="HNSCC","T","N")
pheno_GSE58911$t_stage = substr(pheno_GSE58911[,12],8,9)
AA = as.data.frame(str_split_fixed(pheno_GSE58911[,13],":",2))
BB = as.data.frame(str_split_fixed(pheno_GSE58911[,14],":",2))
pheno_GSE58911 = pheno_GSE58911[,c(1,2,10,17,18)]
pheno_GSE58911 = cbind(pheno_GSE58911,AA[,2],BB[,2])
colnames(pheno_GSE58911) = c("Sample_title","ID","patient","sample_type","t_stage","hpv_state","age")
pheno_GSE58911$tumor_type = "SCC"
pheno_GSE58911$organ = "head and neck"
save(pheno_GSE58911,file = "pheno_GSE58911.RData")



###GSE172120
phe_GSE172120 = read.csv("GSE172120_series_matrix.csv", header = F)
phe_GSE172120 = phe_GSE172120[c(31:41),]
phe_GSE172120 = phe_GSE172120 %>% t()
colnames(phe_GSE172120) = phe_GSE172120[1,]
phe_GSE172120 = phe_GSE172120[-1,]
rownames(phe_GSE172120) = phe_GSE172120[,2]
pheno_GSE172120 = as.data.frame(phe_GSE172120)
pheno_GSE172120 = pheno_GSE172120[c(1:6),] ##剔除细胞系样本

pheno_GSE172120$age = substr(pheno_GSE172120[,10],6,7)
pheno_GSE172120 = pheno_GSE172120[,c(1,2,12)]

colnames(pheno_GSE172120) = c("Sample_title","ID","age")
pheno_GSE172120$sample_type = c("N","N","N","T","T","T")
pheno_GSE172120$tumor_type = "SCC"
pheno_GSE172120$organ = "head and neck"
save(pheno_GSE172120,file = "pheno_GSE172120.RData")




###GSE201777
phe_GSE201777 = read.csv("GSE201777_series_matrix.csv", header = F)
phe_GSE201777 = phe_GSE201777[c(27:38),]
phe_GSE201777 = phe_GSE201777 %>% t()
colnames(phe_GSE201777) = phe_GSE201777[1,]
phe_GSE201777 = phe_GSE201777[-1,]
rownames(phe_GSE201777) = phe_GSE201777[,2]
pheno_GSE201777 = as.data.frame(phe_GSE201777)
pheno_GSE201777 = subset(pheno_GSE201777,pheno_GSE201777[,11] !="tissue: Lymph node") ##剔除淋巴结样本
pheno_GSE201777$sample_type = ifelse(pheno_GSE201777$`!Sample_characteristics_ch1.1`=="tissue: Tumor","T","N")
pheno_GSE201777 = pheno_GSE201777[,c(1,2,12)]
colnames(pheno_GSE201777) = c("Sample_title","ID","sample_type")

pheno_GSE201777$tumor_type = "SCC"
pheno_GSE201777$organ = "head and neck"
save(pheno_GSE201777,file = "pheno_GSE201777.RData")













