##字符串操作
library(tidyr)
library(dplyr)
library(stringr)

##提取特定长度的字符串
pheno_GSE108010_GPL16686$age = substr(pheno_GSE108010_GPL16686$`!Sample_characteristics_ch1.1`,6,7)##提取！sample_characteristic 列字符串的第六个和第七个作为age列。


#剔除特定的行
new_test_df_subset <- subset(test_df, job != 'Unknown')#将job列中的Unknow行都剔除掉


#更新BiocManager,有些包只有更高级时才可以安装
#更新
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18") #version 是多少可以直接在Bioconductor网站上找