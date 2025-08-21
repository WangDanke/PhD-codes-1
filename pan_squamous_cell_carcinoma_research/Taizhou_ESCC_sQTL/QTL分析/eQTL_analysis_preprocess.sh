#1. 首先确定 以转录本起始位点（TSS）上下游 1Mb作为筛选eQTL位点的窗口范围；准备注释基因的TSS位置注释文件
#.sh
#基于Homo_sapiens.GRCh38.110.gtf找到各个转录本的编码id、TSS位点、正负链以及染色体信息。(这个文件也是最开始识别基因表达量的时候使用的注释文件)
awk -F'\t' '
    $3 == "gene" {
        # 提取 gene_id（去掉引号）
        split($9, attr, ";");
        gene_id = "";
        for (i in attr) {
            if (attr[i] ~ /gene_id/) {
                gsub(/"/, "", attr[i]);      # 去掉引号
                gsub(/gene_id /, "", attr[i]); # 去掉字段名
                gene_id = attr[i];
                break;
            }
        }
        # 确定 TSS（正链=start，负链=end）
        tss = ($7 == "+") ? $4 : $5;
        # 输出：chr, gene_id, TSS, strand
        print $1, gene_id, tss, $7;
    }
' Homo_sapiens.GRCh38.110.gtf > gene_tss_strand.txt

#2. 整理基因表达量数据框;分别按照肿瘤组织和癌旁组织进行文件整理
#.R
conda activtae R4.3
R 

library(data.table)
library(dplyr)

ID_list = read.table("./normal_list",header = F,row.names = NULL)[,1]
#ID_list = read.table("./tumor_list",header = F,row.names = NULL)[,1]

Normal_ESCC_express <- data.frame(matrix(nrow = 62754, ncol = 44))
colnames(Normal_ESCC_express) = ID_list

# 循环遍历ID列表
for (file_id in ID_list) {
  # 构建文件路径
  file_path <- paste("/mnt/data/wangdanke/PSCC_RNAseq/ESCC_Taizhou/RNA_seq/Counts/", file_id,"/",file_id,".txt", sep = "")
  # 使用read.table()函数读取文本文件
  data<- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  Normal_ESCC_express[,file_id] = data[,7]
  rownames(Normal_ESCC_express) = data$Geneid
  Normal_ESCC_express$gene_length = data$Length
}

save(Normal_ESCC_express,file = "Normal_ESCC_express.RData")


for (file_id in ID_list) {
  # 构建文件路径
  file_path <- paste("/mnt/data/wangdanke/PSCC_RNAseq/ESCC_Taizhou/RNA_seq/Counts/", file_id,"/",file_id,".txt", sep = "")
  # 使用read.table()函数读取文本文件
  data<- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  Tumor_ESCC_express[,file_id] = data[,7]
  rownames(Tumor_ESCC_express) = data$Geneid
  Normal_ESCC_express$gene_length = data$Length
}

save(Tumor_ESCC_express,file = "Tumor_ESCC_express.RData")


#标准化基因表达量 TPM + 对数转换（log2);分别对癌旁组织和癌组织进行
load("./Normal_ESCC_express.RData")
#load("./Tumor_ESCC_express.RData")

kb = Normal_ESCC_express$gene_length/1000
counts <- Normal_ESCC_express %>%
  select(-gene_length)
counts = as.matrix(counts)

rpk = counts/kb
tpm = t(t(rpk)/colSums(rpk)*1000000)

#在超过一半的样本中的tpm值大于0.1的才保留下来
keep <- rowSums(tpm > 0.1) >= 0.5 * ncol(tpm)
tpm_filter = tpm[keep,]
#对数转换
tpm_filter_log = log2(tpm_filter+1)

Normal_ESCC_express_normalized = tpm_filter_log
save(Normal_ESCC_express_normalized,file="Normal_ESCC_express_normalized.RData")


#3. 计算前十个主成分作为协变量文件；分为癌组织和癌旁组织分别进行主成分分析和计算
load("Normal_ESCC_express_normalized")
#load("Tumor_ESCC_express_normalized")

#  转置矩阵（确保行是样本，列是基因）
exp_matrix <- t(Normal_ESCC_express_normalized)

# 计算PCA
pca_result <- prcomp(exp_matrix, center = TRUE, scale. = TRUE)

# 提取前10个主成分
top10_pcs <- pca_result$x[, 1:10]

# 保存为协变量文件（可选）
write.csv(top10_pcs, file = "top10_PCs_covariates_Normal_ESCC.csv", row.names = TRUE)


#4. 合并TSS注释文件和基因表达信息；分为癌组织和癌旁组织分别进行
load("Normal_ESCC_express_normalized.RData")
load("Tumor_ESCC_express_normalized.RData")

anno = read.table("gene_tss_strand.txt")
colnames(anno) = c("chromosome","feature_id","TSS","chain")
Normal_ESCC_express_normalized = as.data.frame(Normal_ESCC_express_normalized)
Normal_ESCC_express_normalized$feature_id = rownames(Normal_ESCC_express_normalized)

bed_Normal <- left_join(Normal_ESCC_express_normalized, anno%>%select(chromosome, TSS, feature_id)) %>%
  mutate(end=TSS)%>% 
  mutate(start=TSS-1)%>%
  select(-TSS)%>%
  select(Chr = chromosome, start, end, ID = feature_id, everything()) %>%
  arrange(Chr, start) %>%
  filter(!is.na(Chr)) %>%
  rename("#Chr" = Chr)

fwrite(bed_Normal, sep = "\t", file = "ESCC_Normal_TPMNormalised_Counts_for_eQTL.bed")

##由于前面的基因组VCF文件中染色体的格式是chr* ,所以要将我们得到的bed文件中的染色体改成chr+数字 而不是纯数字格式
#sh
awk 'BEGIN{OFS="\t"}
{
    if ($1 != "#Chr" && NR > 1) {
        $1 = "chr" $1
    }
    print
}' ESCC_Tumor_TPMNormalised_Counts_for_eQTL.bed > ESCC_Tumor_TPMNormalised_Counts_for_eQTL.bed


#5. 将基因表达信息按照染色体区分开；Output per-chromosome phenotype files for cis and trans analysis；排除性染色体
$conda activate R4.3
$R

bed_Normal = read.delim
for(i in as.character(1:22)) {
  bedChr <- bed_Normal %>%
    rename(Chr = "#Chr") %>%
    filter(Chr == paste0("chr", i)) %>%  # 修改这里
    rename("#Chr" = Chr)
  fwrite(bedChr, sep = "\t", file = paste0("ESCC_Normal_TPMNormalised_Counts_for_eQTL_chr",i,".bed"))
  rm(bedChr)
}

#将每一个染色体bed文件构建索引
#.sh脚本
#!/bin/bash
for i in {1..22}
do
  bgzip ESCC_Normal_TPMNormalised_Counts_for_eQTL_chr${i}.bed && tabix -p bed ESCC_Normal_TPMNormalised_Counts_for_eQTL_chr${i}.bed.gz
done



