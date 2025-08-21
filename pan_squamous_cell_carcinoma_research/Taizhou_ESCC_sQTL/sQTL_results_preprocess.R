##对结果进行FDR校正
##并保留FDR < 0.05 的sQTL事件

library(data.table)

chr1 <- fread("chr1.txt.gz")
colnames(chr1) = c("Intron_id","snp_id","TSS_diatance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se")
chr1$FDR_BH = p.adjust(chr1$pval_nominal, method = "BH")
chr1_sig = subset(chr1,chr1$FDR_BH < 0.000001)



