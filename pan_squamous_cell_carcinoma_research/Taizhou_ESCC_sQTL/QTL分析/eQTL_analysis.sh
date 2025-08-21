#!/bin/sh 
#$ -S /bin/sh
#$ -N N_sQTL      #任务名
#$ -V                
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./		            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q normal.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq normal.q@ibnode14    # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 4-4        #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR
set -e

CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
      chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
      chr21 chr22)

FastQTL=/mnt/data/wangdanke/fastqtl/bin/fastQTL

# 参数设置

for chr in "${CHRS[@]}"; do

  $FastQTL \
    --vcf normal.snps.withid.vcf.gz \
    --bed ./RNA_expr/ESCC_Normal_TPMNormalised_Counts_for_eQTL_"$chr".bed.gz \
    --cov covariates_N_expr.txt \
    --out results/"$chr".txt.gz \
    --window 1000000 \
    --region $chr
    
done







