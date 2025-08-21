###由于一开始的RNA_seq的样本名是以TZ开头的，而WES的样本名是以TX开头的，因此将样本名统一一下

##将基因型文件的列样本名统一成TZ***_align


###过滤出SNP信息
bcftools view -v ./Tumour/joint_called_allchr.filtered.vcf.gz -Oz -o tumor.snps.vcf.gz
bcftools view -v ./Normal/joint_called_allchr.filtered.vcf.gz -Oz -o normal.snps.vcf.gz

###snpde vcf文件中注释上id列
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT'   -o normal.snps.withid.vcf.gz -O z normal.snps.vcf.gz
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT'   -o tumor.snps.withid.vcf.gz -O z tumor.snps.vcf.gz

##进行QTL分析
#安装fastqtl, ##也是安装了好几天才安装好的
FastQTL = /mnt/data/wangdanke/fastqtl/bin/fastQTL 

###nominal格式 即没有进行置换检验
#!/bin/sh 
#$ -S /bin/sh
#$ -N T_sQTL      #任务名
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
      chr21 chr22 chrX chrY)

FastQTL=/mnt/data/wangdanke/fastqtl/bin/fastQTL

# 参数设置

for chr in "${CHRS[@]}"; do

  $FastQTL \
    --vcf tumor.snps.withid.vcf.gz \
    --bed ./RNA_spli/TaizhouESCC_T_perind.counts.gz.qqnorm_"$chr".gz \
    --cov covariates_T.txt \
    --grp TaizhouESCC_T.groups.txt \
    --out results/chrX.txt.gz \
    --window 100000 \
    --region $chr
    
done

##进行置换检验 permute 1000 10000,意思是是进行1000次置换检验，遇到特别显著的就进行10000 次检验

#!/bin/sh 
#$ -S /bin/sh
#$ -N N_sQTL      #任务名
#$ -V                
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./		            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q normal.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq normal.q@ibnode15    # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 10-10        #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR
set -e

CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
      chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
      chr21 chr22 chrX chrY)

FastQTL=/mnt/data/wangdanke/fastqtl/bin/fastQTL

# 参数设置

for chr in "${CHRS[@]}"; do

  $FastQTL \
    --vcf normal.snps.withid.vcf.gz \
    --bed ./RNA_spli/TaizhouESCC_N_perind.counts.gz.qqnorm_"$chr".gz \
    --cov covariates_N.txt \
    --grp TaizhouESCC_N.groups.txt \
    --out results_permute/"$chr".txt.gz \
    --window 100000 \
    --region $chr \
    --permute 1000 
    
done



##单个分别进行 试试看的哦
/mnt/data/wangdanke/fastqtl/bin/fastQTL \
    --vcf tumor.snps.withid.vcf.gz \
    --bed ./RNA_spli/TaizhouESCC_T_perind.counts.gz.qqnorm_chr7.gz \
    --cov covariates_T.txt \
    --grp TaizhouESCC_T.groups.txt \
    --out results_permute/chr7.txt.gz \
    --region chr7 \
    --permute 1000


/mnt/data/wangdanke/fastqtl/bin/fastQTL \
    --vcf tumor.snps.withid.vcf.gz \
    --bed ./RNA_spli/TaizhouESCC_T_perind.counts.gz.qqnorm_chrX.gz \
    --cov covariates_T.txt \
    --grp TaizhouESCC_T.groups.txt \
    --out results/chrX.txt.gz \
    --region chrX





####z注释permutation之后出来的结果
#!/bin/sh 
#$ -S /bin/sh
#$ -N T_anno      #任务名
#$ -V                
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./		            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q normal.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq normal.q@ibnode15    # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 10-10        #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR
set -e

annotate_qtl_results=/mnt/data/wangdanke/fastqtl/python/annotate_outputs.py
snp_info=/mnt/data/wangdanke/PSCC_RNAseq/ESCC_Taizhou/SQTL_analysis/Tumor/snp_info.tsv
permutation_results=/mnt/data/wangdanke/PSCC_RNAseq/ESCC_Taizhou/SQTL_analysis/Tumor/results_permuta/merged_T_permute.txt.gz
annotation_gtf=/mnt/data/wangdanke/PSCC_RNAseq/sQTL_anno_reference_genom/gencode.v38.annotation.gtf
nominal_results=/mnt/data/wangdanke/PSCC_RNAseq/ESCC_Taizhou/SQTL_analysis/Tumor/results/all_nominal_results.txt.gz
python=~/anaconda3/envs/s_qtl_env/bin/python


mkdir ./annotated_output/
$python $annotate_qtl_results \
  $permutation_results \
  0.05 \
  --snp_lookup $snp_info \
  --nominal_results $nominal_results \
  -o ./annotated_output/
done



###先用一个染色体的的结果注释一下看看
python /mnt/data/wangdanke/fastqtl/python/annotate_outputs.py \
  ./results_permute/chr1.txt.gz \
  0.05 \
  /mnt/data/wangdanke/PSCC_RNAseq/sQTL_anno_reference_genom/gencode.v38.annotation.gtf \
  --nominal_results ./results/chr1.txt.gz \
  -o annotated_output/
































