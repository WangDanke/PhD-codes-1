###两个脚本，进行23条染色体jointcalling的并行处理

#第一个脚本：run_joint_calling.sh 主调度脚本
#!/bin/bash
# 染色体列表（可根据实际情况修改）
CHRS=(chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
      chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
      chr21 chr22 chrX chrY)

for chr in "${CHRS[@]}"; do
  qsub -N JointCall_${chr} -v CHR=${chr} run_joint_calling_by_chr.sh
done

#第二个脚本：run_joint_calling_by_chr.sh 单个染色体的joint calling脚本

#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o logs/
#$ -e logs/
#$ -q normal.q
#$ -pe thread 4
#$ -l h_vmem=8G

source ~/.bashrc
export JAVA_HOME=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env
export PATH=$JAVA_HOME/bin:$PATH

gatk=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env/bin/gatk
REF=/mnt/data/Biodatabase/ref_human/hg38/variant_calling/Homo_sapiens_assembly38.fasta

# 确保染色体名称由环境变量 CHR 提供
if [ -z "$CHR" ]; then
  echo "Error: CHR variable is not set"
  exit 1
fi

echo "▶▶▶ Start GenomicsDBImport for $CHR"

# 1. GenomicsDBImport
$gatk GenomicsDBImport \
  --sample-name-map sample_map.txt \
  --genomicsdb-workspace-path gendb_${CHR} \
  -L $CHR \
  -R $REF

# 2. GenotypeGVCFs
$gatk GenotypeGVCFs \
  -R $REF \
  -V gendb://gendb_chr1 \
  -O joint_called_chr1.vcf

bgzip -f joint_called_chr1.vcf
bcftools index joint_called_chr1.vcf.gz

echo "✅ Completed joint calling for $CHR"

####在运行时，
#先创建 log文件夹
mkdir -p logs/

##提交运行脚本
#bash
./run_joint_calling.sh


##################################################################################################################################################
##再使用两个脚本，分别对每个染色体进行高质量的变异过滤
##主调度脚本：提交每个染色体的过滤任务 
#submit_variant_filtering.sh

#!/bin/bash
CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
      chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
      chr21 chr22 chrX chrY)

for chr in "${CHRS[@]}"; do
  qsub -N Filter_${chr} -v CHR=${chr} run_variant_filtering_by_chr.sh
done

##单个染色体的过滤脚本：run_variant_filtering_by_chr.sh

#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o logs/
#$ -e logs/
#$ -q normal.q
#$ -pe thread 2
#$ -l h_vmem=8G

source ~/.bashrc
export JAVA_HOME=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env
export PATH=$JAVA_HOME/bin:$PATH

gatk=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env/bin/gatk
REF=/mnt/data/Biodatabase/ref_human/hg38/variant_calling/Homo_sapiens_assembly38.fasta

# 变异过滤标准
FILTER_EXPR="QD < 2.0 || FS > 60.0 || MQ < 40.0"
FILTER_NAME="BasicFilter"

# 确保变量存在
if [ -z "$CHR" ]; then
  echo "❌ Error: CHR variable is not set"
  exit 1
fi

INPUT_VCF=./joint_called_${CHR}.vcf.gz
OUTPUT_VCF=./joint_called_${CHR}.filtered.vcf.gz

if [ ! -f "$INPUT_VCF" ]; then
  echo "❌ Error: Input file $INPUT_VCF not found"
  exit 1
fi

echo "▶▶▶ Start VariantFiltration for $CHR"

$gatk VariantFiltration \
  -R $REF \
  -V $INPUT_VCF \
  -O $OUTPUT_VCF \
  --filter-expression "$FILTER_EXPR" \
  --filter-name "$FILTER_NAME"

echo "✅ Completed VariantFiltration for $CHR"

##运行方法
#bash
bash submit_variant_filtering.sh


########################################################################################################################################
#合并过滤后的多个染色体 VCF文件
bcftools concat -Oz -o joint_called_allchr.filtered.vcf.gz \
  joint_called_chr1.filtered.vcf.gz \
  joint_called_chr2.filtered.vcf.gz \
  joint_called_chr3.filtered.vcf.gz \
  joint_called_chr4.filtered.vcf.gz \
  joint_called_chr5.filtered.vcf.gz \
  joint_called_chr6.filtered.vcf.gz \
  joint_called_chr7.filtered.vcf.gz \
  joint_called_chr8.filtered.vcf.gz \
  joint_called_chr9.filtered.vcf.gz \
  joint_called_chr10.filtered.vcf.gz \
  joint_called_chr11.filtered.vcf.gz \
  joint_called_chr12.filtered.vcf.gz \
  joint_called_chr13.filtered.vcf.gz \
  joint_called_chr14.filtered.vcf.gz \
  joint_called_chr15.filtered.vcf.gz \
  joint_called_chr16.filtered.vcf.gz \
  joint_called_chr17.filtered.vcf.gz \
  joint_called_chr18.filtered.vcf.gz \
  joint_called_chr19.filtered.vcf.gz \
  joint_called_chr20.filtered.vcf.gz \
  joint_called_chr21.filtered.vcf.gz \
  joint_called_chr22.filtered.vcf.gz \
  joint_called_chrX.filtered.vcf.gz \
  joint_called_chrY.filtered.vcf.gz

###提取SNP信息用于后续的sQTL分析
bcftools view -v snps -m2 -M2 -i 'MAF>0.05' joint_called_allchr.filtered.vcf.gz -Oz -o normal.filtered.snps.vcf.gz
bcftools view -v snps -m2 -M2 -i 'MAF>0.05' joint_called_allchr.filtered.vcf.gz -Oz -o tumor.filtered.snps.vcf.gz