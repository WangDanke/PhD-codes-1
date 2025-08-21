#sQTL分析的整体步骤
#RNA-seq (fastq) → 剪接事件定量 (LeafCutter)
#WES (fastq) → 基因型矩阵 (GATK/SNP-calling)
#→ sQTL 关联分析 (FastQTL/QTLtools)

##工具包安装
conda create -n s_qtl_env python=3.9 -y
conda activate s_qtl_env
conda install -c bioconda samtools fastqc multiqc star hisat2 \
    subread bcftools gatk4 bedtools tabix htslib -y
conda install -c bioconda leafcutter fastqtl

##还要安装部分R包
conda install -c conda-forge r-base r-data.table r-optparse r-ggplot2 r-ggrepel -y

###RNA-seq (fastq) → 剪接事件定量 (LeafCutter)
#step1：质控和比对
# 质控
fastqc *.fastq.gz
multiqc .

# 比对（推荐使用 STAR 或 HISAT2）
STAR --runThreadN 12 --genomeDir /path/to/genome_index \
     --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix sample_ \
     --outSAMtype BAM SortedByCoordinate


##step2:提取剪接事件（使用LeafCutter)
# 提取junction
bam2junc.sh sample.bam > sample.junc

# 聚合所有样本的 junction 文件
ls *.junc > juncfiles.txt

# 聚类剪接位点
python3 /path/to/leafcutter/clustering/leafcutter_cluster.py \
    -j juncfiles.txt -m 50 -o leafcutter_output


##这个part是

#WES (fastq) → 基因型矩阵 (GATK/SNP-calling)
# 比对（推荐bwa）
bwa mem ref.fa sample_R1.fastq.gz sample_R2.fastq.gz | samtools sort -o sample.bam
samtools index sample.bam

# 标记重复、BQSR、SNP Calling（用GATK）
gatk MarkDuplicates -I sample.bam -O sample.marked.bam -M sample.metrics.txt

#SNP calling
gatk HaplotypeCaller -R ref.fa -I sample.marked.bam -O sample.g.vcf.gz -ERC GVCF

# 合并gVCF并进行 joint calling
gatk CombineGVCFs -R ref.fa --variant sample1.g.vcf.gz --variant sample2.g.vcf.gz ... -O cohort.g.vcf.gz
gatk GenotypeGVCFs -R ref.fa -V cohort.g.vcf.gz -O cohort.vcf.gz

# 过滤变异
gatk VariantFiltration -V cohort.vcf.gz -O cohort.filtered.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "BasicFilter"






