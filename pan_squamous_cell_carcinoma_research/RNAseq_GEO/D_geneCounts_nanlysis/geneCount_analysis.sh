########get the gene expression data from the bamfiles
##step 1 : 使用samtools 将bam 文件转换成sam文件
#BAM文件：Binary Alignment/Map文件，用于存储比对到参考基因组上的测序数据，是一种二进制格式的文件，通常较小且更快速地处理。
#SAM文件：Sequence Alignment/Map文件，是BAM文件的文本表示形式，易于阅读和处理

##step 2 : 使用samtools 对sam文件进行排序，
##上述两个步骤在hisat align 的过程中已经进行过了，于是直接进行第三个步骤，找到表达量


#!/bin/bash
set -e
# path to annotation in GTF format
anno=/mnt/data/wangdanke/PSCC_RNAseq/reference_genome/Homo_sapiens.GRCh38.110.gtf
# path to the spladder executable
featureCounts=/mnt/data/wangdanke/Programes/subread-2.0.6-Linux-x86_64/bin/featureCounts
## text file containing the absolute paths to the alignment files in BAM format
# alignments.txt

cat ~/PSCC_RNAseq/OSCC/GSE176221/GSE176221.txt | while read id
do
mkdir ./Counts/$id/
$featureCounts  -a $anno -p --countReadPairs -t exon -g gene_id -o ./Counts/$id/"$id".txt ./del_unmap/$id/*.bam 
done







