###首先看，直接下载来的数据的质量如何
#!/bin/bash

##首先将压缩文件解压,fastqc 不解压缩也可以
gunzip ./fq/*/*


###试了很久，上面的循环代码也还是跑不通 于是自己用命令行搞定了
#首先先到目标文件夹里面（cd ./GSE....)，并且构建一个存储fastqc之后结果的文件夹（fastqc_result) #> fastqc.log 2>&1 这部分是将QC过程中的日志记录下来，后面可以更好的去看过程中出了什么错

fastqc -o ./fastqc_report/ ./fq/*/* > fastqc.log 2>&1

###使用multiqc将生成的多个报告合并到一起，方便后续查阅，并作出进一步的质控步骤
##到达目标文件夹，上述便是cd ./fastqc_report/
multiqc .


########
##使用trim_galore进行切除接头和剔除低质量碱基
##循环脚本
set -e
cat ~/PSCC_RNAseq/OSCC/GSE176221/GSE176221.txt | while read id
do
mkdir ./fq_qc/$id/
echo " trim_galore cut adapters started at $(date)"
trim_galore -q 25 --phred33 --stringency 3 --length 20 -e 0.1 --fastqc --illumina \
            --paired ./fq/$id/*_1.fastq.gz ./fq/$id/*_2.fastq.gz \
            --gzip -o ./fq_qc/$id/
echo "trim_galore cut adapters finished at $(date)"
done

##调度系统做trim_galore
脚本如下：
#!/bin/sh 
#$ -S /bin/sh
#$ -N qc_trim_galore      #任务名
#$ -V                
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./		            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q fat.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq fat.q@ibfat1    # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 64-64         #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR:$path

set -e
cat ~/PSCC_RNAseq/OSCC/GSE176221/GSE176221.txt | while read id
do
mkdir ./fq_qc/$id/
echo " trim_galore cut adapters started at $(date)"
trim_galore -q 25 --phred33 --stringency 3 --length 20 -e 0.1 --fastqc --illumina \
            --paired ./fq/$id/*_1.fastq.gz ./fq/$id/*_2.fastq.gz \
            --gzip -o ./fq_qc/$id/
echo "trim_galore cut adapters finished at $(date)"
done


#看过滤之后的报告，发现序列前端有一部分的序列basecontent 对不上，因此切除前14个bp碱基，-o 输出文件名，-u 目标序列数
#单个样本
cutadapt -u 14 -o trim_1.fastq SRR14740626_1_val_1.fq.gz

#循环脚本
##去除reads 的前14bp,因为Per Base Sequence Content不对
set -e
cat ~/PSCC_RNAseq/OSCC/GSE186775/GSE186775.txt | while read id
do
mkdir ./fq_qc_cutadapt/$id/
cutadapt -u 14 -o ./fq_qc_cutadapt/$id/"$id"_1_cut14.fq.gz ./fq_qc/$id/*_1.fq.gz
cutadapt -u 14 -o ./fq_qc_cutadapt/$id/"$id"_2_cut14.fq.gz ./fq_qc/$id/*_2.fq.gz
done

###GC 含量per sequence GC content好像不太对，于是去除其中的rRNA
###去除RNA-seq原始数据中的rRNA,很奇怪的一点是，最近conda好像坏了，只能用绝对路径来调用软件
##循环脚本 使用bowtie2 和hisat2 都可以
#原理是：下载需要的人类rRNA参考基因组，将我们的reads 比对上去，将其中未比对上去的reads保留下来作为我们后续RNA-seq分析的数据
#使用Bowtie2，不知道为什么没办法直接调用bowtie2,需要通过绝对路径将bowtie调用上
#下载、构建rRNA索引
/mnt/data/wangdanke/Programes/bowtie2-2.5.2/bowtie2-build sequence.fasta rRNA
##比对脚本
set -e

bowtie2=/mnt/data/wangdanke/Programes/bowtie2-2.5.2/bowtie2

cat ~/PSCC_RNAseq/OSCC/GSE186775/GSE186775.txt | while read id
do
mkdir ./fq_qc_rm_rRNA/$id/
$bowtie2 --very-sensitive-local --no-unal -I 1 -X 1000 -p 6 -x rRNA -1 ./fq_qc/$id/*1_cut14.fq.gz -2 ./fq_qc/$id/*2_cut14.fq.gz --un-conc-gz ./fq_qc_rm_rRNA/$id/"$id"_rRNAremoved.fq.gz 2>./fq_qc_rm_rRNA/$id/"$id"_Map2rRNAStat.xls
done  

##bowtie2参数解析：
##--

#再次使用fastqc 看最终的数据质量


###对于hisat2比对之后的bam文件，对其去除其中未比对到人类基因组上的reads,保证后续定量更完美

#!/bin/bash
set -e 

samtools=/mnt/data/wangdanke/Programes/samtools-1.19/samtools
align_bam=/mnt/data/wangdanke/PSCC_RNAseq/OSCC/GSE186775/Align_1st

echo "samtools to delete unmapped sequences start at $(date)"
cat ~/PSCC_RNAseq/OSCC/GSE186775/GSE186775.txt | while read id
do
mkdir ./fq_qc_del_unmap/$id/
$samtools view -F 4 -b ${align_bam}/$id/*_align.bam > ./fq_qc_del_unmap/$id/${id}_del.bam
echo "delete unmapped sequences finished at $(date)"
done











##############去除PCR过程中产生的重复，使用picard，好像RNA-seq数据不需要进行这一步的质控
##an安装
#wgey

java -jar picard-tools-1.119/MarkDuplicates.jar REMOVE_DUPLICATES=true I=H1hesc_Input_Rep1_chr12_aln.bam O=H1hesc_Input_Rep1_chr12_aln.dedup.bam M=H1hesc.duplicates.log
java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=output.bam \
      M=marked_dup_metrics.txt \
      REMOVE_DUPLICATES=true
##去除线粒体RNA
samtools index H1hesc_Input_Rep1_chr12_aln.dedup.bam
samtools idxstats H1hesc_Input_Rep1_chr12_aln.dedup.bam > H1hesc_Input_Rep1_chr12_aln.dedup.mitochondrial.stats
samtools view -h H1hesc_Input_Rep1_chr12_aln.dedup.bam | grep -v 'chrM' | samtools view -bS -o H1hesc.final.bam

