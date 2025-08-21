
#!/bin/sh 
#$ -S /bin/sh
#$ -N hisat2      #任务名
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
hisat2=/mnt/data/wangdanke/anaconda3/bin/hisat2
samtools=/mnt/data/wangdanke/Programes/samtools-1.19/samtools


echo "$(date) -- Alignment start"
cat ~/PSCC_RNAseq/HNSC/GSE184616/GSE184616_30 | while read id
do
mkdir ./Align_1st/$id/
$hisat2 -p 10 -x ~/PSCC_RNAseq/GRCh38_reference/GCA_000001405.15_GRCh38_full_analysis_set ./fq_qc_rm_rRNA/$id/*_1.fq.gz -2 ./fq_qc_rm_rRNA/$id/*_2.fq.gz | $samtools view -bS - | $samtools sort -o ./Align_1st/$id/"$id"_align.bam -
echo "$(date) -- Alignment done"
done


##运行脚本
bash hisat2_alignment.sh >hisat2.log 2>&1

#单个样本
/mnt/data/wangdanke/anaconda3/bin/hisat2 -p 10 -x ~/PSCC_RNAseq/grch38_snp_tran/genome_snp_tran -1 ./fq_qc_rm_rRNA/SRR16612817/*_1.fq.gz -2 ./fq_qc_rm_rRNA/SRR16612817/*_2.fq.gz | /mnt/data/wangdanke/Programes/samtools-1.19/samtools view -bS - | /mnt/data/wangdanke/Programes/samtools-1.19/samtools sort -o ./Align_1st/SRR16612817/SRR16612817_align.bam -



