###RNA-seq数据的比对

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
hisat2=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env/bin/hisat2
samtools=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env/bin/samtools

echo "$(date) -- Alignment start"
cat ~/PSCC_RNAseq/ESCC_Taizhou/sample_list | while read id
do
    echo "Processing $id"
    mkdir -p ./Align_hisat2/$id/

    # 获取read1和read2的fq文件名（确保唯一）
    read1=$(ls ./fq_cut10/$id/*R1_cut.fq.gz)
    read2=$(ls ./fq_cut10/$id/*R2_cut.fq.gz)

    # 对每个样本进行比对并排序
    $hisat2 -p 10 -x ~/PSCC_RNAseq/GRCh38_reference/GCA_000001405.15_GRCh38_full_analysis_set \
        -1 $read1 -2 $read2 | \
    $samtools view -bS - | \
    $samtools sort -o ./Align_hisat2/$id/"$id"_align.bam -

    echo "$(date) -- $id alignment done"
done

##########WES数据的比对
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
bwa=/mnt/data/wangdanke/Programes/bwa-0.7.17/bwa
samtools=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env/bin/samtools
# 比对（推荐bwa）

$bwa mem ref.fa sample_R1.fastq.gz sample_R2.fastq.gz | samtools sort -o sample.bam
$samtools index sample.bam



