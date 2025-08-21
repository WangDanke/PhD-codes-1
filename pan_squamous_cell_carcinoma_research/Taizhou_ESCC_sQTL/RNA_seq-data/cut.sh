#!/bin/sh
#$ -S /bin/sh
#$ -N cutadapt      #任务名
#$ -V
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./                            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q normal.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq normal.q@ibnode13    # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 4-4        #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR
set -e

cat ~/PSCC_RNAseq/ESCC_taizhou/sample_list | while read id
do
mkdir ./fq_cut10/$id/
cutadapt -u 10 -o ./fq_cut10/$id/"$id"_R1_cut.fq.gz /mnt/data/hbyuan/project/TaixinESCC_45/RNA/2.clean_fq/"$id"_R1_val_1.fq.gz
cutadapt -u 10 -o ./fq_cut10/$id/"$id"_R2_cut.fq.gz /mnt/data/hbyuan/project/TaixinESCC_45/RNA/2.clean_fq/"$id"_R2_val_2.fq.gz
done




