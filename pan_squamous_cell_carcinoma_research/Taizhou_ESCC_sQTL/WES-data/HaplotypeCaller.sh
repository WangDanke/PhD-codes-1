###使用bwa比对和gatk进行变异检测
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
export JAVA_HOME=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env
export PATH=$JAVA_HOME/bin:$PATH
export path=$TMPDIR
set -e
 
gatk=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env/bin/gatk
cat ~/PSCC_RNAseq/ESCC_Taizhou/ESCC_WES/sample_list_WES | while read id
do
##SNP calling
$gatk HaplotypeCaller \
    -R /mnt/data/Biodatabase/ref_human/hg38/variant_calling/Homo_sapiens_assembly38.fasta \
    -I /mnt/data/hbyuan/project/TaixinESCC_45/WES_run1/5.gatk/"$id".bqsr.bam \
    -O ~/PSCC_RNAseq/ESCC_taizhou/ESCC_WES/haplotypeCaller/"$id".g.vcf \
    -ERC GVCF
done

