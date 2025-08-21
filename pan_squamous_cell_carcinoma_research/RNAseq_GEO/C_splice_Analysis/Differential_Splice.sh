#step 1 :为每一个独立数据集构建可变剪切图谱
#step 2 ：每一个单独的数据集进行差异可变剪切
#step 3 : 对于单个癌种的差异可变剪切取并集，



##step 2 具体步骤，建立在step 1之上
###分别为normal样本和tumor样本构建index of location 
bash index.sh 

#####跑差异剪切的脚本
#!/bin/sh 
#$ -S /bin/sh
#$ -N GSE139505      #任务名
#$ -V                
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./		            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q normal.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq normal.q@ibnode18     # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 28-28         #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR
set -e

spladder=/mnt/data/wangdanke/anaconda3/bin/spladder
## text file containing the absolute paths to the alignment files in BAM format
# alignments.txt

$spladder test --conditionA tumor_location.txt --conditionB normal_location.txt --outdir ./Spladder_out -c 2 --labelA Tumor --labelB Normal --diagnose-plots --plot-format pdf  --parallel 8
done
