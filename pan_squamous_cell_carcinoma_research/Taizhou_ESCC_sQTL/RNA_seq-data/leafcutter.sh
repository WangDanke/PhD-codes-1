##使用leafcutter做可变剪接分析

#!/bin/sh 
#$ -S /bin/sh
#$ -N leafcutter      #任务名
#$ -V                
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./		            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q normal.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq normal.q@ibnode14    # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 8-8        #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR
set -e

samtools=/mnt/data/wangdanke/anaconda3/envs/s_qtl_env/bin/samtools
regtools=/mnt/data/wangdanke/Programes/regtools/build/regtools 

cat bam_paths.txt |while read id
do
file=$(basename $id )
sample=${file%%.*}
    echo Converting $id to ./leafCutter/$sample.junc
    $samtools index $id
    $regtools junctions extract -a 8 -m 50 -s 1 -M 500000 $id -o ./leafCutter/$sample.junc
    echo $file.junc >> test_juncfiles.txt
done


#step2:intron clustering

#!/bin/sh 
#$ -S /bin/sh
#$ -N clustering      #任务名
#$ -V                
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./		            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q normal.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq normal.q@ibnode14    # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 8-8        #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR
set -e
leafcutter_cluster_regtools=/mnt/data/wangdanke/Programes/leafcutter/clustering/leafcutter_cluster_regtools.py
python $leafcutter_cluster_regtools -j test_juncfiles.txt -m 50 -o testYRIvsEU -l 500000


##计算PCA，前十个主成分，准备各个染色体的标准化后的可变剪接事件,用于后续的sQTL分析，分成肿瘤组织和癌旁组织分别计算

python /mnt/data/wangdanke/leafcutter/scripts/prepare_phenotype_table.py ./testYRIvsEU_perind.counts.gz -p 10




##################################################################差异可变剪接，使用的是R脚本 因此必须要安装R包 leafcutter, 但是我安装了很久都失败了
##安装了大半个月的leafcutter R包 都失败了，于是借用佳琪师姐的服务器账号，跑了差异可变剪接 后面可能还要借用可视化
/home/public/myspace/jqzhou/leafcutter/scripts/leafcutter_ds.R --num_threads 4 /home/public/myspace/jqzhou/DK_splicing/GSE164158/testYRIvsEU_perind_numers.counts.gz /home/public/myspace/jqzhou/DK_splicing/GSE164158/sample_group.txt


###注释差异可变剪接的结果
Rscript ./prepare_ds_results.R \
  -o my_data_leafviz.RData \
  -m sample_group.txt \
  -f 0.05 \
  -c /home/public/myspace/jqzhou/leafcutter/scripts/leafcutter_ds.R \
  testYRIvsEU_perind_numers.counts \
  leafcutter_ds_cluster_significance.txt \
  leafcutter_ds_effect_sizes.txt \
  annotation_codes/gencode_hg38/gencode_hg38

Rscript ./prepare_ds_results.R \
  -o GSE149609_leafviz.RData \
  -m sample_group.txt \
  -f 0.05 \
  -c /home/public/myspace/jqzhou/leafcutter/scripts/leafcutter_ds.R \
  testYRIvsEU_perind_numers.counts.gz \
  leafcutter_ds_cluster_significance.txt \
  leafcutter_ds_effect_sizes.txt \
  /home/public/myspace/jqzhou/DK_splicing/annotation_codes/gencode_hg38/gencode_hg38


##可视化
Rscript run_leafviz.R my_data_leafviz.RData 






