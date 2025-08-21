#创建这个文件的初衷就是希望将自己写的每一条垃圾命令和代码都详细的记录下来，本人记性不好。


#进入目标文件夹
cd ./目标文件夹/ #其中./代表的是当前文件夹（相对位置）

##查看当前文件夹的内容
ls -lh ## -l 显示文件的长名字信息，-h 以人类能看懂的方式显示
ll -h ./文件夹1/子文件夹/*/* ##显示子文件夹中的所有文件夹中的所有文件


###目录操作
mkdir ./新目录 # 在当前文件夹下创建新目录
rmdir ./目录/文件 #移除空目录或者文件
pwd #显示当前文目录
cp -rp ./原位置原目录 ./新位置原目录 #复制文件或目录 -r 代表的是复制目录，如果是文件就可以不加；-p 保留文件的属性
mv ./原位置原目录 ./新位置原目录 #剪切文件或目录 
rm -rf #删除文件或者目录 -r 代表的是删除目录（必须加）；-f 代表的是删除文件（强制删除）




#脚本命令
#bash 脚本开头加上
#！/bin/bash 告诉服务器这是一个bash脚本
set -e ##放在脚本开头，遇到错误脚本便会停止运行


#查看当前的目录各文件夹占用多少内存
du -sh ./*

##杀死相关进程的代码
killall -9 ascp ##杀死commond名称为ascp的所有进程
kill PID ##杀死特定PID的进程

###终止脚本运行
#首先以脚本名称查找进程
ps -ef | grep cut14.sh
##其中grep --color-auto是查找的进程本身，另外的则是你要找的脚本进程



###下载和安装ascp
##下载ascp 和解压
wget -c https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_4.2.7.445_linux_x86_64.tar.gz
tar -vaf ibm-aspera-connect_4.2.7.445_linux_x86_64.tar.gz
##运行安装
bash ibm-aspera-connect-3.8.3.170430-linux-g2.12-64.sh
#添加环境变量
export PATH = $PATH:/mnt/data/wangdanke/.aspera/connect/bin
##验证是否安装成功
ascp --help


###直接使用conda安装ascp
conda install -c rpetit3 aspera-connect



##同步两个服务器之间的文件和文件夹
scp -r -P 22 -l 2048 wangdanke@122.194.116.207:/share/home/wangdanke/PSCC_RNAseq/OSCC/GSE176221/fq/SRR14740626/ /mnt/data/wangdanke/PSCC_RNAseq/OSCC/GSE176221/fq/ 
##-r:递归该文件夹下的所有文件和文件夹；-P：老服务器的端口号;-l 2048:限速2048KB传输

rsync -av wangdanke@122.194.116.207:/share/home/wangdanke/PSCC_RNAseq/OSCC/GSE186775/ /mnt/data/wangdanke/PSCC_RNAseq/OSCC/ ##这个的优势是可以断点续传，但是代码是否正确仍不清楚



##配置环境变量（一件我觉得很难的事情，所以在最开始的时候我失败过很多次
##下载和配置fastqc的环境变量
wget -c 下载的链接地址 #下载
unzip fastqc_v0.12.1.zip #解压
cd ./FastQC/ #到达刚解压之后的文件夹
chmod u+x fastqc #给fastqc赋予执行权限
##配置环境
vim ~/.bashrc #进入记载环境变量的文件
export PATH="mnt/data/wangdanke/Programes/FastQC:$PATH" #将此行代码加入.bashrc文件中（大概的意思是，fastqc所在的位置，加入了环境变量，于是可以随时调用）
source ~./bashrc #使环境变量生效
fastqc --help ##查看是否配置成功


###安装multiqc
pip install multiqc 

##2024年1月2日
#由于之前重新安装了anaconda，于是需要重新安装multiqc,
pip install multiqc #发现有版本对不上的问题，于是将pip升级到更高的版本，再运行此代码
pip install --upgrade pip


#2023年11月13日
#这几天突然遇到了不能使用conda命令的情况，发现也是直接添加环境变量就可以搞定的
vim ~/.bashrc #发现里面已经有conda的位置环境变量
source ~/.bashrc #激活虚拟环境


####虚拟终端：虚拟终端中在执行的命令不会随着你关闭终端软件（MOBAXterm）而中断
##创建
screen -R Download #构建了一个名为Download的虚拟终端
###退出
screen -R Download #去到Download 终端
exit #先进去，再exit 即可退出该虚拟终端
screen -ls #列出当前的所有虚拟终端
##对于意外退出，且仍保持Attached的状态的虚拟终端，要先踢掉前一次的登录，之后再登录
#踢掉
screen -D -r GSE139505
#再进入
screen -R GSE139505

#  新建screen
screen -S your_screen_name 
# 进入screen
screen -r your_screen_name
 
Ctrl+D  # 在当前screen下，输入Ctrl+D，删除该screen
Ctrl+A，Ctrl+D  # 在当前screen下，输入先后Ctrl+A，Ctrl+D，退出该screen 
#  显示screen list
​​​​​​​screen -ls
# 连接状态为【Attached】的screen
screen -D  -r your_screen_name  # 解释：-D -r 先踢掉前一用户，再登陆
# 判断当前是否在screen中断下,Ubuntu系统,可以这样:
sudo vim /etc/screenrc
# 文件末尾追加一行即可允许设置screen标题
caption always "%{.bW}%-w%{.rW}%n %t%{-}%+w %=%H %Y/%m/%d "
# 删除指定screen, your_screen_name为待删除的screen name
​​​​​​​screen -S your_screen_name -X quit



####定时执行脚本或程序软件等
crontab -e #打开定时配置文件
###在配置文件中写入定时任务操作
0 12 * * * source ~/.bashrc; sh ~/PSCC_RNAseq/HNSC/download.sh > download.log #每天的12点0分执行~/mnt/data/wangdanke/PSCC_RNAseq/HNSC/文件夹下的download.sh 脚本
#并且把日志写入download.log中
：wq ##保存并退出即可
#写入定时任务的基本格式是：
* * * * * command
分 时 日 月 周
*/1 * * * *  sh test.sh #表示每分钟执行一次任务test.sh 脚本
0 1,2,3,4,5,6 * * * sh test.sh #表示每天1点0分，2点0分，3点0分......各执行一次test.sh 脚本

#  ~/表示的就是本人的主目录；./表示当前目录


##创建和配置R环境
conda create -n R-4.3
conda activate R-4.3
conda install r-base=4.3.0


##统计数量
ls | wc -w #查看有多少个文件及文件夹
ls | wc -c #查看目录下有多少个字节数
ls -l |grep "^-"|wc -l #统计某文件夹下文件的个数
ls -l |grep "^-"|wc -l #统计某文件夹下目录的个数
ls -lR|grep "^-"|wc -l #统计文件夹下文件的个数，包括子文件夹里的


#安装trim-galore
pip install trim-galore
#配套使用cutadapt
pip install cutadapt








###转移特定名称的文件夹
#首先构建你需要的文件夹名字的list
$ cat >idname
SRR4299084
SRR4427172
SRR4427185
SRR4427195
SRR4427197
SRR4427205
SRR4427207
SRR4427208
SRR4427209
SRR4427211
SRR4427212
SRR4427217
SRR4427222
SRR4880462
SRR4880464
SRR4880465
SRR4880466
SRR4427198
SRR4299085
SRR4427168
SRR4427184
SRR4427186
SRR4427206
SRR4427210
SRR4427218
SRR4427221
SRR4880467
SRR4880468
SRR4427169
SRR4427173
SRR4427196
SRR4427199
SRR4427200
SRR4427213
SRR4427214
SRR4880463
^C 
#运行脚本
$ bash mv.sh


##!/bin/bash
set -e 
#循环转移脚本mv.sh
cat ~/PSCC_RNAseq/LUSC/GSE87410_LUSC/idname_LUSC | while read id
do
mv ~/PSCC_RNAseq/GSE87410/fq_del_unmap/$id/ ~/PSCC_RNAseq//LUSC/GSE87410_LUSC/fq_del_unmap/
done


##下载和配置miniconda 
#使用wget 下载
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh



###使用调度系统，在需要较大的算力的时候，可以使用ibfat1、ibfat2、ibfat3 节点，脚本需要的格式和环境配置 #187-200行

#!/bin/sh 
#$ -S /bin/sh
#$ -N GSE186775      #任务名
#$ -V                
#$ -j y                     #任务的stdout信息和stderr信息的文件名一样
#$ -o ./                    #定义调度系统stdout文件存放目录
#$ -e ./		            #定义调度系统stderr文件存放目录
#$ -cwd                     #在节点上运行任务时，会先进入执行qsub时的路径
#$ -q normal.q                 #选项用于指定作业提交到哪个队列中。在这种情况下，作业将被提交到名为 fat.q 的队列中
#$ -masterq normal.q@ibnode12    # 选项用于指定作业提交到哪个主机上。在这种情况下，作业将被提交到名为 ibfat1 的主机上
#$ -pe thread 4-4        #mpi是并行环境名，64-64表示申请64核
source ~/.bashrc
hash -r
export path=$TMPDIR
set -e
cat ~/PSCC_RNAseq/OSCC/GSE186775/idname | while read id
do
mkdir ./fq_qc_cutadapt/$id/
cutadapt -u 14 -o ./fq_qc_cutadapt/$id/"$id"_1_cut14.fq.gz ./fq_qc/$id/*_1.fq.gz
cutadapt -u 14 -o ./fq_qc_cutadapt/$id/"$id"_2_cut14.fq.gz ./fq_qc/$id/*_2.fq.gz
done

##切换节点
ssh ibfat1/ibfat3 #password:fudanwangdanke

##为服务器配置conda 环境，之前使用的是miniconda，觉得不是很好，于是卸载了miniconda，重新安装了anaconda（anaconda内部就配置了python，因此不需要本人再重新下载python了
#删除miniconda
rm -rvf miniconda3
#删除相关的环境变量配置
#编辑.bashrc ,将其中与conda相关的环境变量配置都删除掉
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/opt/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/opt/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
#并且删除相关的文件
rm -rfv ~/.condarc ~/.conda ~/.continuum

#重新下载anaconda，不要下载最新版本，因为可能会不是很稳定，吸取之前的教训。此次下载的是2020年的anaconda
wget -c https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
#运行安装脚本,安装过程中一直yes、enter就行
chmod +x Ananconda3-2020.11-linux-x86_64.sh
bash Ananconda3-2020.11-linux-x86_64.sh
#激活环境变量
source ~/.bashrc

####在conda 中安装R语言
conda info --envs  # 查看目前的conda环境
conda search r-base #查看当前conda 适配的R版本
conda create -n R4.3  # 创建名为R3.5的环境
source activate R4.3   #激活R3.5环境
conda install r-base=4.3.1  #安装R 指定为R版本
conda deactivate # 退出当前环境
conda remove --name R3.5 --all #移除conda 的虚拟环境

##安装hisat2
conda install -c bioconda hisat2


####修改文件名
cat ~/PSCC_RNAseq/OSCC/GSE176221/GSE176221.txt | while read id
do
mv ./fq_qc_rm_rRNA/$id/*.fq.2.gz ./fq_qc_rm_rRNA/$id/"$id"_2.fq.gz
mv ./fq_qc_rm_rRNA/$id/*.fq.1.gz ./fq_qc_rm_rRNA/$id/"$id"_1.fq.gz
done 

####修改文件名
cat ~/PSCC_RNAseq/OSCC/GSE186775/GSE186775.txt | while read id
do
mv ./fq_fastqc/after_rmrRNA/*_rRNAremoved.fq.2_fastqc.html ./fq_fastqc/after_rmrRNA/mv_2/"$id"_2_fastqc.html
mv ./fq_fastqc/after_rmrRNA/*_rRNAremoved.fq.2_fastqc.zip ./fq_fastqc/after_rmrRNA/mv_2/"$id"_2_fastqc.zip
mv ./fq_fastqc/after_rmrRNA/*_rRNAremoved.fq.1_fastqc.zip ./fq_fastqc/after_rmrRNA/mv_2/"$id"_1_fastqc.zip
mv ./fq_fastqc/after_rmrRNA/*_rRNAremoved.fq.1_fastqc.html ./fq_fastqc/after_rmrRNA/mv_2/"$id"_1_fastqc.html
done



###在下载软件之前，默认处于base环境下，最好使用conda新建一个小环境，防止软件因为依赖关系而无法安装。
# 查看conda中所有的环境
conda env list
# 查看当前环境中安装的软件
conda list

# 新建环境，python版本可以指定
conda create -n rna -y python=3.6
# 激活环境
conda activate rna

#以下是Jimmy老师2018年推荐安装的生信软件：
#质控：fastqc、trimmomatic、cutadapt、trim-galore
#比对：star、hisat2、bowtie2、bwa、subread
#计数：htseq、bedtools、deeptools、salmon
##从这个角度看，后面分析可以新建一个专门做RNA-seq分析的小环境，这样可以帮助更好地安装软件等


##########转移具有特定名称的文件夹

cat  ~/PSCC_RNAseq/LUSC/GSE159857/GSE159857_SCC_list | while read id
do
mv ./fq/$id/ ./fq_sub_SCC/
done


#####输出具有特定名称的文件的路径
data_dir="./Align_hisat2/"
output="bam_paths.txt"

> "$output"  # 清空输出文件
while read id; do
    if [ -d "$data_dir/$id" ]; then
        find "$data_dir/$id" -type f -name "*.bam" >> "$output"
    fi
done < sample_list