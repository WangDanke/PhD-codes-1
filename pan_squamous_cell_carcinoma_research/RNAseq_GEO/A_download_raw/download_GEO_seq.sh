##!/bin/bash
###方法一：将SRR编号存在一个文件里，利用循环批量下载 只针对SRR后有七个数字的
###下载
cat GSE87410.txt|while read id ###GSE87410.txt记载了SRR编号
do
x=$(echo $id | cut -b1-6)
y=$(echo $id | cut -b10-10)
echo $id
ascp -QT -l 500m -P33001  -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/$x/00$y/$id/ ./
done

####运行时
nohup bash dowloading_GEO_seq.sh >log_dowloading_GEO_seq 2>&1 &

####方法二：直接利用下载地址一个一个的下载
#####一个非常直接的办法，直接将每一个下载地址都写一条命令，放在一个脚本里，去下载数据。
ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR036/SRR036752/SRR036752.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR036/SRR036756/SRR036756.fastq.gz  ./


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR036/SRR036753/SRR036753.fastq.gz  ./


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR036/SRR036757/SRR036757.fastq.gz  ./



ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/S10/SRR036754/SRR036754.fastq.gz  ./



####方法一改进：循环下载，自行判断SRR后面有几个数字编号
echo -e "\n\n\n 000# ascp Download FASTQ !!! \n\n\n"
date
## 目录设置
mkdir -p  ~/PSCC_RNAseq/GSE20116/fq/
cd  ~/PSCC_RNAseq/GSE20116/fq/
pwd


## 密匙路径
openssh=~/.aspera/connect/etc/asperaweb_id_dsa.openssh

cat ~/PSCC_RNAseq/GSE20116.txt | while read id
do
num=`echo $id | wc -m ` #这里注意，wc -m 会把每行结尾$也算为一个字符
echo "SRRnum+1 is $num" #因此SRRnum + 1才等于$num值
#如果样本编号为SRR+8位数 #
if [ $num -eq 12 ]
then
        echo "SRR + 8"
        x=$(echo $id | cut -b 1-6)
        y=$(echo $id | cut -b 10-11)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001  -k 1 -i $openssh \
                era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/0$y/$id/   ./ & )
#如果样本编号为SRR+7位数 #
elif [ $num -eq 11 ]
then
        echo  "SRR + 7"
        x=$(echo $id | cut -b 1-6)
        y=$(echo $id | cut -b 10-10)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001  -k 1 -i $openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/00$y/$id/   ./ & )
#如果样本编号为SRR+6位数 #
elif [ $num -eq 10 ]
then
        echo  "SRR + 6"
        x=$(echo $id |cut -b 1-6)
        echo "Downloading $id "
        ( ascp  -QT -l 500m -P33001 -k 1 -i  $openssh \
                        era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$x/$id/   ./ & )
fi
done

####后台运行脚本，并且把日志输出到log_download里面
nohup bash download.sh >log_download 2>&1 &


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/041/SRR10357041/SRR10357041_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/041/SRR10357041/SRR10357041_2.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/041/SRR10357041/SRR10357042_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/041/SRR10357041/SRR10357042_2.fastq.gz  ./



###发现有文件受损，于是重新下载一下，看是不是没什么问题
#GSE87410
ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR429/004/SRR4299084/SRR4299084_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR429/005/SRR4299085/SRR4299085_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/008/SRR4427168/SRR4427168_2.fastq.gz  ./           #下载好了20231212

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/002/SRR4427172/SRR4427172_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/004/SRR4427184/SRR4427184_1.fastq.gz  ./           ##下载好了20231213

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/005/SRR4427185/SRR4427185_1.fastq.gz  ./ 

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/006/SRR4427186/SRR4427186_1.fastq.gz  ./ 

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/006/SRR4427186/SRR4427186_2.fastq.gz  ./ 

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/003/SRR4427173/SRR4427173_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/003/SRR4427173/SRR4427173_2.fastq.gz  ./


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/000/SRR4427200/SRR4427200_1.fastq.gz  ./    #

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/005/SRR4427205/SRR4427205_1.fastq.gz  ./   # 下载好了20231213

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/006/SRR4427206/SRR4427206_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/000/SRR4427210/SRR4427210_1.fastq.gz  ./    #

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/001/SRR4427211/SRR4427211_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/008/SRR4427218/SRR4427218_2.fastq.gz  ./    #

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR442/002/SRR4427222/SRR4427222_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR488/002/SRR4880462/SRR4880462_1.fastq.gz  ./  ##

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR488/003/SRR4880463/SRR4880463_1.fastq.gz  ./  #

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR488/005/SRR4880465/SRR4880465_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR488/006/SRR4880466/SRR4880466_1.fastq.gz  ./ #

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR488/007/SRR4880467/SRR4880467_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR488/008/SRR4880468/SRR4880468_1.fastq.gz  ./





#GSE186775
ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR166/022/SRR16612822/SRR16612822_2.fastq.gz  ./      #下载好了

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR166/027/SRR16612827/SRR16612827_2.fastq.gz  ./    #下载好了，至此GSE186775数据全部下载完毕，并且数据完整


#GSE184616

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR160/047/SRR16013047/SRR16013047_1.fastq.gz  ./  #


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR160/050/SRR16013050/SRR16013050_1.fastq.gz  ./ #


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR160/053/SRR16013053/SRR16013053_1.fastq.gz  ./ #


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR160/056/SRR16013056/SRR16013056_1.fastq.gz  ./


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR160/063/SRR16013063/SRR16013063_1.fastq.gz  ./ #


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR160/086/SRR16013086/SRR16013086_2.fastq.gz  ./


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR160/090/SRR16013090/SRR16013090_1.fastq.gz  ./

######################cSCC
##GSE139505

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/041/SRR10357041/SRR10357041_1.fastq.gz  ./ #

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/041/SRR10357041/SRR10357041_2.fastq.gz  ./ #



ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/031/SRR10357031/SRR10357031_2.fastq.gz  ./


ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/041/SRR10357041/SRR10357041_1.fastq.gz  ./

##
ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/042/SRR10357042/SRR10357042_2.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR103/042/SRR10357042/SRR10357042_1.fastq.gz  ./

##GSE191335
ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR172/030/SRR17283430/SRR17283430_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR172/031/SRR17283431/SRR17283431_1.fastq.gz  ./

##GSE84293

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR387/003/SRR3877303/SRR3877303_1.fastq.gz  ./ #

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR387/008/SRR3877308/SRR3877308_1.fastq.gz  ./ #

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR387/004/SRR3877314/SRR3877314_1.fastq.gz  ./

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR387/005/SRR3877315/SRR3877315_1.fastq.gz  ./ 

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR387/005/SRR3877315/SRR3877315_2.fastq.gz  ./

###LUSC
#GSE159857

ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR128/078/SRR12879178/SRR12879178_2.fastq.gz  ./

#GSE179879
ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR150/051/SRR15097651/SRR15097651_2.fastq.gz  ./


###Cervical_SCC
##GSE144293
ascp  -vQT -l 500m -P33001 -k 1 -i \
~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR109/093/SRR10969793/SRR10969793_1.fastq.gz  ./

