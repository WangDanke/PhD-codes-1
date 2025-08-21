#做可变剪切之前的准备
#安装Spladder
pip install spladder
# 下载基因注释文件
wget -c https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gzip -d Homo_sapiens.GRCh38.110.gtf.gz

###1.为样本比对后的bam文件构建一个索引，就是每一个bam文件的具体位置，脚本
#!/bin/bash
set -e 
cat ~/PSCC_RNAseq/OSCC/GSE186775/GSE186775.txt | while read id
do
echo ~/PSCC_RNAseq/OSCC/GSE186775/Align_1st/$id/${id}_align.bam >> index_location.txt
done 

##每次到对应的数据集中去运行脚本
bash create_index.sh

##2.为每一个bam文件构建.bai文件
/mnt/data/wangdanke/Programes/samtools-1.19/samtools index -M ./fq_del_unmap/*/*_del.bam

##3.做可变剪切的脚本
#!/bin/bash
set -e

# path to annotation in GTF format
anno=/mnt/data/wangdanke/PSCC_RNAseq/reference_genome/Homo_sapiens.GRCh38.110.gtf
# path to the spladder executable
spladder=/mnt/data/wangdanke/anaconda3/bin/spladder
## text file containing the absolute paths to the alignment files in BAM format
# alignments.txt

$spladder build -o ./Spladder_out -a $anno -b  index_location.txt -v -c 2 --merge-strat merge_graphs 
done

####参数解析：-c 2,置信度，3是最严格的，这里取2，借鉴那篇《Cancer cell》

###单个样本
spladder build -o Alt_splice/ -b /mnt/data/wangdanke/PSCC_RNAseq/OSCC/GSE176221/Align_1st/SRR14740626/SRR14740626_align.bam -a ~/PSCC_RNAseq/Homo_sapiens.GRCh38.110.gtf



###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###网上借鉴
cat ../case | while read id
do
spladder build -o ./ -a ~/reference/gtf/hg38.gtf -b ../4.bam/${id}.bam --merge-strat single --no-extract-ase
echo ../4.bam/${id}.bam >> alignments.txt
done
sed -i ':t;N;s/\n/,/;b t' alignments.txt

spladder build -o ./ -a ~/reference/gtf/hg38.gtf -b `cat alignments.txt` --merge-strat merge_graphs --no-extract-ase

cat ../case | while read id
do
spladder build -o ./ -a ~/reference/gtf/hg38.gtf -b ../4.bam/${id}.bam --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode single
done

spladder build -o ./ -a ~/reference/gtf/hg38.gtf -b `cat alignments.txt` --merge-strat merge_graphs --no-extract-ase --quantify-graph --qmode collect

for type in exon_skip intron_retention alt_3prime alt_5prime mult_exon_skip mutex_exons
do
spladder build -o ./ -a ~/reference/gtf/hg38.gtf -b `cat alignments.txt` --event-types ${type}
done




