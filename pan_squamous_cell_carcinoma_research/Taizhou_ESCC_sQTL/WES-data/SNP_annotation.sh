##对SNP进行注释
#第一种方法
##使用annovar对SNP进行注释

#进入目录存放工具 
cd ~/Programes

##下载ANNOVAR,需要现在官网注册 随后才会发给你下载链接，估计每个人都不一样

wget -c http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz


# 下载基因注释数据库（RefGene）
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/

# 下载常用变异数据库（dbSNP）
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 humandb/

# 下载 ClinVar 数据库（临床意义）,这一步当时下载失败了 所以放弃了 后续没有使用这个数据库
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240318 humandb/

# 下载 gnomAD 频率库
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad211_exome humandb/

##进到这个annovar的文件夹下，进行注释
#直接进行功能注释
perl table_annovar.pl tumor.snps.withid.vcf.gz humandb/ \
    -buildver hg38 \
    -out result_annotion_tumor \
    -remove \
    -protocol refGene,avsnp150,gnomad211_exome \
    -operation g,f,f \
    -nastring . \
    -vcfinput


