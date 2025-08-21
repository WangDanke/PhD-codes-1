##!/bin/bash
##this is a script used to dowload data from GEO database
if [ -z "$1" ]
then 
echo "please provide a GSE ID"
exit $1 
fi 
##创建文件夹，以GSE编号命名
mkdir $1
##递归下载该GSE号下的有数据
wget -c -nH ftp://ftp.ncbi.nlm.nih.gov/geo/series/${1: 0:-3}nnn/$1/matrix/