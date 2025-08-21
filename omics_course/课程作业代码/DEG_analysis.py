#导入需要的包
import omicverse as ov 
import pandas as pd 
import numpy as np 
import scanpy as sc 
import matplotlib.pyplot as plt
import saeborn as sns 

##设定绘图格式
ov.utils.ov_plotset()

##下载基因集，当我们转换基因ID时，我们需要准备一个映射文件，在这里我们预处理了六个基因组gtf文件和生成的映射对，包括T2T-CHM13，GRCh38，GRCh37，GRCm39，danRer7和 danRer11。
# 如果需要转换其他ID，可以使用高铁费将文件放在genesets目录中生成立自己的映射
ov.utils.download_geneid_annotation_pair()

##读取数据##可以从线上的网站直接读取数据
data=pd.read_csv(https://raw.githubusercontent.com/Starlitnightly/ov/master/sample/counts.txt',index_col=0,sep='\t',header=1)
data.columns=[i.split('/'[-1]).replace('.bam','')for i in data.columns]
data.head()

###我们可以非常简单地通过omicverse进行差异表达基因分析，只需要提供一个表达式矩阵。我们首先创建一个 pyDEG 对象，并使用drop_duplicates_index去除重复的基因。
# 由于部分基因名相同，我们的去除保留了表达量最大的基因名
dds=ov.bulk.pyDEG(data)
dds.drop_duplicates_index()
print('... drop_duplicates_index success')

##我们还需要去除表达矩阵的批次效应 (batch effect)，我们使用DEseq2的的 SizeFactor 来对我们的矩阵计算归一化因子来去除批次效应。
dds.normalize()
print('... estimateSizeFactors and normalize success')

##现在我们可以从表达矩阵中计算差异表达基因，在计算前我们需要输入实验组和对照组。在这里，我们指定 4-3和4-4为实验组，1--1, 1--2为对照组，使用ttest进行差异表达分析计算。
# 当然你也可以使用wilcox来计算。此外deseq2也是支持的，不过流程可能会有一些区别，我们放到下一期讲。
treatment_groups=['4-3','4-4']
control_groups=['1--1','1--2']
result=dds.deg_analysis(treatment_groups,control_groups,methond='ttest')
result.head()

##在计算完差异表达基因后，我们会发现一个重要的事情，就是低表达基因有很多，如果我们不对其进行过滤，会影响后续火山图的绘制，我们设定基因的平均表达量大于1作为阈值，
# 将平均表达量低于1的基因全部过滤掉。
print(result.shape)
result=result.loc[result['log2(BaseMean)']>1]
print(result.shape)

#我们还需要设置 Foldchange 的阈值，我们准备了一个名为 foldchange_set 的方法函数来完成。此函数根据 log2FC 分布自动计算适当的阈值，但您也可以手动输入阈值。该函数有三个参数：
#fc_threshold: 差异表达倍数的阈值，-1为自动计算
#pval_threshold: 差异表达基因的p-value过滤值，默认为0.05，在有些情况下可以设定为0.1，意味着统计学差异不显著。
#logp_max: p值的最大值，由于部分p值过小，甚至为0，取对数后火山图绘制较为困难，我们可以设定一个上限，高于这个上限的p值全部统一。
# -1 means automatically calculates
dds.foldchange_set(fc_threshold=-1,
                   pval_threshold=0.05,
                   logp_max=6)
                   
#差异表达结果的可视化
# #omicverse除了有较为完善的分析能力外，还有极强的可视化能力。首先是火山图，我们使用 plot_volcano函数来实现。该函数可以绘制你感兴趣的基因或高表达的基因。您需要输入一些参数:
#title: 火山图的标题
#figsize: 图像大小
#plot_genes: 需要绘制的基因，格式为list。如['Gm8925','Snorc']
#plot_genes_num: 需要绘制的基因数，该参数与plot_genes互斥，如果我们没有指定需要绘制的基因，可以自动绘制前n个高差异表达倍数的基因。
#此外还可以制定绘制的颜色等，具体的参数可以使用help(dds.plot_volcano)

dds.plot_volcano(title('DEG analysis'),figsize=(4,4),plot_genes_num=8,plot_genes_fontsize=2,)

##还可以绘制特定基因的箱线图，我们可以使用plot_boxplot函数来完成该任务
dds.plotboxplot(genes=["Ckap2","Lef1"],treatment_groups=treatment_groups,control_groups=control_groups,figsize=(2,3),fontsize=12,legend_bbox=（2,0.55）)

#富集分析
#gseapy包被封进了omicverse,其中包含了GSEA的富集分析相关功能，类似的还有一些通路数据库，可以使用ov.utils.download_pathway_database()进行自动下载（里面有5个基因集），
#除此之外还可以在以下网站上下载自己感兴趣的基因集https://maayanlab.cloud/enrichr
ov.utils.download_pathway_database()
#读取通路基因集，我们读取Wiki通路数据库
pathway_dict=ov.utils.geneset_prepare('genesets/WikiPathways_2019_Mouse.txt',organism='Mouse')
#差异表达基因提取
deg_genes=dds.result.loc[dds.result['sig']!='normal'].index.tolist()
#通路富集分析
enr=ov.bulk.geneset_enrichment(gene_list=deg_genes,
                                pathways_dict=pathway_dict,
                                pvalue_type='auto',
                                organism='mouse')
##可视化
ov.bulk.geneset_plot(enr,figsize=(2,5),fig_title='Wiki Pathway enrichment',
                        cmap='Reds')
#如果需要保存的话,使用`plt.savefig`来保存图像
plt.savefig("enr_pathway.png",dpi=300,bbox_inches = 'tight')