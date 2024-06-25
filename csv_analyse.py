# 1.安装包
pip install pandas scanpy anndata matplotlib seaborn

# 2.数据导入
# 假设你有以下两个CSV文件：
# gene_expression.csv：基因表达矩阵，行是细胞，列是基因。
# metadata.csv：元数据，包含每个细胞的附加信息（如坐标等）。
import pandas as pd
import scanpy as sc

# 读取基因表达数据
expression_data = pd.read_csv('path/to/gene_expression.csv', index_col=0)

# 读取元数据
metadata = pd.read_csv('path/to/metadata.csv', index_col=0)

# 创建AnnData对象
adata = sc.AnnData(X=expression_data.values)
adata.obs = metadata
adata.var['gene_symbols'] = expression_data.columns
adata.var_names = adata.var['gene_symbols']

# 查看数据概况
print(adata)

# 3.质量控制
# 计算线粒体基因比例
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 绘制质量控制图
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# 过滤低质量细胞
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# 过滤低表达基因
sc.pp.filter_genes(adata, min_cells=3)

# 4.数据标准化和高变基因选择
# 数据标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 识别高变基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# 数据缩放
sc.pp.scale(adata, max_value=10)

# 5.PCA降维与聚类分析
# PCA降维
sc.tl.pca(adata, svd_solver='arpack')

# 绘制PCA图
sc.pl.pca(adata, color='CST3')

# 聚类分析
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata)

# UMAP降维
sc.tl.umap(adata)

# 绘制UMAP图
sc.pl.umap(adata, color=['leiden', 'CST3'])

# 6.数据注释与可视化
# 找到聚类标记基因
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

# 绘制聚类标记基因
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# 查看某个基因在空间上的表达
sc.pl.spatial(adata, img_key='tissue_hires', color=['CST3', 'leiden'])

# 7.功能分析
import gseapy as gp
from gseapy.plot import barplot, dotplot

# 对某个聚类的标记基因进行GO富集分析
cluster_genes = adata.uns['rank_genes_groups']['names'][0]
go_results = gp.enrichr(gene_list=cluster_genes, gene_sets='GO_Biological_Process_2021', organism='Human')

# 绘制富集分析结果
barplot(go_results.res2d, title='GO_Biological_Process_2021', ofname='go_enrichment.png')
