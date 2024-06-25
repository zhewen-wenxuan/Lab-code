import numpy as np
import anndata as ad
import scanpy as sc

# 加载 .npy 文件
expression_matrix = np.load('/path/to/expression_matrix.npy')
genes = np.load('/path/to/genes.npy')
barcodes = np.load('/path/to/barcodes.npy')

# 创建 AnnData 对象
adata = ad.AnnData(X=expression_matrix, var={'gene_names': genes}, obs={'cell_barcodes': barcodes})

# 计算质量控制指标
sc.pp.calculate_qc_metrics(adata, inplace=True)

# 过滤掉低于一定阈值的基因和细胞
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 对数转换和标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 选择高度变异的基因
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']]

# 数据缩放
sc.pp.scale(adata, max_value=10)

# PCA降维
sc.tl.pca(adata, svd_solver='arpack')

# 可视化 PCA 结果
sc.pl.pca(adata, color='gene_names')

# 保存处理后的数据
adata.write('/path/to/processed_expression_matrix.h5ad')
