import numpy as np
import anndata as ad
import scanpy as sc

# 加载 .npy 文件
expression_matrix = np.load('/path/to/expression_matrix.npy')
genes = np.load('/path/to/genes.npy')
barcodes = np.load('/path/to/barcodes.npy')

# 创建 AnnData 对象
adata = ad.AnnData(X=expression_matrix, var={'gene_names': genes}, obs={'cell_barcodes': barcodes})

# 查看 AnnData 对象
print(adata)

# 保存为 .h5ad 文件
adata.write('/path/to/expression_matrix.h5ad')

# 以后可以直接加载 .h5ad 文件
adata = sc.read_h5ad('/path/to/expression_matrix.h5ad')

# 查看数据
print(adata)
