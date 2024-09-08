import scanpy as sc  
import pandas as pd  
import numpy as np  
  
# 假设你有多个 matrix 和 barcodes 文件，这里以列表的形式给出它们的文件名  
matrix_files = ['matrix1.csv', 'matrix2.csv', 'matrix3.csv']  # 根据实际情况修改文件名  
barcodes_files = ['barcodes1.tsv', 'barcodes2.tsv', 'barcodes3.tsv']  # 根据实际情况修改文件名  
  
# 初始化一个空的 AnnData 对象  
adata = sc.AnnData()  
  
# 循环遍历所有的 matrix 和 barcodes 文件，并将它们整合到 AnnData 对象中  
for matrix_file, barcodes_file in zip(matrix_files, barcodes_files):  
    # 读取 matrix 文件（假设第一列是基因名，其他列是表达值）  
    matrix = pd.read_csv(matrix_file, index_col=0)  
      
    # 读取 barcodes 文件  
    barcodes = pd.read_csv(barcodes_file, header=None, names=['barcode'])  
      
    # 创建一个新的 AnnData 对象（临时）  
    adata_tmp = sc.AnnData(X=matrix.values, obs=barcodes, var=pd.DataFrame(index=matrix.index))  
      
    # 如果 adata 是空的，则直接用 adata_tmp 初始化它  
    if adata.shape == (0, 0):  
        adata = adata_tmp  
    else:  
        # 否则，将 adata_tmp 与 adata 合并  
        # 注意：这里我们假设所有 matrix 的基因顺序和名称都是一致的  
        adata = adata.concatenate(*[adata, adata_tmp], axis=0)  
  
# 在合并后的 AnnData 对象中，确保 obs（观测值/样本）的索引是连续的  
adata.obs_names = np.arange(adata.n_obs).astype(str)  
  
# 如果你想要根据某些标准对基因进行过滤（例如，去除表达量为0的基因），你可以在这里做  
# 例如：sc.pp.filter_genes(adata, min_counts=1)  # 去除在所有样本中表达量都为0的基因  
  
# 最后，将整合后的 AnnData 对象保存为 .h5ad 格式的文件  
adata.write('merged_data.h5ad')  
  
print("Merged data has been saved to merged_data.h5ad")