import scanpy as sc  
import pandas as pd  
import numpy as np  
  
# 假设你有多个 matrix 和 barcodes 文件，这里以列表的形式给出它们的文件名  
matrix_files = ['matrix1.csv', 'matrix2.csv', 'matrix3.csv']  # 根据实际情况修改文件名  
barcodes_files = ['barcodes1.tsv', 'barcodes2.tsv', 'barcodes3.tsv']  # 根据实际情况修改文件名  
batch_keys = ['batch1', 'batch2', 'batch3']  # 对应的批次标签，长度应与 matrix_files 相同  
  
# 初始化一个空的 AnnData 对象  
adata = sc.AnnData()  
  
# 循环遍历所有的 matrix 和 barcodes 文件，并将它们整合到 AnnData 对象中  
for matrix_file, barcodes_file, batch_key in zip(matrix_files, barcodes_files, batch_keys):  
    # 读取 matrix 文件（假设第一列是基因名，其他列是表达值）  
    matrix = pd.read_csv(matrix_file, index_col=0)  
      
    # 读取 barcodes 文件  
    barcodes = pd.read_csv(barcodes_file, header=None, names=['barcode'])  
      
    # 添加一个批次标签列  
    barcodes['batch'] = batch_key  
      
    # 创建一个新的 AnnData 对象（临时）  
    adata_tmp = sc.AnnData(X=matrix.values, obs=barcodes, var=pd.DataFrame(index=matrix.index))  
      
    # 如果 adata 是空的，则直接用 adata_tmp 初始化它  
    if adata.shape == (0, 0):  
        adata = adata_tmp  
    else:  
        # 否则，将 adata_tmp 与 adata 合并  
        adata = adata.concatenate(*[adata, adata_tmp], axis=0)  
  
# 确保 obs（观测值/样本）的索引是连续的，并重置批次标签为分类变量  
adata.obs_names = np.arange(adata.n_obs).astype(str)  
adata.obs['batch'] = adata.obs['batch'].astype('category')  
  
# 数据预处理：归一化和方差稳定化  
sc.pp.normalize_total(adata, target_sum=1e6)  # 归一化到总和为 1e6  
sc.pp.log1p(adata)  # 对数变换（加 1 后取对数）  
sc.pp.pca(adata, svd_solver='arpack')  # PCA 降维，用于后续分析  
  
# 使用 bbknn 减少批次效应  
sc.external.pp.bbknn(adata, batch_key='batch', n_pcs=50)  # 使用前 50 个 PCA 成分来计算 bbknn 图  
  
# 你可以选择进一步处理 adata，例如聚类、降维可视化等  
# ...  
  
# 最后，将整合并处理后的 AnnData 对象保存为 .h5ad 格式的文件  
adata.write('merged_data_corrected.h5ad')  
  
print("Merged and batch-corrected data has been saved to merged_data_corrected.h5ad")