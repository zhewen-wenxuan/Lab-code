import scanpy as sc  
import pandas as pd  
import numpy as np  
from scipy.io import mmread  
  
# 设定mtx、barcodes及features文件的文件名（请根据实际情况修改）  
mtx_files = ['matrix1.mtx', 'matrix2.mtx', 'matrix3.mtx']  
barcodes_files = ['barcodes1.tsv', 'barcodes2.tsv', 'barcodes3.tsv']  
features_file = 'features.tsv'  
metadata_file = 'metadata.csv'  # 新增metadata文件  
  
# 读取features文件（我们假设第一列是特征名）  
features = pd.read_csv(features_file, sep='\t', header=None, names=['feature'], usecols=[0])  
  
# 初始化一个空的AnnData对象  
adata = sc.AnnData()  
  
# 遍历所有的mtx和barcodes文件，整合数据到AnnData对象中  
for mtx_file, barcodes_file in zip(mtx_files, barcodes_files):  
    # 读取mtx文件（稀疏矩阵），并转为密集矩阵  
    matrix = mmread(mtx_file).toarray()  
      
    # 读取barcodes文件  
    barcodes = pd.read_csv(barcodes_file, sep='\t', header=None, names=['barcode'])  
      
    # 若AnnData对象为空，则进行初始化  
    if adata.shape == (0, 0):  
        adata = sc.AnnData(X=matrix, obs=barcodes, var=features)  
    else:  
        # 若非空，则合并数据（假设所有mtx的特征顺序和名称一致）  
        adata_tmp = sc.AnnData(X=matrix, obs=barcodes, var=features)  
        adata = adata.concatenate(*[adata, adata_tmp], axis=0)  
  
# 确保观测值（obs）的索引连续  
adata.obs_names = np.arange(adata.n_obs).astype(str)  
  
# 读取metadata文件  
metadata_df = pd.read_csv(metadata_file)  
  
# 假设metadata_df和adata.obs有一个共同的列'sample_id'用于数据合并  
# 如果metadata_df中的样本ID列名不是'sample_id'，请做相应修改  
# 同时，我们也需要将adata.obs中的样本ID列名设为'sample_id'（如果原本不是）  
adata.obs['sample_id'] = adata.obs_names  # 这里我们简单地将obs_names作为sample_id  
  
# 合并数据（使用'outer'合并确保所有样本都被保留）  
merged_df = pd.merge(adata.obs, metadata_df, on='sample_id', how='outer')  
  
# 根据合并后的DataFrame更新adata的obs属性  
adata.obs = merged_df  
  
# 如果需要，可以根据某些标准对基因进行过滤  
# 例如：sc.pp.filter_genes(adata, min_counts=1)  
  
# 保存整合后的AnnData对象为.h5ad格式文件  
adata.write('merged_data.h5ad')  
  
print("Merged data has been saved to merged_data.h5ad")