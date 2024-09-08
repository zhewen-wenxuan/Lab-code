import scanpy as sc  
import pandas as pd  
  
# 假设您的matrix文件是一个CSV或TSV文件，第一列是基因名，后续列是表达值  
# 且列名（除了第一列）与barcodes文件中的行对应  
matrix_filename = 'matrix.csv'  # 或者matrix.tsv，根据您的文件类型  
barcodes_filename = 'merged_barcodes.tsv'  # 这是您之前生成的包含标识的barcodes文件  
output_filename = 'output.h5ad'  
  
# 读取matrix  
matrix = pd.read_csv(matrix_filename, index_col=0)  # 假设第一列是基因名  
  
# 读取barcodes  
barcodes = pd.read_csv(barcodes_filename, header=None, names=['barcode'])  
  
# 确保barcodes的顺序与matrix的列顺序一致  
# 如果您的barcodes文件已经按照matrix的列顺序排列，这一步可能不是必需的  
# 但如果顺序不一致，您需要根据实际情况调整  
# 例如，如果matrix的列名存储在某个列表中，您可以使用这个列表来重新排序barcodes  
# 这里我们假设它们已经是一致的  
  
# 创建AnnData对象  
adata = sc.AnnData(matrix, obs=barcodes, obs_names=barcodes['barcode'])  
  
# 添加批次信息（这里我们假设批次信息存储在barcodes的某个部分，例如'_sample1'或'_sample2'）  
# 您需要根据您的实际情况调整这一部分  
# 下面的代码是一个示例，它提取barcode中的批次信息并添加到adata的obsm字典中  
# 如果批次信息不是这样存储的，请相应修改  
adata.obsm['X_batch'] = adata.obs['barcode'].str.extract('_(sample\d+)')[0].astype('category').cat.codes  
  
# 减弱批次效应（使用scanpy的bbknn或mutual nearest neighbors等方法）  
# 这里我们使用scanpy的内置函数regress_out，但请注意，这只是一个简单的线性回归  
# 对于更复杂的批次效应校正，您可能需要使用更高级的方法，如Harmony、Seurat的Mutual Nearest Neighbors等  
# regress_out会直接在adata对象上原地修改表达矩阵  
sc.tl.pca(adata, svd_solver='arpack')  # 首先运行PCA，这是许多后续步骤的基础  
sc.pp.regress_out(adata, ['X_batch'])  # 回归掉批次效应  
  
# 另外，如果您想使用bbknn进行批次效应校正，可以尝试以下代码（但这需要更多的设置和理解）：  
# from sklearn.neighbors import NearestNeighbors  
# sc.tl.ppca(adata)  # 使用PPCA进行降维，这是bbknn的推荐步骤  
# bbknn = sc.tl.bbknn(adata, batch_key='X_batch', n_pcs=50)  # 调整n_pcs以匹配您使用的PCA组件数  
# adata.obsp['connectivities'] = bbknn  # 将bbknn结果存储在AnnData对象中  
# 然后，您可以使用这些连接性进行后续的聚类、降维和可视化，但请注意，bbknn本身并不直接修改表达矩阵  
  
# 保存AnnData对象为.h5ad文件  
adata.write(output_filename)  
  
print(f"AnnData object has been saved to {output_filename}")