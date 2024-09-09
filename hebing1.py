import scanpy as sc  
import pandas as pd  
import numpy as np  
from scipy.io import mmread  
import os  
import re  
  
# 设定文件夹路径  
data_folder = '/mnt/merged-data'  # 请替换为您的数据文件夹路径  
save_folder = '/mnt/merged-data'  

# 获取文件夹中所有文件的列表  
file_list = os.listdir(data_folder)  
  
# 分离出matrix文件和barcodes文件  
matrix_files = [f for f in file_list if f.endswith('.mtx')]  
barcodes_files = [f for f in file_list if f.endswith('_barcodes_modified.tsv')]  # 假设barcodes文件以_barcodes.tsv结尾  
  
# 确保matrix文件和barcodes文件数量相同  
if len(matrix_files) != len(barcodes_files):  
    raise ValueError("Number of matrix files does not match number of barcodes files.")  
  
# 定义一个函数，用于从文件名中提取排序关键字（假设是数字）  
def extract_sort_key(filename):  
    # 使用正则表达式提取文件名中的数字部分作为排序关键字  
    # 假设数字部分是由一个或多个数字组成，并且前面和后面可能有其他字符  
    match = re.search(r'sample(\d+)', filename) 
    if match:  
        return int(match.group())  
    else:  
        # 如果没有找到数字，返回一个很大的数字以确保这些文件排在最后（或者抛出一个错误）  
        return float('inf')  
    
  
# 对matrix文件和barcodes文件进行排序  
matrix_files_sorted = sorted(matrix_files, key=extract_sort_key)  
barcodes_files_sorted = sorted(barcodes_files, key=extract_sort_key)  
  
# 初始化一个空的列表来存储所有的矩阵和条形码数据  
matrices = []  
barcodes_data = []  
  
# 遍历排序后的文件，加载数据  
for mtx_file, barcodes_file in zip(matrix_files_sorted, barcodes_files_sorted):  
    # 读取mtx文件（稀疏矩阵）  
    matrix = mmread(os.path.join(data_folder, mtx_file))  
      
    # 读取barcodes文件  
    barcodes = pd.read_csv(os.path.join(data_folder, barcodes_file), sep='\t', header=None, names=['barcode'])  
      
    # 存储矩阵和条形码数据  
    matrices.append(matrix)  
    barcodes_data.append(barcodes)  
  
# 合并所有的稀疏矩阵（使用scipy.sparse的vstack方法）  
# 注意：这里我们保持矩阵为稀疏状态  
from scipy.sparse import vstack  
merged_matrix = vstack(matrices)  
  
# 使用pd.concat合并条形码数据（按行拼接）  
merged_barcodes = pd.concat(barcodes_data, ignore_index=True) 
  
# 假设只有一个features文件（或者您已经知道要使用哪个）  
# 这里我们假设features文件以_features.tsv结尾，并且只有一个这样的文件  
features_file = [f for f in file_list if f.endswith('_features.tsv')]  
if len(features_file) != 1:  
    raise ValueError("Expected exactly one features file.")  
features_file = features_file[0]  
  
# 读取features文件  
features = pd.read_csv(os.path.join(data_folder, features_file), sep='\t', header=None, names=['feature'])  

# 保存merged_matrix为稀疏矩阵格式（.mtx）  
matrix_save_path = os.path.join(save_folder, 'merged_matrix.mtx')  
mmwrite(matrix_save_path, merged_matrix)  
  
# 保存merged_barcodes为.tsv格式  
barcodes_save_path = os.path.join(save_folder, 'merged_barcodes.tsv')  
merged_barcodes.to_csv(barcodes_save_path, sep='\t', index=False, header=False)  

# 创建AnnData对象  
adata = sc.AnnData(X=merged_matrix, obs=merged_barcodes, var=features)  
  
# 确保观测值（obs）的索引连续（如果需要的话）  
adata.obs_names = np.arange(adata.n_obs).astype(str)  
  
# 打印确认信息  
print(f"Merged data has {adata.n_obs} observations and {adata.n_vars} variables.")  
print(f"AnnData object has been created successfully.")  
  
# 在此处可以继续进行下一步操作，例如保存adata对象或进行后续分析  
# ...