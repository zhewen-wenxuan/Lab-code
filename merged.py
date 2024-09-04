import scipy.io
import scipy.sparse as sp
import numpy as np

# 假设你有3个样本：sample1, sample2, sample3
matrix_files = ['sample1/matrix.mtx', 'sample2/matrix.mtx', 'sample3/matrix.mtx']
matrices = []

for file in matrix_files:
    matrix = scipy.io.mmread(file).tocsc()  # 读取并转换为稀疏矩阵格式
    matrices.append(matrix)

# 将稀疏矩阵按列合并
merged_matrix = sp.hstack(matrices)

# 保存合并后的矩阵
scipy.io.mmwrite("merged_matrix.mtx", merged_matrix)

barcodes_files = ['sample1/barcodes.tsv', 'sample2/barcodes.tsv', 'sample3/barcodes.tsv']
merged_barcodes = []

for i, file in enumerate(barcodes_files):
    with open(file, 'r') as f:
        barcodes = f.read().splitlines()
    # 为每个样本的条形码添加样本ID前缀
    barcodes = [f'sample{i+1}_{barcode}' for barcode in barcodes]
    merged_barcodes.extend(barcodes)

# 保存合并后的条形码
with open("merged_barcodes.tsv", 'w') as f:
    f.write('\n'.join(merged_barcodes))

import shutil

# 假设所有样本的features文件是相同的，你可以直接复制一个
shutil.copyfile("sample1/features.tsv", "merged_features.tsv")
