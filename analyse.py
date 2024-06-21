import scanpy as sc

# 设置文件路径
data_dir = "/path/to/space-ranger/output/outs"

# 读取10x Genomics数据
adata = sc.read_10x_mtx(data_dir + "/filtered_feature_bc_matrix/", var_names='gene_symbols', cache=True)

# 查看数据概况
print(adata)


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

