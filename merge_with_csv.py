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