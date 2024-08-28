# the code is to prepare files for inferCNV

# 1. load scATAC data
# 2. extract gene activity table
# 3. extract cell cluster annotation table
# 4. identify the malignant cells cluster number

# 1. load scATAC data
x1=readRDS(".//Primary_tumor/scATAC/1_measurement/peak_score_raw_matrix/hg38/filter_cell_on_sample_individual_peak/MGH240_filtered_all_cell.rds")

# 2. extract gene activity table
mm=as.matrix(x1@assays$RNA@counts)
ss(mm,"~/Desktop/MGH240_atac_gene_activity.txt")

# 3. extract cell cluster annotation table
x2=x1@meta.data[,"seurat_clusters",drop=F]
x2$seurat_clusters=paste("c",x2$seurat_clusters,sep = "")
write.table(x2,"~/Desktop/MGH240_atac_anno.txt",quote = F,sep = "\t", col.names = F)

# 4. identify the malignant cells cluster number
all_m=subset(
  x = all_filter,
  subset = seurat_clusters %in% c("1","0","10","13","9")
)

