# input gene matrix and annotation
q2=readRDS("./Primary_tumor/scRNA/0_QC/raw_data/malignant_cells/MGH_smart_1_m.rds")
q3=q2@meta.data
col=unique(q3$predicted.id)
col_list=list(predicted.id=brewer.pal(length(col),"Set3")) # pick the # of colors that matches to # of categories
names(col_list$predicted.id)=unique(col);col_list # give the name to the color
q4=q2@assays$RNA@data

# sample 200 cells from each categories
t5=sample(rownames(q3[which(q3$predicted.id=="OPC"),]),size = 200)
t6=sample(rownames(q3[which(q3$predicted.id=="ASC"),]),size = 200)
t7=sample(rownames(q3[which(q3$predicted.id=="GPC"),]),size = 200)
t8=sample(rownames(q3[which(q3$predicted.id=="GluN" | q3$predicted.id=="EN" | q3$predicted.id=="nIPC"),]),200)
cell1=c(t6,t7,t5,t8);length(cell1)

# in gene list
q5=read.table("./Annotation/gene_list/NMF/4_major_program_n.txt")

p=pheatmap(
  mat               = q4[q5[,1],cell1],
  #show_rownames     = FALSE, 
  show_colnames    = FALSE, 
  annotation_col    = q3[,"predicted.id",drop=F],
  annotation_colors = col_list,
  cluster_rows      = F, # boolean values determining if rows should be clustered or hclust object,
  cluster_cols      = F, # boolean values determining if rows should be clustered or hclust object,
  scale            = "row" # "row", "column" and "none"
);p

