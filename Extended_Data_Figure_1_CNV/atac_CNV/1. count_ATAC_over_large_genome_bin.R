# the code is to calculate the scATAC-seq signal over the large genome bins 

# 1. prepare genome bins
# 2. load scATAC-seq data
# 3. create object from fragments and granges
# 4. extract malignant and normal cells

library(SingleCellExperiment)
library(GenomicRanges)
library(Signac)
library(future)
library(Seurat)
plan("future::multisession")
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM


# 1. prepare genome bins
bed="~/Desktop/Gliomas//Revision/240221_CNV_analysis/chr_large.bed" # these are genome bins that have blacklist regions removed
peak=read.table(bed)[,1:3]
colnames(peak)=c("chr","start","end");

# 2. load scATAC-seq data
h5="./Primary_tumor/scATAC/0_QC/fragment_files/mgh240/filtered_peak_bc_matrix.h5"
me="./Primary_tumor/scATAC/0_QC/fragment_files/mgh240/singlecell.csv"
fragpath="./Primary_tumor/scATAC/0_QC/fragment_files/mgh240/fragments.tsv.gz"

# 3. create object from fragments and granges
fragcounts <- CountFragments(fragments = fragpath)
atac.cells <- fragcounts[fragcounts$frequency_count > 2000, "CB"]
atac.frags <- CreateFragmentObject(path = fragpath, cells = atac.cells)
tumor_counts <- FeatureMatrix(fragments = atac.frags,features = makeGRangesFromDataFrame(peak),cells = atac.cells) # 
atac.assay <- CreateChromatinAssay(counts = tumor_counts,fragments = atac.frags)
tumor_atac<- CreateSeuratObject(counts = atac.assay, assay = "ATAC")
tumor_atac<- RunTFIDF(tumor_atac)

# 4. extract malignant and normal cells
all_filter=readRDS("Primary_tumor/scATAC/1_measurement/peak_score_raw_matrix/hg38/filter_cell_on_sample_individual_peak/MGH240_filtered_all_cell.rds")
x1=rr("./Primary_tumor/scATAC/0_QC/MGH240_QC/MGH240_filtered_malignant_cell.txt")
filter_cell=rownames(all_filter@meta.data);filter_cell
all=tumor_atac@assays$ATAC@data
all_m=all[,intersect(colnames(all),x1$x)];dim(all_m)
all_n=all[,intersect(colnames(all),filter_cell[!filter_cell%in%x1$x])];dim(all_n)
ss(all_m,"mgh240_tumor_nt.txt")
ss(all_n,"mgh240_normal_nt.txt")