# The code is to project malignant cells onto the reference made of integrated fetal and adult normal brain cells

# 1. sample cells from different data sets
# 2. integrate cells from different data sets
# 3. project tumor on integrated datasets

# 1. sample cells from different data sets

# sample cells from fetal brain dataset1
w1=readRDS("./fetal_brain_scRNA_1.rds");table(w1$V2)
c=c()
for (n in c("CGE_IN","Cyc_Prog","Early_RG","GluN1","GluN2","GluN3","GluN4",
            "GluN5","Late_RG","MGE_IN","MG","mGPC","nIPC","OPC_Oligo","tRG")){
  w2=w1@meta.data[which(w1$V2==n),]
  if(dim(w2)[1]>100){b=rownames(head(w2[order(w2$nFeature_RNA,decreasing = T),],n=100))}
  else{b=rownames(w2)}
  c=c(c,b)
}
length(c);table(w1@meta.data[c,"V2"])
w3=w1[,c]

# sample cells from fetal brain dataset2
e1=readRDS("./fetal_brain_scRNA_2.rds");table(e1$cell_type)
c=c()
for (n in c("tRG","Astro","GPC","INP","RG","OPC","EN","IN","ENP") ){
  e2=e1@meta.data[which(e1$cell_type==n),]
  if(dim(e2)[1]>100){b=rownames(head(e2[order(e2$nFeature_RNA,decreasing = T),],n=100))}
  else{b=rownames(e2)}
  c=c(c,b)
}
length(c);table(e1@meta.data[c,"cell_type"])
e3=e1[,c];e3

# sample cells from adult brain dataset
a1=readRDS("./adult_brain_scRNA.rds");table(a1$Cell.Type)
c=c()
for (n in c("ODC","ASC","EX","OPC","MG","INH") ){
  a2=a1@meta.data[which(a1$Cell.Type==n),]
  if(dim(a2)[1]>100){b=rownames(head(a2[order(a2$nFeature_RNA,decreasing = T),],n=100))}
  else{b=rownames(a2)}
  c=c(c,b)
}
length(c);table(a1@meta.data[c,"Cell.Type"])
a3=a1[,c];a3

# 2. integrate cells from different data sets
list=list(w3,e3,a3)
features <- SelectIntegrationFeatures(object.list = list)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
rna.combined <- IntegrateData(anchorset = anchors)
DefaultAssay(rna.combined) <- "integrated"

# recluster it
rna.combined.n=subset(x = rna.combined, subset = Cell.Type != "MG")
rna.combined.n=subset(x = rna.combined.n, subset = Cell.Type != "tRG")
rna.combined.n =FindVariableFeatures(rna.combined.n)
rna.combined.n <- ScaleData(rna.combined.n, verbose = FALSE)
rna.combined.n <- RunPCA(rna.combined.n, npcs = 20, verbose = FALSE)
rna.combined.n <- RunUMAP(rna.combined.n, reduction = "pca", dims = 1:20,return.model=TRUE)
rna.combined.n<- FindNeighbors(rna.combined.n, dims = 1:15)
rna.combined.n<- FindClusters(rna.combined.n, resolution = 0.5)
DimPlot(rna.combined.n, reduction = "umap",group.by = "seurat_clusters",label = TRUE)

# add metadata to it
me1=w1@meta.data;colnames(me1)[4]="Cell.Type";hed(me1);me1$stage="fetal"
me2=a1@meta.data[,c(1,2,3,7)];me2$stage="adult";hed(me2)
me3=e1@meta.data;colnames(me3)[4]="Cell.Type";me3$stage="fetal";hed(me3)
me=rbind(me1,me2,me3);hed(me)
old_m=rna.combined.n@meta.data
rna.combined.n=AddMetaData(rna.combined.n,metadata = me)

# simplify the cell type 
b=rna.combined.n$Cell.Type
b=gsub(pattern = "OPC_Oligo",replacement = "OPC",b)
b=gsub(pattern = "GluN\\d",replacement = "GluN",b)
b=gsub(pattern = "MGE_IN",replacement = "IN",b)
b=gsub(pattern = "CGE_IN",replacement = "IN",b)
b=gsub(pattern = "PER.END",replacement = "Peric",b)
b=gsub(pattern = "EC",replacement = "Peric",b)
b=gsub(pattern = "Early_RG",replacement = "ASC",b)
b=gsub(pattern = "Late_RG",replacement = "ASC",b)
#b=gsub(pattern = "tRG",replacement = "RG",b)
b=gsub(pattern = "Astro",replacement = "ASC",b)
b=gsub(pattern = "mGPC",replacement = "GPC",b)
b=gsub(pattern = "INP",replacement = "Cyc_Prog",b)
b=gsub(pattern = "ENP",replacement = "Cyc_Prog",b)
b=gsub(pattern = "EN",replacement = "GluN",b)
b=gsub(pattern = "EX",replacement = "GluN",b)
b=gsub(pattern = "INH",replacement = "IN",b)
rna.combined.n$new_type=b

# 3. project tumor on integrated datasets
t1=readRDS("./tumor_malignant_cell.rds")
anchors <- FindTransferAnchors(reference = rna.combined.n, query = t1, dims = 1:20,reference.reduction="pca",k.filter = NA)
predictions <- TransferData(anchorset = anchors, refdata = rna.combined.n$new_type,dims = 1:20) # included in the mapQuery
t1 <- AddMetaData(t1, metadata = predictions)
t1 <- MapQuery(anchorset = anchors, reference = rna.combined.n, query = t1,
               refdata = list(celltype = "new_type"), reference.reduction = "pca", reduction.model = "umap")
