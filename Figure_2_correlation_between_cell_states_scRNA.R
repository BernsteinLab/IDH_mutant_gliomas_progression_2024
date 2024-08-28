# the code is to quantify the similarity between different cell types/states 
# for normal brain cells and glioma malignant cells based on scRNA-seq data


# 1. identify the dynamic genes
# 2. load tumor data
# 3. random sample 30 cells from each cell state to make pseudo-bulk and repeat 200 times
# 4. calculate the correlation based on pseudo-bulk data [figure 2b]
# 5. load normal brain data
# 6. random sample 30 cells from each cell state to make pseudo-bulk and repeat 200 times
# 7. calculate the correlation based on pseudo-bulk data
# 8. quantify correlation within/across cell states [figure 2c]
# 9. calculate correlation/pvalue between tumor vs normal [figure 2d]

# 1. identify the dynamic genes
gene1=read.table("4_major_program_n.txt") # dynamic genes across all the tumor cells
gene2=read.table("cell_type_gene_specific.txt") # dynamic genes across all the normal cells
gene=unique(c(gene1[,1],gene2[,1]));gene # combine them

# 2. load tumor data
x1=readRDS("MGH_smart_1_m.rds")
q3=x1@meta.data;table(q3$predicted.id)
x33=rbind(q3[which(q3$predicted.id=="ASC"),],head(q3[which(q3$predicted.id=="GluN"),],243),
          head(q3[which(q3$predicted.id=="GPC"),],243),head(q3[which(q3$predicted.id=="OPC"),],243)) # select same amount of cells for each cell state
x3_dd=x1@assays$RNA@data
x3_d=x3_dd[intersect(rownames(x3_dd),gene),] # generate cell x gene matrix for correlation calculation

# 3. random sample 30 cells from each cell state to make pseudo-bulk and repeat 200 times
b=list()
for (name in c("ASC","GPC","OPC","GluN")){
  for (n in c(1:200)) {
    a=list(sample(rownames(x33[which(x33$predicted.id==name),]),30))
    names(a)=paste(name,n,sep = "");a
    b=c(b,a)
  }
};b
c=data.frame(do.call(cbind,lapply(b, function(i)rowMeans(as.matrix(x3_d[,i])))));head(c)

# 4. calculate the correlation based on pseudo-bulk data
corr=cor(c,m="s") # figure 2b

# 5. load normal brain data
x1=readRDS("fetal_brain_scRNA_CP.rds")
q3=x1@meta.data;table(q3$cell_type)
x33=rbind(q3[which(q3$cell_type=="Astro"),],head(q3[which(q3$cell_type=="EN"),],243),
          head(q3[which(q3$cell_type=="GPC"),],243), q3[which(q3$cell_type=="OPC"),])  # select same amount of cells for each cell state
x3_dd=x1@assays$RNA@data
x3_d=x3_dd[intersect(rownames(x3_dd),gene),] # generate cell x gene matrix for correlation calculation

# 6. random sample 30 cells from each cell state to make pseudo-bulk and repeat 200 times
b=list()
for (name in c("Astro","GPC","OPC","EN")){
  for (n in c(1:200)) {
    a=list(sample(rownames(x33[which(x33$cell_type==name),]),30))
    names(a)=paste(name,n,sep = "");a
    b=c(b,a)
  }
};b
c1=data.frame(do.call(cbind,lapply(b, function(i)rowMeans(as.matrix(x3_d[,i])))));head(c1)

# 7. calculate the correlation based on pseudo-bulk data
corr1=cor(c1,m="s") # figure 2b

# 8.quantify correlation within/across cell states -- figure 2c: show plasticity in tumor 
self_tumor_rna=c(mean(unlist(corr[0:200,0:200])),mean(unlist(corr[200:400,200:400])),
                 mean(unlist(corr[400:600,400:600])),mean(unlist(corr[600:800,600:800])))
cross_tumor_rna=c(mean(unlist(corr[0:200,200:800])),mean(unlist(corr[200:400,400:800])),
                  mean(unlist(corr[400:600,600:800])))
self_nor_rna=c(unlist(corr1[0:200,0:200]),unlist(corr1[200:400,200:400]),unlist(corr1[400:600,400:600]),unlist(corr1[600:800,600:800]))
cross_nor_rna=c(unlist(corr1[0:200,200:800]),unlist(corr1[200:400,400:800]),unlist(corr1[400:600,600:800]))
boxplot(self_nor_rna,cross_nor_rna,self_tumor_rna,cross_tumor_rna,outline = F,ylim=c(0.35,1)) 
t.test(self_tumor_rna,cross_tumor_rna,alternative = "greater")
t.test(self_nor_rna,cross_nor_rna,alternative = "greater")

# 9.calculate correlation/pvalue between tumor vs normal
d=mapa(c,c1) # merge tumor and normal pseudo-bulk
x5=cor(d,m="s") 
boxplot(unlist(x5[401:600,1201:1400]),unlist(x5[201:400,1001:1200]),
        unlist(x5[601:800,1401:1600]),outline=F,ylim=c(0.35,0.75),names = c("OPC","GPC","NPC"),main="scRNA")  # figure 2d
t.test(unlist(x5[401:600,1201:1400]),unlist(x5[201:400,1001:1200]),alternative = "greater")
t.test(unlist(x5[401:600,1201:1400]),unlist(x5[601:800,1401:1600]),alternative = "greater")
