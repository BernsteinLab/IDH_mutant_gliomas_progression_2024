# the code is to cal the CNV number of each single cell

# 1. load scATAC-seq signal over the large genome bins
# 2. cap the data between the quantile 0.2-0.8
# 3. calculate cnv for tumor
# 4. plot and save [Extended Data Figure 1d]
# 5. calculate cnv for normal
# 6. plot and save [Extended Data Figure 1d]

library(pheatmap)
# 1. load scATAC-seq signal over the large genome bins 
pre="MGH240"
t_mat=paste(pre,"_tumor_nt.txt",sep = "")
n_mat=paste(pre,"_normal_nt.txt",sep = "")
tumor1=read.table(t_mat,h=T,row.names = 1);tumor1[1:5,1:5];dim(tumor1)
normal1=read.table(n_mat,h=T,row.names = 1);normal1[1:5,1:5];dim(normal1)

# 2. cap the data between the quantile 0.2-0.8
res=function(x){
  up=0.8;down=0.2
  l=quantile(x,down)
  h=quantile(x,up)
  x[x<l]=l
  x[x>h]=h
  return(x)
}
tumor=apply(tumor1, 1, res);normal=apply(normal1, 1, res)

# 3. calculate cnv for tumor
upc=0.8;downc=0.2 # set the cutoff
total=data.frame(row.names =  rownames(tumor));n=dim(tumor)[2]
for (i in c(1:n)){
  cnv.blmax=quantile(normal[,i],upc) # use normal quantile 0.8 as max cutoff
  cnv.blmin=quantile(normal[,i],downc) # use normal quantile 0.2 as min cutoff
  cnv.tum.bl <- rep(0, length(tumor[,i]))
  m.max <- tumor[,i] > cnv.blmax # positions of the cells with number more than normal max
  m.min <- tumor[,i] < cnv.blmin # positions of the cells with number less than normal min
  cnv.tum.bl[m.max] <- (tumor[,i] - cnv.blmax)[m.max] # only give a value if tumor value is more than normal quantile 0.8
  cnv.tum.bl[m.min] <- (tumor[,i] - cnv.blmin)[m.min] # only give a value if tumor value is less than normal quantile 0.2
  a=data.frame(cnv.tum.bl)
  total=data.frame(total,a)
}
colnames(total)=colnames(tumor);dim(total);total[1:5,1:5]

# 4. plot and save
tumor_pdf=paste(pre,"_tumor_cnv.pdf",sep="")
scale_li=unique(c(seq(-0.6,0.6,length=100)))
col=c(colorRampPalette(colors = c("blue","white"))(40),
      colorRampPalette(colors = c("white"))(20),
      colorRampPalette(colors = c("white","red"))(40))
p=pheatmap(total,cluster_cols = F,color = col,treeheight_row=0,
           show_rownames=F,show_colnames = T,
           breaks = scale_li,clustering_method= "average",fontsize_col = 4,
           clustering_distance_rows = "correlation",border_color = "white");p
pdf(tumor_pdf);p;dev.off()

tumor_txt=paste(pre,"_tumor_cnv.txt",sep = "")
a=total[p$tree_row$order,]
write.table("gene",tumor_txt,quote=F,sep="\t",col.names=F,row.names=F,eol="\t")
write.table(a,tumor_txt,sep="\t",quote = F,append = T)

# 5. calculate cnv for normal
total1=data.frame(row.names =  rownames(normal));n=dim(normal)[2]
for (i in c(1:n)){
  cnv.blmax=quantile(normal[,i],upc) # use normal quantile 0.8 as max cutoff
  cnv.blmin=quantile(normal[,i],downc) # use normal quantile 0.2 as min cutoff
  cnv.tum.bl <- rep(0, length(normal[,i]))
  m.max <- normal[,i] > cnv.blmax 
  m.min <- normal[,i] < cnv.blmin 
  cnv.tum.bl[m.max] <- (normal[,i] - cnv.blmax)[m.max]
  cnv.tum.bl[m.min] <- (normal[,i] - cnv.blmin)[m.min]
  a=data.frame(cnv.tum.bl)
  total1=data.frame(total1,a)
}
colnames(total1)=colnames(normal);dim(total1)

# 6. plot and save
normal_pdf=paste(pre,"_normal_cnv.pdf",sep = "")
p1=pheatmap(total1,cluster_cols = F,color = col,
            treeheight_row=0,show_rownames=F,show_colnames = F,
            breaks = scale_li,
            clustering_method="average",
            border_color="NA");p1
pdf(normal_pdf);p1;dev.off()

normal_txt=paste(pre,"_normal_cnv.txt",sep = "")
b=total1[p1$tree_row$order,]
write.table("gene",normal_txt,quote=F,sep="\t",col.names=F,row.names=F,eol="\t")
write.table(b,normal_txt,sep="\t",quote = F,append = T)