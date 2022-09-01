library(infercnv)
library(Seurat)
library(AnnoProbe)
set.seed(2022)
#setwd("./test")
combined<-readRDS("combined_pc20_random_res08")
obj<-subset(combined,subset = type=="Epithelial"|type=="B")
rm(combined)
gc()
dfcount<-as.data.frame(obj@assays$RNA@counts)
groupinfo<-data.frame(cellID=colnames(dfcount),cellType=obj@meta.data$type)

geneInfor<-annoGene(rownames(dfcount),"SYMBOL","human")
geneInfor<-geneInfor[with(geneInfor,order(chr,start)),c(1,4:6)]
geneInfor<-geneInfor[!duplicated(geneInfor[,1]),]
dfcount<-dfcount[rownames(dfcount) %in% geneInfor[,1],]
dfcount<-dfcount[match(geneInfor[,1], rownames(dfcount)),]

write.table(dfcount, file = "expFile.txt",sep = '\t',quote = F)
write.table(groupinfo,file='metaFile.txt',sep='\t',quote=F, col.names = F,row.names=F)
write.table(geneInfor,file='geneFile.txt',sep='\t',quote=F, col.names = F,row.names=F)


cnvobj<-CreateInfercnvObject(raw_counts_matrix = "./expFile.txt",annotations_file = "./metaFile.txt",
                             delim = "\t",gene_order_file = "./geneFile.txt",ref_group_names = c("B"))

cnvobj<-infercnv::run(cnvobj,cutoff = 0.1, 
                      out_dir = "./cnvout",cluster_by_groups = F,
                      denoise = TRUE,HMM = TRUE,
                      output_format = "pdf",
                      analysis_mode = "subclusters",
                      )

#change chromosomes order
new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
  new_gene_order = rbind(new_gene_order, cnvobj@gene_order[which(cnvobj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
}
names(new_gene_order) <- c("chr", "start", "stop")
copy_infercnv_obj<-cnvobj
copy_infercnv_obj@gene_order = new_gene_order
copy_infercnv_obj@expr.data = cnvobj@expr.data[rownames(new_gene_order), , drop=FALSE]
plot_cnv(copy_infercnv_obj,
         cluster_by_groups = F,
         cluster_references=TRUE,
         out_dir="./",
         output_filename="infercnv.png")
