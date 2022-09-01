library(Seurat)
library(dplyr)
library(patchwork)
library(dittoSeq)
set.seed(2022)


pt1<-Read10X(data.dir = "sc_data/PT1")
pt2<-Read10X(data.dir = "sc_data/PT2")
pt3<-Read10X(data.dir = "sc_data/PT3")
li1<-Read10X(data.dir = "sc_data/Li1")
li2<-Read10X(data.dir = "sc_data/Li2")
ln1<-Read10X(data.dir = "sc_data/LN1")
ln2<-Read10X(data.dir = "sc_data/LN2")
nt1<-Read10X(data.dir = "sc_data/NT1")
p1<-Read10X(data.dir = "sc_data/P1")
o1<-Read10X(data.dir = "sc_data/O1")

pt1<-CreateSeuratObject(pt1,project = "PT1")
pt2<-CreateSeuratObject(pt2,project = "PT2")
pt3<-CreateSeuratObject(pt3,project = "PT3")
li1<-CreateSeuratObject(li1,project = "Li1")
li2<-CreateSeuratObject(li2,project = "Li2")
ln1<-CreateSeuratObject(ln1,project = "LN1")
ln2<-CreateSeuratObject(ln2,project = "LN2")
nt1<-CreateSeuratObject(nt1,project = "NT1")
p1<-CreateSeuratObject(p1,project = "P1")
o1<-CreateSeuratObject(o1,project = "O1")

quality_mesures<- function(a) {
  a[['percent.mt']]<-PercentageFeatureSet(a,pattern = '^MT-')
  return(a)
}

samples<-list(pt1,pt2,pt3,li1,li2,ln1,ln2,nt1,p1,o1)
samples<-lapply(samples,FUN = quality_mesures)
rm(pt1,pt2,pt3,li1,li2,ln1,ln2,nt1,p1,o1)
gc()

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
for (i in samples){
  print(VlnPlot(i,features = feats,pt.size = 0.1, ncol = 3))
}
samples<-lapply(samples, FUN = function(x){
  x<-subset(x,subset = percent.mt<20 & nFeature_RNA>200&nFeature_RNA<5000)
  return(x)
})
samples<-lapply(samples,FUN=function(x){
  x<-subset(x,downsample=1295)
  return(x)
})


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

samples <-lapply(samples, FUN = function(x){ 
  x<-NormalizeData(x)
  x<-CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  return(x)})

samples<-lapply(samples,FUN = function(x){
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)})
features<-SelectIntegrationFeatures(object.list = samples)
anchors<-FindIntegrationAnchors(object.list = samples,anchor.features = features)
combined<-IntegrateData(anchorset = anchors)

rm(samples,anchors)
gc()

DefaultAssay(combined)<-"integrated"
combined<-ScaleData(combined)
combined<-RunPCA(combined)
combined<-RunTSNE(combined,dims = 1:20)
combined<-FindNeighbors(combined,dims = 1:20)
combined<-FindClusters(combined, resolution = 0.8)
DimPlot(combined, reduction = "tsne", label = TRUE)





DefaultAssay(combined)<-"RNA"
combined_markers<- FindAllMarkers(object = combined, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25, verbose=FALSE
                                   )
combined_markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)->top30
write.table(top30,"markers_res08.txt",row.names=TRUE,col.names=TRUE)

old.ids<-c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
new.ids<-c("T_CD8","T_CD4","T_CD4","B","T_CD4","Neutro","Epithelial","T_CD4","Neutro","B","MaMo","NK","Treg","Plasma","Endothelial","Fibroblast","MaMo","Plasma","Proliferative","MAST")
combined@meta.data$type<-plyr::mapvalues(x=combined@meta.data$integrated_snn_res.0.8,from = old.ids,to=new.ids)
names(new.ids)<-levels(combined)
combined<-RenameIdents(combined,new.ids)

DimPlot(combined,reduction = "tsne",split.by = "Phase",ncol = 2)

saveRDS(combined,"combined_pc20_random_res08")



# Re Epithelial cells
epi<-subset(combined,subset = type=="Epithelial")
DefaultAssay(epi)<-"RNA"
epi<-FindVariableFeatures(epi)
DefaultAssay(epi)<-"integrated"
epi <- ScaleData(epi, verbose = FALSE)
epi<- RunPCA(epi, verbose = FALSE)
ElbowPlot(epi)
epi <- RunUMAP(epi, reduction = "pca", dims = 1:20,verbose = FALSE)
epi <- FindNeighbors(epi, reduction = "pca", dims = 1:20,verbose = FALSE)
epi <- FindClusters(epi, resolution = 0.2,verbose=FALSE)
DimPlot(epi, reduction = "umap",label = T)


malignant_genes<- read.delim("./malignant_genes",sep = " ")
malignant_genes<-malignant_genes[,2]

non_malingnat<-read.delim("./normal_tissue",sep=" ")
non_malingnat<-non_malingnat[,2]


DefaultAssay(epi)<-"RNA"
epi<-AddModuleScore(epi,list(malignant_genes, non_malingnat),name="malignant")
epi@meta.data$malignant_score<-epi@meta.data$malignant1-epi@meta.data$malignant2
RidgePlot(epi,features = "malignant_score")
epi@meta.data$malignancy<-epi@meta.data$malignant_score>-0.02
epi@meta.data$malignancy<-plyr::mapvalues(x=epi@meta.data$malignancy,from = c(T,F),to=c("Malignant","Normal"))
saveRDS(epi,"Epithelial_seurat")
epi<-SetIdent(epi,value="malignancy")


cells.num<-table(combined@active.ident)
cells.num

dittoBarPlot(combined,"type",group.by = "orig.ident",scale = "percent")
dittoBarPlot(epi,"malignancy",group.by = "orig.ident")

