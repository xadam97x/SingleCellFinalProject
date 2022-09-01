library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(CellChat)
library(ggalluvial)


combined<-readRDS("combined_pc20_random_res08")
levels(combined@active.ident)
new.idents<-c("TCD8",'TCD4','B','Neutro','Epi','MaMo','NK','Treg','Plasma','Endo','Fib','Proli','MAST')
old<-c("T_CD8","T_CD4",'B','Neutro','Epithelial','MaMo','NK',"Treg","Plasma","Endothelial","Fibroblast","Proliferative","MAST")
combined@active.ident<- plyr::mapvalues(combined@active.ident,from = old, to= new.idents)
cellchat <- createCellChat(combined@assays$RNA@data)
dir.create("./cellchat_output")
setwd("./cellchat_output")
identity = data.frame(group =combined@active.ident, row.names = names(combined@active.ident)) #  cell labels
unique(identity$group)

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels") ##add lables
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) 
DB<-CellChatDB.human
DB.use<-subsetDB(DB,search = "Secreted Signaling")
cellchat@DB<-DB.use
cellchat<-subsetData(cellchat)
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-projectData(cellchat,PPI.human)
cellchat<-computeCommunProb(cellchat,raw.use=T)
cellchat<-filterCommunication(cellchat,min.cells = 10)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))
pathways.show<-cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = seq(1,4)
cellchat<-netAnalysis_computeCentrality(cellchat,slot.name = "netP")
pathways.show

for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show[i], layout = "circle", vertex.weight = groupSize,pt.title=20,vertex.label.cex = 1.7)
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 2, units = 'in', dpi = 300)
  # Visualize signaling roles of cell groups
  pdf(file = paste0(pathways.show[i], "_signalRole.pdf"))
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show[i], width = 8, height = 2.5, font.size = 10)
  dev.off()
  pdf(file =paste0(pathways.show[i], "_chord_adjust.pdf"))
  netVisual_chord_cell(cellchat, signaling = pathways.show[i], title.name = paste0(pathways.show[i], " signaling network")) 
  dev.off()
  pdf(file = paste0(pathways.show[i], "_heatmap.pdf"))
  netVisual_heatmap(cellchat, signaling = pathways.show[i], color.heatmap = "Reds") 
  dev.off()
  pdf(file = paste0(pathways.show[i], "_gene.pdf"))
  print(plotGeneExpression(cellchat, signaling =pathways.show[i]))
  dev.off()
}


ident <-c(1:13)
for (i in 1:length(ident)) {
  pdf(file =paste0(levels(cellchat@idents)[i], "_LR.pdf"))
  print(netVisual_bubble(cellchat, sources.use = i, targets.use = ident[-i], remove.isolate = FALSE))
  dev.off()
  pdf(file =paste0(levels(cellchat@idents)[i], "_netVisual_chord.pdf"))
  netVisual_chord_gene(cellchat, sources.use = i, targets.use = ident[-i], lab.cex = 0.5,legend.pos.y = 30) 
  dev.off()
}




pdf(file = "Communication_strength_weight.pdf")
par(mfrow = c(1,2), xpd=TRUE) 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions") 
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


pdf(file = "Communication_weight.pdf")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE) 
for (i in 1:nrow(mat)) { 
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)) 
  mat2[i, ] <- mat[i, ] 
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]) 
}
dev.off()

nPatterns = 5
pdf(file = paste0("CommunicationPatterns_sender_heatmap.pdf"))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()


pdf(file = "patternAnalysis_sender_river.pdf", width = 7, height = 4)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()


gg <- netAnalysis_dot(cellchat, pattern = "outgoing")
ggsave(filename="patternAnalysis_sender_dot.pdf", plot=gg, width = 5.5, height = 4, units = 'in', dpi = 300)
pdf(file = paste0("CommunicationPatterns_receiver_heatmap.pdf"))
dev.off()

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
pdf(file = "patternAnalysis_receiver_river.pdf", width = 7, height = 4)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
gg <- netAnalysis_dot(cellchat, pattern = "incoming")
ggsave(filename="patternAnalysis_receiver_dot.pdf", plot=gg, width = 5.5, height = 4, units = 'in', dpi = 300)
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "functional", thresh = 0.25)    
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional", k = 4)
gg <- netVisual_embedding(cellchat, type = "functional", pathway.remove.show = F)
cowplot::save_plot("2Dmanifold_FunctionalSimilarity_signalingPathways.pdf", gg, base_height = 3, base_width = 4)
pdf(file = "2Dmanifold_FunctionalSimilarity_signalingPathways_zoomIn.pdf", width = 2, height = 2.5*3)
netVisual_embeddingZoomIn(cellchat, type = "functional")
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "structural", thresh = 0.25)
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
gg <- netVisual_embedding(cellchat, type = "structural")
cowplot::save_plot("2Dmanifold_StructureSimilarity_signalingPathways.pdf", gg, base_height = 3, base_width = 4)
pdf(file = "2Dmanifold_StructureSimilarity_signalingPathways_zoomIn.pdf", width = 2, height = 2.5*3)
netVisual_embeddingZoomIn(cellchat, type = "structural")
dev.off()

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing") 
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming") 
pdf(file = "heatmap_outgoing_incoming.pdf", width =10, height = 5)
print(ht1 + ht2)
dev.off()
gc()


saveRDS(cellchat,"cellchat_final")
