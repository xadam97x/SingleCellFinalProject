library(monocle3)
library(Seurat)
epi<-readRDS("Epithelial_seurat")

DefaultAssay(epi)<-"integrated"

gene_annotation <- as.data.frame(rownames(epi@reductions[["pca"]]@feature.loadings), row.names = rownames(epi@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
cells_info <- epi@meta.data
New_matrix <- epi@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(epi@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cells_info,
                                     gene_metadata = gene_annotation)


### Construct and assign the made up partition

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


### Assign the cluster info

list_cluster <- epi@meta.data[["integrated_snn_res.0.2"]]
names(list_cluster) <- epi@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_from_seurat@int_colData@listData$reducedDims$UMAP<-epi@reductions[["umap"]]@cell.embeddings


cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)

get_earliest_principal_node <- function(cds){
  cell_ids <- which(colData(cds)[, "malignant_score"] < -0.2)
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
}

cds_from_seurat<- order_cells(cds_from_seurat,reduction_method = "UMAP",root_pr_nodes=get_earliest_principal_node(cds_from_seurat))
plot_cells(cds_from_seurat,color_cells_by = "pseudotime",label_groups_by_cluster = TRUE, cell_size = 1)

deg<-graph_test(cds_from_seurat,neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status=="OK")->deg
head(deg,10)
saveRDS(cds_from_seurat,"epi_monocle")





