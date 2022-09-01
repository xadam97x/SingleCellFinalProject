# SingleCellFinalProject

The final project for the Modeling of Complex Biological Systems 2022 course.
The project includes single-cell RNA-seq analysis of gastric cancer samples - cell clustering, pseudo-time trajectory inference, copy number variation analysis, and cell-cell communication analysis. 
Precise description in [raport file](https://github.com/xadam97x/SingleCellFinalProject/blob/main/raport%20.pdf).


# Required packages
Seurat_4.1.1  
monocle3_1.2.9    
AnnoProbe_0.1.6  
ggalluvial_0.12.3   
cowplot_1.1.1  
infercnv_1.13.0   
CellChat_1.4.0  
patchwork_1.1.2  
dittoSeq_1.9.3  
ggplot2_3.3.6  
dplyr_1.0.9  

# Running the analysis

At first run seurat_analysis.R script for cell clustering. This step must be done before any further analysis. After that cnv.R (copy number variation analysis), Monocle_trajectory.R (pseudo-time trajectory analysis), and cellchat_analysis.R (signaling pathways analysis) can be used.
