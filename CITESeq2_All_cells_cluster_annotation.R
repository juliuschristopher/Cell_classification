## ##CITE-Seq-2 all cells cluster identification#####
#Julius Christopher Baeck

####Set up####
##Set working firectory##
etwd("/Volumes/GoogleDrive/Shared drives/Okkengroup/Experiments/Julius/Experiments/CITE-Sequencing/CITE-Seq (2)/Overall_analysis/CITE-Seq2_all_cells/Cell_classification")

##Packages##
library(clustifyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ArrayExpress)
library(GEOquery)
library(org.Mm.eg.db)
library(viridis)
library(stringr)
library(clustifyrdatahub)
library(devtools)
library(scTyper)
library(rjags)
library(infercnv)
library(Polychrome)


##Colour paneles##
col_con1 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
col_con2 <- createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))
col_con2 <-as.character(col_con2)

## ##Loading data test#####
All_cells <- LoadH5Seurat("CITE-Seq2_all_cells.h5seurat")
All_cells$seurat_clusters <- All_cells$wsnn_res.0.7
Idents(All_cells) <- All_cells$seurat_clusters

##Set levels##
my_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
               "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
               "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",
               "31", "32", "33", "34", "35")
Idents(All_cells) <- factor(Idents(All_cells), levels = my_levels)


UMAP1 <- DimPlot(All_cells, reduction = "wnn.umap", cols = col_con2, pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("All Clusters CITE-Seq-2") +
  theme(plot.title = element_text(size=16, face = "bold"))
print(UMAP1)
ggsave("Original_cluster_1.pdf", width = 30, height = 20, units = "cm")


####Using canonical markers####
##B cells##
DefaultAssay(All_cells) <- "ADT"

Cd19 <- FeaturePlot(All_cells, features = "Cd19", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD19 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Cd19)
ggsave("Cd19_ab.pdf", width = 30, height = 20, units = "cm")


DefaultAssay(All_cells) <- "RNA"
Cd19_rna <- FeaturePlot(All_cells, features = "Cd19", reduction = "wnn.umap", cols = mako(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD19 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Cd19_rna)
ggsave("Cd19_rna.pdf", width = 30, height = 20, units = "cm")

DefaultAssay(All_cells) <- "ADT"
Vln_Cd19 <- VlnPlot(All_cells, features = "Cd19", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD19 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_Cd19)
ggsave("Vln_Cd19.pdf", width = 30, height = 20, units = "cm")
B_cell_clus_adt <- c("0", "1", "2", "10", "15", "18", "19", "20", "21", "23", "27", "28", "30", "31", "35")

DefaultAssay(All_cells) <- "RNA"
Vln_Cd19_rna <- VlnPlot(All_cells, features = "Cd19", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD19 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_Cd19_rna)
ggsave("Vln_Cd19_rna.pdf", width = 30, height = 20, units = "cm")
B_cell_clus_rna <- c("0", "1", "2", "10", "15", "18", "19", "20", "21", "23", "27", "28", "30", "35")



