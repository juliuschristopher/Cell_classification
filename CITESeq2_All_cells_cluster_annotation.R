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

##CD4+ cells##
DefaultAssay(All_cells) <- "ADT"
Cd4 <- FeaturePlot(All_cells, features = "Cd4", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Cd4)
ggsave("Cd4_ab.pdf", width = 30, height = 20, units = "cm")


DefaultAssay(All_cells) <- "RNA"
Cd4_rna <- FeaturePlot(All_cells, features = "Cd4", reduction = "wnn.umap", cols = mako(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Cd4_rna)
ggsave("Cd4_rna.pdf", width = 30, height = 20, units = "cm")

DefaultAssay(All_cells) <- "ADT"
Vln_Cd4 <- VlnPlot(All_cells, features = "Cd4", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD4 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_Cd4)
ggsave("Vln_Cd4.pdf", width = 30, height = 20, units = "cm")
CD4_cell_clus_adt <- c("4", "5", "6", "7", "17", "26", "27", "29", "33")

DefaultAssay(All_cells) <- "RNA"
Vln_Cd4_rna <- VlnPlot(All_cells, features = "Cd4", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD4 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_Cd4_rna)
ggsave("Vln_Cd4_rna.pdf", width = 30, height = 20, units = "cm")
CD4_cell_clus_rna <- c("4", "5", "6", "7", "17", "27", "33")

##CD8+ cells##
DefaultAssay(All_cells) <- "ADT"
Cd8a <- FeaturePlot(All_cells, features = "Cd8a", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD8 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Cd8a)
ggsave("Cd8a_ab.pdf", width = 30, height = 20, units = "cm")


DefaultAssay(All_cells) <- "RNA"
Cd8a_rna <- FeaturePlot(All_cells, features = "Cd8a", reduction = "wnn.umap", cols = mako(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD8 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Cd8a_rna)
ggsave("Cd8a_rna.pdf", width = 30, height = 20, units = "cm")

DefaultAssay(All_cells) <- "ADT"
Vln_Cd8a <- VlnPlot(All_cells, features = "Cd8a", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD8 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_Cd8a)
ggsave("Vln_Cd8a.pdf", width = 30, height = 20, units = "cm")
CD8_cell_clus_adt <- c("4", "5", "9", "22", "28", "29")

DefaultAssay(All_cells) <- "RNA"
Vln_Cd8a_rna <- VlnPlot(All_cells, features = "Cd8a", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD8 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_Cd8a_rna)
ggsave("Vln_Cd8_rna.pdf", width = 30, height = 20, units = "cm")
CD8_cell_clus_rna <- c("4", "5", "9", "22", "28", "29", "33")

##CD3+ expression to identify "true" T cells##
DefaultAssay(All_cells) <- "RNA"
Cd3d_rna <- FeaturePlot(All_cells, features = "Cd3d", reduction = "wnn.umap", cols = mako(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD3d mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Cd3d_rna)
ggsave("Cd3d_rna.pdf", width = 30, height = 20, units = "cm")

DefaultAssay(All_cells) <- "RNA"
Vln_Cd3d_rna <- VlnPlot(All_cells, features = "Cd3d", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD3d mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_Cd3d_rna)
ggsave("Vln_Cd3d_rna.pdf", width = 30, height = 20, units = "cm")
CD3d_cell_clus_rna <- c("3", "4", "5", "6", "7", "9", "26", "27", "28", "29")

DefaultAssay(All_cells) <- "RNA"
Cd3e_rna <- FeaturePlot(All_cells, features = "Cd3e", reduction = "wnn.umap", cols = mako(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD3e mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Cd3e_rna)
ggsave("Cd3e_rna.pdf", width = 30, height = 20, units = "cm")

DefaultAssay(All_cells) <- "RNA"
Vln_Cd3e_rna <- VlnPlot(All_cells, features = "Cd3e", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD3e mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_Cd3e_rna)
ggsave("Vln_Cd3e_rna.pdf", width = 30, height = 20, units = "cm")
CD3e_cell_clus_rna <- c("3", "4", "5", "6", "7", "9", "26", "27", "28", "29")

##True CD4 and CD8 T cells##
B_cell <- c("0", "1", "2", "10", "15", "18", "19", "20", "21", "23", "30", "35")
CD4_Tcell <- c("4", "6", "7", "26")
CD8_Tcell <- c("3", "9", "29")
CD4_cell <- c("17")
CD8_cell <- c("22")
CD4_CD8_Tcell <- c("5")
CD4_CD8_cell <-c("33")
CD4_Tcell_Bcell <-c("27")
CD4_CD8_Tcell_Bcell <- c("28")
Other_cells <- c("8", "11", "12", "13", "14", "16", "24", "25", "31", "32", "34")

####B cell characterisation####
setwd("/Volumes/GoogleDrive/Shared drives/Okkengroup/Experiments/Julius/Experiments/CITE-Sequencing/CITE-Seq (2)/Overall_analysis/CITE-Seq2_all_cells/Cell_classification/Bcells")

B_cell_gr <- WhichCells(All_cells, idents = B_cell)

B_cell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
        label = TRUE, cells.highlight = B_cell_gr, cols.highlight = "steelblue3"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("B cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(B_cell_gr_plot)
ggsave("B_cell_gr_plot.pdf", width = 30, height = 20, units = "cm")

##Make a UMAP plot with only B cell clusters##
Bcell_clus <- subset(All_cells, idents = B_cell)

my_levels2 <- c("0", "1", "2", "10", "15", "18", "19", "20", "21", "23", "30", "35")
Idents(Bcell_clus) <- factor(Idents(Bcell_clus), levels = my_levels2)

Bcell_clus_plot <- DimPlot(Bcell_clus, reduction = "wnn.umap", cols = col_con2) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("B cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Bcell_clus_plot)
ggsave("Bcell_clus_plot.pdf", width = 30, height = 20, units = "cm")

##Define by common B cell markers##
#ADT UMAP
DefaultAssay(Bcell_clus) <- "ADT"
B220 <- FeaturePlot(Bcell_clus, features = "B220", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("B220 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(B220)
ggsave("B220_ab.pdf", width = 30, height = 20, units = "cm")

CD93 <- FeaturePlot(Bcell_clus, features = "Cd93", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD93 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD93)
ggsave("CD93_ab.pdf", width = 30, height = 20, units = "cm")

CD21 <- FeaturePlot(Bcell_clus, features = "Cd21", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD21 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD21)
ggsave("CD21_ab.pdf", width = 30, height = 20, units = "cm")

CD23 <- FeaturePlot(Bcell_clus, features = "Cd23", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD23 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD23)
ggsave("CD23_ab.pdf", width = 30, height = 20, units = "cm")

IgM <- FeaturePlot(Bcell_clus, features = "Igm", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("IgM cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(IgM)
ggsave("IgM_ab.pdf", width = 30, height = 20, units = "cm")

IgD <- FeaturePlot(Bcell_clus, features = "Igd", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("IgD cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(IgD)
ggsave("IgD_ab.pdf", width = 30, height = 20, units = "cm")

CD38 <- FeaturePlot(Bcell_clus, features = "Cd38", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD38 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD38)
ggsave("CD38_ab.pdf", width = 30, height = 20, units = "cm")

CD95 <- FeaturePlot(Bcell_clus, features = "Cd95", reduction = "wnn.umap", cols = magma(10), pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD95 (FAS) cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD95)
ggsave("CD95_ab.pdf", width = 30, height = 20, units = "cm")

CTLA4 <- FeaturePlot(Bcell_clus, features = "Ctla4", reduction = "wnn.umap", cols = magma(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CTLA4 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CTLA4)
ggsave("CTLA4_ab.pdf", width = 30, height = 20, units = "cm")

PD_1 <- FeaturePlot(Bcell_clus, features = "Pd-1", reduction = "wnn.umap", cols = magma(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("Pd-1 (FAS) cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(PD_1)
ggsave("PD_1_ab.pdf", width = 30, height = 20, units = "cm")

PD_L1 <- FeaturePlot(Bcell_clus, features = "Pd-L1", reduction = "wnn.umap", cols = magma(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("PD-L1 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(PD_L1)
ggsave("PD_L1_ab.pdf", width = 30, height = 20, units = "cm")

PD_L2 <- FeaturePlot(Bcell_clus, features = "Pd-L2", reduction = "wnn.umap", cols = magma(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("PD-L1 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(PD_L2)
ggsave("PD_L2_ab.pdf", width = 30, height = 20, units = "cm")

CD80 <- FeaturePlot(Bcell_clus, features = "Cd80", reduction = "wnn.umap", cols = magma(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD80 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD80)
ggsave("CD80_ab.pdf", width = 30, height = 20, units = "cm")

CD83 <- FeaturePlot(Bcell_clus, features = "Cd83", reduction = "wnn.umap", cols = magma(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD83 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD83)
ggsave("CD83_ab.pdf", width = 30, height = 20, units = "cm")

CD86 <- FeaturePlot(Bcell_clus, features = "Cd86", reduction = "wnn.umap", cols = magma(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD86 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD86)
ggsave("CD86_ab.pdf", width = 30, height = 20, units = "cm")


#RNA UMAP
DefaultAssay(Bcell_clus) <- "RNA"
CD93_rna <- FeaturePlot(Bcell_clus, features = "Cd93", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD93 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD93_rna)
ggsave("CD93_rna.pdf", width = 30, height = 20, units = "cm")

CD21_rna <- FeaturePlot(Bcell_clus, features = "Cr2", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD21 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD21_rna)
ggsave("CD21_rna.pdf", width = 30, height = 20, units = "cm")

CD23_rna <- FeaturePlot(Bcell_clus, features = "Fcer2a", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD23 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD23_rna)
ggsave("CD23_rna.pdf", width = 30, height = 20, units = "cm")

IgM_rna <- FeaturePlot(Bcell_clus, features = "Ighm", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("IgM mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(IgM_rna)
ggsave("IgM_rna.pdf", width = 30, height = 20, units = "cm")

IgD_rna <- FeaturePlot(Bcell_clus, features = "Ighd", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("IgD mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(IgD_rna)
ggsave("IgD_rna.pdf", width = 30, height = 20, units = "cm")

CD38_rna <- FeaturePlot(Bcell_clus, features = "Cd38", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = FALSE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD38 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD38_rna)
ggsave("CD38_rna.pdf", width = 30, height = 20, units = "cm")

CD95_rna <- FeaturePlot(Bcell_clus, features = "Fas", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD95 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD95_rna)
ggsave("CD95_rna.pdf", width = 30, height = 20, units = "cm")

CTLA4_rna <- FeaturePlot(Bcell_clus, features = "Ctla4", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CTLA4 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CTLA4_rna)
ggsave("CTLA4_rna.pdf", width = 30, height = 20, units = "cm")

PD_1_rna <- FeaturePlot(Bcell_clus, features = "Pdcd1", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("PD-1 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(PD_1_rna)
ggsave("PD_1_rna.pdf", width = 30, height = 20, units = "cm")

PD_L1_rna <- FeaturePlot(Bcell_clus, features = "Cd274", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("PD-L1 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(PD_L1_rna)
ggsave("PD_L1_rna.pdf", width = 30, height = 20, units = "cm")

PD_L2_rna <- FeaturePlot(Bcell_clus, features = "Pdcd1lg2", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("PD-L2 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(PD_L2_rna)
ggsave("PD_L2_rna.pdf", width = 30, height = 20, units = "cm")

CD80_rna <- FeaturePlot(Bcell_clus, features = "Cd80", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD80 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD80_rna)
ggsave("CD80_rna.pdf", width = 30, height = 20, units = "cm")

CD83_rna <- FeaturePlot(Bcell_clus, features = "Cd83", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD83 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD83_rna)
ggsave("CD83_rna.pdf", width = 30, height = 20, units = "cm")

CD86_rna <- FeaturePlot(Bcell_clus, features = "Cd86", reduction = "wnn.umap", cols = mako(10), pt.size = 1, order = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD86 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(CD86_rna)
ggsave("CD86_rna.pdf", width = 30, height = 20, units = "cm")

#ADT VlnPlot
DefaultAssay(Bcell_clus) <- "ADT"
Vln_B220 <- VlnPlot(Bcell_clus, features = "B220", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("B220 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_B220)
ggsave("Vln_B220.pdf", width = 30, height = 20, units = "cm")

Vln_CD93 <- VlnPlot(Bcell_clus, features = "Cd93", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD93 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD93)
ggsave("Vln_CD93.pdf", width = 30, height = 20, units = "cm")

Transitional_Bcells <- c("15", "30")

Vln_CD21 <- VlnPlot(Bcell_clus, features = "Cd21", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD21 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD21)
ggsave("Vln_CD21.pdf", width = 30, height = 20, units = "cm")

MarginalZone_Bcells <- c("1", "20")

Vln_CD23 <- VlnPlot(Bcell_clus, features = "Cd23", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD23 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD23)
ggsave("Vln_CD23.pdf", width = 30, height = 20, units = "cm")

Follicular_Bcells <- c("0", "18", "19", "21", "35")

Vln_IgM <- VlnPlot(Bcell_clus, features = "Igm", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("IgM cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_IgM)
ggsave("Vln_IgM.pdf", width = 30, height = 20, units = "cm")

IgMhi_Bcells <- c("0", "1", "15", "20")

Vln_IgD <- VlnPlot(Bcell_clus, features = "Igd", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("IgD cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_IgD)
ggsave("Vln_IgD.pdf", width = 30, height = 20, units = "cm")

IgDhi_Bcells <- c("0", "1", "15", "18", "19", "21", "35")

Vln_CD38 <- VlnPlot(Bcell_clus, features = "Cd38", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD38 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD38)
ggsave("Vln_CD38.pdf", width = 30, height = 20, units = "cm")

Vln_CD95 <- VlnPlot(Bcell_clus, features = "Cd95", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD95 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD95)
ggsave("Vln_CD95.pdf", width = 30, height = 20, units = "cm")

GC_Bcells <- c("23")

Vln_CTLA4 <- VlnPlot(Bcell_clus, features = "Ctla4", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CTLA4 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CTLA4)
ggsave("Vln_CTLA4.pdf", width = 30, height = 20, units = "cm")

Vln_PD_1 <- VlnPlot(Bcell_clus, features = "Pd-1", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("PD-1 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_PD_1)
ggsave("Vln_PD_1.pdf", width = 30, height = 20, units = "cm")

Vln_PD_L1 <- VlnPlot(Bcell_clus, features = "Pd-L1", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("PD-L1 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_PD_L1)
ggsave("Vln_PD_L1.pdf", width = 30, height = 20, units = "cm")

Vln_PD_L2 <- VlnPlot(Bcell_clus, features = "Pd-L2", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("PD-L2 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_PD_L2)
ggsave("Vln_PD_L2.pdf", width = 30, height = 20, units = "cm")

Vln_CD80 <- VlnPlot(Bcell_clus, features = "Cd80", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD80 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD80)
ggsave("Vln_CD80.pdf", width = 30, height = 20, units = "cm")

Memory_Bcells <- c("2", "10", "23")

Vln_CD83 <- VlnPlot(Bcell_clus, features = "Cd83", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD83 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD83)
ggsave("Vln_CD83.pdf", width = 30, height = 20, units = "cm")

Vln_CD86 <- VlnPlot(Bcell_clus, features = "Cd86", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD86 cell surface expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD86)
ggsave("Vln_CD86.pdf", width = 30, height = 20, units = "cm")

CD86hi_Bcells <- c("2", "23")

#RNA VlnPlot
DefaultAssay(Bcell_clus) <- "RNA"
Vln_CD93_rna <- VlnPlot(Bcell_clus, features = "Cd93", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD93 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD93_rna)
ggsave("Vln_CD93_rna.pdf", width = 30, height = 20, units = "cm")

Vln_CD21_rna <- VlnPlot(Bcell_clus, features = "Cr2", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD21 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD21_rna)
ggsave("Vln_CD21_rna.pdf", width = 30, height = 20, units = "cm")

Vln_CD23_rna <- VlnPlot(Bcell_clus, features = "Fcer2a", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD23 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD23_rna)
ggsave("Vln_CD23_rna.pdf", width = 30, height = 20, units = "cm")

Vln_IgM_rna <- VlnPlot(Bcell_clus, features = "Ighm", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("IgM mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_IgM_rna)
ggsave("Vln_IgM_rna.pdf", width = 30, height = 20, units = "cm")

Vln_IgD_rna <- VlnPlot(Bcell_clus, features = "Ighd", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("IgD mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_IgD_rna)
ggsave("Vln_IgD_rna.pdf", width = 30, height = 20, units = "cm")

Vln_CD38_rna <- VlnPlot(Bcell_clus, features = "Cd38", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD38 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD38_rna)
ggsave("Vln_CD38_rna.pdf", width = 30, height = 20, units = "cm")

Vln_CD95_rna <- VlnPlot(Bcell_clus, features = "Fas", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD95 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD95_rna)
ggsave("Vln_CD95_rna.pdf", width = 30, height = 20, units = "cm")

Vln_CTLA4_rna <- VlnPlot(Bcell_clus, features = "Ctla4", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CTLA4 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CTLA4_rna)
ggsave("Vln_CTLA4_rna.pdf", width = 30, height = 20, units = "cm")

CTLA4hi_Bcells <- c("2", "18")

Vln_PD_1_rna <- VlnPlot(Bcell_clus, features = "Pdcd1", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("PD-1 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_PD_1_rna)
ggsave("Vln_PD_1_rna.pdf", width = 30, height = 20, units = "cm")

Vln_PD_L1_rna <- VlnPlot(Bcell_clus, features = "Cd274", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("PD-L1 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_PD_L1_rna)
ggsave("Vln_PD_L1_rna.pdf", width = 30, height = 20, units = "cm")

Vln_PD_L2_rna <- VlnPlot(Bcell_clus, features = "Pdcd1lg2", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("PD-L2 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_PD_L2_rna)
ggsave("Vln_PD_L2_rna.pdf", width = 30, height = 20, units = "cm")

Vln_CD80_rna <- VlnPlot(Bcell_clus, features = "Cd80", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD80 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD80_rna)
ggsave("Vln_CD80_rna.pdf", width = 30, height = 20, units = "cm")

Vln_CD83_rna <- VlnPlot(Bcell_clus, features = "Cd83", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD83 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD83_rna)
ggsave("Vln_CD83_rna.pdf", width = 30, height = 20, units = "cm")

Vln_CD86_rna <- VlnPlot(Bcell_clus, features = "Cd86", cols = col_con2) +
  theme_bw() + NoLegend() + ggtitle("CD86 mRNA expression") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))
print(Vln_CD86_rna)
ggsave("Vln_CD86_rna.pdf", width = 30, height = 20, units = "cm")

##Summary B cell canonical markers##

Transitional_Bcells <- c("15", "30")
MarginalZone_Bcells <- c("1", "20")
Follicular_Bcells <- c("0", "18", "19", "21", "35")
IgMhi_Bcells <- c("0", "1", "15", "20")
IgDhi_Bcells <- c("0", "1", "15", "18", "19", "21", "35")
CTLA4hi_Bcells <- c("2", "18")
GC_Bcells <- c("23")
Memory_Bcells <- c("2", "10", "23")
CD86hi_Bcells <- c("2", "23")

#Change Idents on Bcell_clus#
Idents(Bcell_clus)
Bcell_clus <- RenameIdents(Bcell_clus, `0` = "IgM+ IgD+ Follicular B cells", `1` ="IgM+ IgD+ Marginal Zone B cells",
                           `2` = "CTLA4+ CD86+ Memory B cells", `10` = "Memory B cells", `15` = "IgM+ IgD+ Transitional B cells",
                           `18` = "IgD+ CTLA4+ Follicular B cells", `19` = "IgD+ Follicular B cells (1)", `20` = "IgM+ Marginal Zone B cells",
                           `21` = "IgD+ Follicular B cells (2)", `23` = "Germinal Centre B cells/CD86+ Memory B cells", `30` = "Transitonal B cells",
                           `35` = "IgD+ Follicular B cells (3)")
Bcell_clus[["Type_assigned"]] <- Idents(Bcell_clus)
Idents(Bcell_clus)

UMAP2 <- DimPlot(Bcell_clus, reduction = "wnn.umap", cols = col_con2, pt.size = 1) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("B cell Clusters CITE-Seq-2") +
  theme(plot.title = element_text(size=16, face = "bold"))
print(UMAP2)
ggsave("Bcellsl_cluster_1.pdf", width = 30, height = 20, units = "cm")

UMAP3 <- DimPlot(Bcell_clus, reduction = "wnn.umap", cols = col_con2, pt.size = 1, label = TRUE, repel = TRUE, label.box = TRUE) +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("B cell Clusters CITE-Seq-2") +
  theme(plot.title = element_text(size=16, face = "bold")) + NoLegend()
print(UMAP3)
ggsave("Bcellsl_cluster_2.pdf", width = 30, height = 20, units = "cm")

##use addmoduel score to test for cell signature enrichment##


##Use scGate to gate out popualtions##




#ADT dotplot
DefaultAssay(Bcell_clus) <- "ADT"
Abs_Bcells <- c("Cxcr5", "Cd86", "Cd80", "Cd83", "Pd-L1", "Pd-L2", "Pd-1", "Ctla4", "Cd40", "Cd44")

DotPlot(Bcell_clus, features = Abs_Bcells, cluster.idents = FALSE) + RotatedAxis() +
  scale_colour_gradient2(low = "#000004FF", mid = "#CD4071FF", high = "#FCFDBFFF")
?DotPlot





####CD4 T cell characterisation####
CD4_Tcell_gr <- WhichCells(All_cells, idents = CD4_Tcell)
CD4_Tcell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
                          label = TRUE, cells.highlight = CD4_Tcell_gr, cols.highlight = "indianred2"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4+ T cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold")) +
  scale_color_manual(labels = c("Other cells", "CD4+ T cells"), values = c("grey", "indianred2"))
print(CD4_Tcell_gr_plot)
ggsave("CD4_Tcell_gr_plot.pdf", width = 30, height = 20, units = "cm")


####CD8 T cell characterisation####
CD8_Tcell_gr <- WhichCells(All_cells, idents = CD8_Tcell)
CD8_Tcell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
                             label = TRUE, cells.highlight = CD8_Tcell_gr, cols.highlight = "turquoise4"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD8+ T cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold")) +
  scale_color_manual(labels = c("Other cells", "CD8+ T cells"), values = c("grey", "turquoise4"))
print(CD8_Tcell_gr_plot)
ggsave("CD8_Tcell_gr_plot.pdf", width = 30, height = 20, units = "cm")


####CD4 cell characterisation####
CD4_cell_gr <- WhichCells(All_cells, idents = CD4_cell)
CD4_cell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
                             label = TRUE, cells.highlight = CD4_cell_gr, cols.highlight = "mediumorchid2"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4+ cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold")) +
  scale_color_manual(labels = c("Other cells", "CD4+ cells"), values = c("grey", "mediumorchid2"))
print(CD4_cell_gr_plot)
ggsave("CD4_cell_gr_plot.pdf", width = 30, height = 20, units = "cm")

####CD8 cell characterisation####
CD8_cell_gr <- WhichCells(All_cells, idents = CD8_cell)
CD8_cell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
                            label = TRUE, cells.highlight = CD8_cell_gr, cols.highlight = "seagreen3"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD8+ cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold")) +
  scale_color_manual(labels = c("Other cells", "CD8+ cells"), values = c("grey", "seagreen3"))
print(CD8_cell_gr_plot)
ggsave("CD8_cell_gr_plot.pdf", width = 30, height = 20, units = "cm")

####CD4 T and CD8 T cell characterisation####
CD4_CD8_Tcell_gr <- WhichCells(All_cells, idents = CD4_CD8_Tcell)
CD4_CD8_Tcell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
                            label = TRUE, cells.highlight = CD4_CD8_Tcell_gr, cols.highlight = "peachpuff2"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4+ and CD8+ T cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold")) +
  scale_color_manual(labels = c("Other cells", " CD4+ and CD8+ T cells"), values = c("grey", "peachpuff2"))
print(CD4_CD8_Tcell_gr_plot)
ggsave("CD4_CD8_Tcell_gr_plot.pdf", width = 30, height = 20, units = "cm")

####CD4 and CD8 cell characterisation####
CD4_CD8_cell_gr <- WhichCells(All_cells, idents = CD4_CD8_cell)
CD4_CD8_cell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
                                 label = TRUE, cells.highlight = CD4_CD8_cell_gr, cols.highlight = "orangered3"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4+ and CD8+ cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold")) +
  scale_color_manual(labels = c("Other cells", " CD4+ and CD8+ cells"), values = c("grey", "orangered3"))
print(CD4_CD8_cell_gr_plot)
ggsave("CD4_CD8_cell_gr_plot.pdf", width = 30, height = 20, units = "cm")

####CD4 T and B cell characterisation####
CD4_Tcell_Bcell_gr <- WhichCells(All_cells, idents = CD4_Tcell_Bcell)
CD4_Tcell_Bcell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
                                label = TRUE, cells.highlight = CD4_Tcell_Bcell_gr, cols.highlight = "lightcoral"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4+ T and B cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold")) +
  scale_color_manual(labels = c("Other cells", " CD4+ T and B cells"), values = c("grey", "lightcoral"))
print(CD4_Tcell_Bcell_gr_plot)
ggsave("CD4_Tcell_Bcell_gr_plot.pdf", width = 30, height = 20, units = "cm")


####CD4 T, CD8 T and B cell characterisation####
CD4_CD8_Tcell_Bcell_gr <- WhichCells(All_cells, idents = CD4_CD8_Tcell_Bcell)
CD4_CD8_Tcell_Bcell_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
                                   label = TRUE, cells.highlight = CD4_CD8_Tcell_Bcell_gr, cols.highlight = "olivedrab3"  ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("CD4+ T, CD8+ T and B cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold")) +
  scale_color_manual(labels = c("Other cells", " CD4+ T, CD8+ T and B cells"), values = c("grey", "olivedrab3"))
print(CD4_CD8_Tcell_Bcell_gr_plot)
ggsave("CD4_CD8_Tcell_Bcell_gr_plot.pdf", width = 30, height = 20, units = "cm")


####All cells highlighted####
all_cells_gr <- c(B_cell_gr, CD4_Tcell_gr, CD8_Tcell_gr, CD4_cell_gr, CD8_cell_gr, CD4_CD8_Tcell_gr,CD4_CD8_cell_gr, CD4_Tcell_Bcell_gr, CD4_CD8_Tcell_Bcell_gr)
all_cells_cols <- c("steelblue3", "indianred2", "turquoise4", "mediumorchid2", "seagreen3", "peachpuff2", "orangered3", "lightcoral", "olivedrab3")
all_cells_cols_2 <- c("grey", "steelblue3", "indianred2", "turquoise4", "mediumorchid2", "seagreen3", "peachpuff2", "orangered3", "lightcoral", "olivedrab3")
all_cells_labels <- c("Other cells", " CD4+ T cells, CD8+ T cells, CD4+ cells, CD8+ cells, CD4+ T and CD8+ T cells, CD4+ cells and CD8+ cell, CD4+ T cells sand B cells", "CD4+, CD8+ T cells and B cells")


all_cells_gr_plot <- DimPlot(All_cells, reduction = "wnn.umap",
  label = TRUE, cells.highlight = all_cells_gr, cols.highlight = all_cells_cols ,cols = "grey") +
  theme_bw() + xlab("UMAP1") + ylab("UMAP2") + ggtitle("All cells") +
  theme(plot.title = element_text(color="black", size=16, face="bold"))





print(all_cells_gr_plot)
ggsave("all_cells_gr_plot.pdf", width = 30, height = 20, units = "cm")






