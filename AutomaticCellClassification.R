## ##Cell type identification####
#Julius Christiopher Baeck

## ##Setup####
##set working directory##
setwd("/Volumes/GoogleDrive/Shared drives/Okkengroup/Experiments/Julius/Experiments/CITE-Sequencing/CITE-Seq (2)/Overall_analysis/CITE-Seq2_all_cells")

getOption('timeout')
options(timeout = 6000)

##Install required packages##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clustifyr")

BiocManager::install("ArrayExpress")

BiocManager::install("clustifyrdatahub")

BiocManager::install("infercnv")

BiocManager::install("devtools")

devtools::install_github("omicsCore/scTyper")

devtools::install_url("http://sourceforge.net/projects/mcmc-jags/files/rjags/3/rjags_3-2.tar.gz",
                      args="--configure-args='--with-jags-include=/Users/casallas/homebrew/opt/jags/include/JAGS        
                                              --with-jags-lib=/Users/casallas/homebrew/opt/jags/lib'
                            "
)

##Load required packages##
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

##Functions##
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  mousex <- unique(genesV2)
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

##Colour panels##
col_con <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


## ##Loading data test#####
All_cells <- LoadH5Seurat("CITE-Seq2_all_cells.h5seurat")
All_cells$seurat_clusters <- All_cells$wsnn_res.0.7
Idents(All_cells) <- All_cells$seurat_clusters

##UMAP of identified clusters in test dataset####
UMAP_1 <- DimPlot(All_cells, reduction = "wnn.umap") + theme_classic() +
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("Clustering of all cells - CITE-Seq-2") +
  theme(plot.title = element_text(size=16, face = "bold"))
ggsave("Original_cluster_1.pdf", width = 30, height = 20, units = "cm")

## ##References: Single Cell Expression Atlas####
####Tabula Muris Datasets####
#Select raw files
expression_matrix <- ReadMtx(
  mtx = "E-ENAD-15.aggregated_filtered_counts.mtx", features = "E-ENAD-15.aggregated_filtered_counts.mtx_rows",
  cells = "E-ENAD-15.aggregated_filtered_counts.mtx_cols",
)
rownames(expression_matrix)
colnames(expression_matrix)

#Convert ENSEMBL IDs to GENE Symbols
expression_df <- as.data.frame(expression_matrix) #Convert to dataframe
expression_df$Gene <- mapIds(org.Mm.eg.db, keys = rownames(expression_df), keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
expression_df_Gene <- expression_df[, "Gene", drop = FALSE] #Generate separate df

#Find NA in DF and replace with "Unknown"
is.na(expression_df_Gene)
colSums(is.na(expression_df_Gene))
expression_df_Gene[is.na(expression_df_Gene)] <- "Unknown"

#Replace rownmaes by GENE Symbols
rownames(expression_matrix) <- expression_df_Gene$Gene #Replace rownmaes by Genes
is.na(expression_matrix)

#Generate Seurat Object
seurat_object <- CreateSeuratObject(counts = expression_matrix)
head(seurat_object[[]])

#Addd metadata to Seurat Object
meta.data <- read_tsv(file.choose("E-ENAD-15.clusters.tsv"))
meta.data <- t(meta.data) #Columns to rows
meta.data <- as.data.frame(meta.data[-(1:2), 5]) #as data frame + remove first two rows
colnames(meta.data) <- c("K47")

head(meta.data)
rownames(meta.data)
seurat_object_comp <- AddMetaData(seurat_object, meta.data)
head(seurat_object_comp[[]])

SaveH5Seurat(seurat_object_comp, filename = "TabulaMuris", overwrite = TRUE)
seurat_object_comp <- LoadH5Seurat(file.choose("TabulaMuris.h5seurat"))

#Normalise and scale seurat object
seurat_object_comp <- NormalizeData(seurat_object_comp)
seurat_object_comp <- FindVariableFeatures(seurat_object_comp)
seurat_object_comp <- ScaleData(seurat_object_comp)
seurat_object_comp
SaveH5Seurat(seurat_object_comp, filename = "TabulaMuris_norm", overwrite = TRUE)
head(seurat_object_comp[[]])


#Generate Reference for clustify
s_ref <- seurat_ref(
  seurat_object = seurat_object_comp,
  cluster_col = "K47"
)
head(s_ref)

#Run clustifyr
DefaultAssay(All_cells) <- "RNA"

res <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = s_ref,
  obj_out = TRUE,
  dr = "wnn.umap"
)
head(res[[]])

#Visualise matching cell type and correlation
DimPlot(res, reduction = "wnn.umap", group.by = "type")
FeaturePlot(res, reduction = "wnn.umap", features = "r", cols = viridis(10))

#Rename clusters accroding to literature
Idents(res) <- res$type
res <- RenameIdents(res, `3` = "B cells", `8` ="T cells", `14` ="Progenitor cells", `15` = "Granulocyte and hematopoetic precursor cells",
                    `26` = "Myeloid cells",`34` ="Myeloid cells 2", `42` = "Monocytes") #Based on Single Cell Expression Atlas website
res[["Type_assigned"]] <- Idents(res)
Idents(res) <- res$seurat_clusters

DimPlot(res, reduction = "wnn.umap", group.by = "Type_assigned")

#Add meta.data_2 from publication
TabMur_droplet <- LoadH5Seurat(file.choose("droplet.seurat.h5Seurat"))
reference_1 <- TabMur_droplet
head(reference_1[[]])

#Normalise and scale data
DefaultAssay(reference_1) <- "originalexp"

reference_1 <- NormalizeData(reference_1)
reference_1 <- FindVariableFeatures(reference_1)
reference_1 <- ScaleData(reference_1)
reference_1

#Make reference matrix
Idents(reference_1) <- reference_1$cell_ontology_class
s_ref_1 <- AverageExpression(object = reference_1, assay = "originalexp")
head(s_ref_1)
s_ref_1 <- as.data.frame(s_ref_1)
head(s_ref_1)

#Chlean up column names
colnames(s_ref_1) <- str_sub(colnames(s_ref_1), 13, -1)
colnames(s_ref_1) <- str_to_title(colnames(s_ref_1))
colnames(s_ref_1) <- gsub(colnames(s_ref_1), pattern = "\\.", replacement = " ")
names(s_ref_1)[names(s_ref_1) == 'Immature t cell'] <- "Immature T cell"



#Run clustifyr
res_1 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = s_ref_1,
  obj_out = TRUE,
  dr = "wnn.umap"
)

#Visualise matching cell type and correlation
TabMur_plot1 <- DimPlot(res_1, reduction = "wnn.umap", group.by = "type", cols = col_con)
TabMur_plot1 <- LabelClusters(TabMur_plot1, id = "type",  fontface = "bold", color = "black") +
  theme_classic() + NoLegend() + ggtitle("Automatic Cell Annotation - Clustifyr - Ref: Tabula Muris") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=22, face = "bold"))
ggsave("TabMur_plot1.pdf", width = 30, height = 20, units = "cm")

TabMur_plot2 <-FeaturePlot(res_1, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Confidence (r) Score - Clustifyr - Ref: Tabula Muris") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=22, face = "bold"))
ggsave("TabMur_plot2.pdf", width = 30, height = 20, units = "cm")


####CD4+ Treg and Tmem####
expression_matrix_2 <- ReadMtx(
  mtx = "E-MTAB-7311.aggregated_filtered_counts.mtx", features = "E-MTAB-7311.aggregated_filtered_counts.mtx_rows",
  cells = "E-MTAB-7311.aggregated_filtered_counts.mtx_cols",
)

rownames(expression_matrix_2)
colnames(expression_matrix_2)

#Generate DF for changing row- and colnames
expression_df_2 <- as.data.frame(expression_matrix_2) #Convert to dataframe

#Convert ENSEMBL IDs to GENE Symbols
expression_df_2$Gene <- mapIds(org.Mm.eg.db, keys = rownames(expression_df_2), keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
expression_df_2_Gene <- expression_df_2[, "Gene", drop = FALSE] #Generate separate df
head(expression_df_2_Gene)

#Find NA in DF and replace with "Unknown"
is.na(expression_df_2_Gene)
colSums(is.na(expression_df_2_Gene)) #2309
expression_df_2_Gene[is.na(expression_df_2_Gene)] <- "Unknown"

#Replace rownmaes by GENE Symbols
rownames(expression_matrix_2) <- expression_df_2_Gene$Gene #Replace rownmaes by Genes
is.na(expression_matrix_2)

#Change colnames to just have barcodes
colnames(expression_matrix_2) <- str_sub(colnames(expression_matrix_2), 14, -1)
sum(duplicated(colnames(expression_matrix_2))) #844 of 35,470 cells

sum(duplicated(rownames(expression_matrix_2))) #2340 duplicated features
sum(rownames(expression_matrix_2) == "Unknown") #2309 are from unknown genes

#Generate Seurat Object
seurat_object_2 <- CreateSeuratObject(counts = expression_matrix_2) #males all row- and colnames unique
head(seurat_object_2[[]])

#Load metadata_1
CD4_meta <- read.csv(file.choose("all_data_cl_metadata.csv"))
head(CD4_meta)

#Change rownames to barcodes and remove duplications
CD4_meta$Column1 <- str_sub(CD4_meta$Column1,17,-3)
sum(duplicated(CD4_meta$Column1)) #801
CD4_meta <- CD4_meta[!duplicated(CD4_meta$Column1), ]
sum(duplicated(CD4_meta$Column1)) #No duplicated barcodes
rownames(CD4_meta) <- CD4_meta$Column1
rownames(CD4_meta)

seurat_object_comp_2 <- AddMetaData(seurat_object_2, CD4_meta)
head(seurat_object_comp_2[[]])

#Load metadata_2
CD4_meta_2 <- read_tsv(file.choose("E-MTAB-7311.clusters.tsv"))
CD4_meta_2 <- t(CD4_meta_2) #Columns to rows
head(CD4_meta_2)
CD4_meta_2 <- as.data.frame(CD4_meta_2[-(1:2), 5]) #as data frame + remove first two rows
colnames(CD4_meta_2) <- c("K21")
CD4_meta_2 <- tibble::rownames_to_column(CD4_meta_2, "Barcode") #Rownames to column
CD4_meta_2$Barcode <- str_sub(CD4_meta_2$Barcode ,14,-1) #Change to barcode
sum(duplicated(CD4_meta_2$Barcode)) #844 duplicates
CD4_meta_2 <- CD4_meta_2[!duplicated(CD4_meta_2$Barcode), ] #Remove duplicates
sum(duplicated(CD4_meta_2$Barcode)) #0 duplicates
rownames(CD4_meta_2) <- CD4_meta_2$Barcode #Make barcodes rownames
head(rownames(CD4_meta_2))

#Add metadata_2 to Seurat object - seurat_object_comp
seurat_object_comp_2 <- AddMetaData(seurat_object_comp_2, CD4_meta_2)
head(seurat_object_comp_2[[]])

#Normalise and scale seurat object
seurat_object_comp_2 <- NormalizeData(seurat_object_comp_2)
seurat_object_comp_2 <- FindVariableFeatures(seurat_object_comp_2)
seurat_object_comp_2 <- ScaleData(seurat_object_comp_2)
seurat_object_comp_2

#Generate Reference for clustify
head(seurat_object_comp_2)
s_ref_2 <- seurat_ref(
  seurat_object = seurat_object_comp_2,
  cluster_col = "cl_annot" #or other column
)
head(s_ref_2)

#Clean up column name
colnames(s_ref_2) <- gsub(colnames(s_ref_2), pattern = "\\_", replacement = " ")

#Run clustifyr
DefaultAssay(All_cells) <- "RNA"

res_2 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = s_ref_2,
  obj_out = TRUE,
  dr = "wnn.umap"
)
head(res_2[[]])

#Visualise matching cell type and correlation
CD4_plot1 <- DimPlot(res_2, reduction = "wnn.umap", group.by = "type", cols = col_con, repel = TRUE)
CD4_plot1 <- LabelClusters(CD4_plot1, id = "type",  fontface = "bold", color = "black") +
  theme_classic() + ggtitle("Automatic Cell Annotation - Clustifyr - Ref: CD4+ Treg and Tmem cells (Miragaia et al., 2019)") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
ggsave("CD4_plot1.pdf", width = 30, height = 20, units = "cm")

CD4_plot2 <- FeaturePlot(res_2, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr - Ref: CD4+ Treg and Tmem cells (Miragaia et al., 2019)") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
ggsave("CD4_plot2.pdf", width = 30, height = 20, units = "cm")


#Labeling of K21 clusters 
marker_genes_21 <- read_tsv(file.choose("E-MTAB-7311,marker_gnes_21"))
head(marker_genes_21)

#Change ENSEMBL IDs of marker genes to GENE symbols
marker_genes_21_df <- as.data.frame(marker_genes_21) #Convert to dataframe
marker_genes_21_df$Symbols <- mapIds(org.Mm.eg.db, keys = marker_genes_21_df$genes, keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
head(marker_genes_21_df)

#Filter out top genes per cluster 
marker_genes_21_df <- marker_genes_21_df %>%
  arrange(cluster, desc(logfoldchanges))
write.csv(marker_genes_21_df, "CD4_marker_21.csv" ,row.names = TRUE)

####Hematopoietic stem and progenitor cells####
expression_matrix_3 <- ReadMtx(
  mtx = "E-GEOD-81682.aggregated_filtered_counts.mtx", features = "E-GEOD-81682.aggregated_filtered_counts.mtx_rows",
  cells = "E-GEOD-81682.aggregated_filtered_counts.mtx_cols",
)

rownames(expression_matrix_3)
colnames(expression_matrix_3)

#Generate DF for changing row- and colnames
expression_df_3 <- as.data.frame(expression_matrix_3) #Convert to dataframe

#Convert ENSEMBL IDs to GENE Symbols
expression_df_3$Gene <- mapIds(org.Mm.eg.db, keys = rownames(expression_df_3), keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
expression_df_3_Gene <- expression_df_3[, "Gene", drop = FALSE] #Generate separate df
head(expression_df_3_Gene)

#Find NA in DF and replace with "Unknown"
is.na(expression_df_3_Gene)
colSums(is.na(expression_df_3_Gene)) #4689
expression_df_3_Gene[is.na(expression_df_3_Gene)] <- "Unknown"

#Replace rownmaes by GENE Symbols
rownames(expression_matrix_3) <- expression_df_3_Gene$Gene #Replace rownmaes by Genes
is.na(expression_matrix_3)

#Generate Seurat Object
seurat_object_3 <- CreateSeuratObject(counts = expression_matrix_3)
head(seurat_object_3[[]])

#Addd metadata to Seurat Object
meta.data_3 <- read_tsv(file.choose("E-GEOD-81682.clusters.tsv"))
meta.data_3 <- t(meta.data_3) #Columns to rows
meta.data_3 <- as.data.frame(meta.data_3[-(1:2), 4]) #as data frame + remove first two rows
colnames(meta.data_3) <- c("K8")

head(meta.data_3)
rownames(meta.data_3)
seurat_object_comp_3 <- AddMetaData(seurat_object_3, meta.data_3)
head(seurat_object_comp_3[[]])

#Normalise and scale seurat object
seurat_object_comp_3 <- NormalizeData(seurat_object_comp_3)
seurat_object_comp_3 <- FindVariableFeatures(seurat_object_comp_3)
seurat_object_comp_3 <- ScaleData(seurat_object_comp_3)
seurat_object_comp_3

#Generate Reference for clustify
s_ref_3 <- seurat_ref(
  seurat_object = seurat_object_comp_3,
  cluster_col = "K8"
)
head(s_ref_3)

#Run clustifyr
DefaultAssay(All_cells) <- "RNA"

res_3 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = s_ref_3,
  obj_out = TRUE,
  dr = "wnn.umap"
)
head(res_3[[]])

#Visualise matching cell type and correlation
DimPlot(res_3, reduction = "wnn.umap", group.by = "type")
FeaturePlot(res, reduction = "wnn.umap", features = "r", cols = viridis(10))

#Labeling of K8 clusters 
marker_genes3_8 <- read_tsv(file.choose("E-GEOD-81682.marker_genes_8.tsv"))
head(marker_genes3_8)

#Change ENSEMBL IDs of marker genes to GENE symbols
marker_genes3_8_df <- as.data.frame(marker_genes3_8) #Convert to dataframe
marker_genes3_8_df$Symbols <- mapIds(org.Mm.eg.db, keys = marker_genes3_8_df$genes, keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
head(marker_genes3_8_df)

#Filter out top genes per cluster 
marker_genes3_8_df <- marker_genes3_8_df %>%
  arrange(cluster, desc(logfoldchanges))
write.csv(marker_genes3_8_df, "marker_genes3_8_df.csv" ,row.names = TRUE)


#Rename clusters according to publication and associated website
Idents(res_3) <- res_3$type
res_3 <- RenameIdents(res_3, `1` ="Hematopoietic Stem and Progenitor - Lymphoid Multipotent Progenitors",`4` = "Progenitor - Common Myeloid Progenitor and/or Granulocyte-Monocyte Progenitor",
                    `5` ="Progenitor - Megakaryocyte-Erythrocyte Progenitor",
                    `8` ="B cell or B cell Precursor (self-determined)")
res_3[["Type_assigned"]] <- Idents(res_3)
Idents(res_3) <- res_3$seurat_clusters

DimPlot(res_3, reduction = "wnn.umap", group.by = "Type_assigned")

HSC_plot1 <- DimPlot(res_3, reduction = "wnn.umap", group.by = "Type_assigned", cols = col_con, repel = TRUE)
HSC_plot1 <- LabelClusters(HSC_plot1, id = "Type_assigned",  fontface = "bold", color = "black") +
  theme_classic() + NoLegend() + ggtitle("Automatic Cell Annotation - Clustifyr - Ref: HSC and HSC Progenitors (Nestorowa et al., 2016)") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
ggsave("HSC_plot1.pdf", width = 30, height = 20, units = "cm")

HSC_plot2 <- FeaturePlot(res_3, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr - Ref: HSC and HSC Progenitors (Nestorowa et al., 2016)") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
ggsave("HSC_plot2.pdf", width = 30, height = 20, units = "cm")


####Denritic cells####
expression_matrix_4 <- ReadMtx(
  mtx = "E-MTAB-10196.aggregated_filtered_counts.mtx", features = "E-MTAB-10196.aggregated_filtered_counts.mtx_rows",
  cells = "E-MTAB-10196.aggregated_filtered_counts.mtx_cols",
)

rownames(expression_matrix_4)
colnames(expression_matrix_4)

#Generate DF for changing row- and colnames
expression_df_4 <- as.data.frame(expression_matrix_4) #Convert to dataframe

#Convert ENSEMBL IDs to GENE Symbols
expression_df_4$Gene <- mapIds(org.Mm.eg.db, keys = rownames(expression_df_4), keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
expression_df_4_Gene <- expression_df_4[, "Gene", drop = FALSE] #Generate separate df
head(expression_df_4_Gene)

#Find NA in DF and replace with "Unknown"
is.na(expression_df_4_Gene)
colSums(is.na(expression_df_4_Gene)) #2162
expression_df_4_Gene[is.na(expression_df_4_Gene)] <- "Unknown"

#Replace rownmaes by GENE Symbols
rownames(expression_matrix_4) <- expression_df_4_Gene$Gene #Replace rownmaes by Genes
is.na(expression_matrix_4)

#Generate Seurat Object
seurat_object_4 <- CreateSeuratObject(counts = expression_matrix_4)
head(seurat_object_4[[]])

#Addd metadata to Seurat Object
meta.data_4 <- read_tsv(file.choose("E-MTAB-10196.clusters.tsv"))
meta.data_4 <- t(meta.data_4) #Columns to rows
meta.data_4 <- as.data.frame(meta.data_4[-(1:2), 5]) #as data frame + remove first two rows
colnames(meta.data_4) <- c("K19")

head(meta.data_4)
rownames(meta.data_4)
seurat_object_comp_4 <- AddMetaData(seurat_object_4, meta.data_4)
head(seurat_object_comp_4[[]])

#Normalise and scale seurat object
seurat_object_comp_4 <- NormalizeData(seurat_object_comp_4)
seurat_object_comp_4 <- FindVariableFeatures(seurat_object_comp_4)
seurat_object_comp_4 <- ScaleData(seurat_object_comp_4)
seurat_object_comp_4

#Generate Reference for clustify
s_ref_4 <- seurat_ref(
  seurat_object = seurat_object_comp_4,
  cluster_col = "K19"
)
head(s_ref_4)

#Run clustifyr
DefaultAssay(All_cells) <- "RNA"

res_4 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = s_ref_4,
  obj_out = TRUE,
  dr = "wnn.umap"
)
head(res_4[[]])

#Visualise matching cell type and correlation
DimPlot(res_4, reduction = "wnn.umap", group.by = "type")
FeaturePlot(res, reduction = "wnn.umap", features = "r", cols = viridis(10))

#Labeling of K19 clusters 
marker_genes4_19 <- read_tsv(file.choose("E-MTAB-10196.marker_genes_19.tsv"))
head(marker_genes4_19)

#Change ENSEMBL IDs of marker genes to GENE symbols
marker_genes4_19_df <- as.data.frame(marker_genes4_19) #Convert to dataframe
marker_genes4_19_df$Symbols <- mapIds(org.Mm.eg.db, keys = marker_genes4_19_df$genes, keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
head(marker_genes4_19_df)

#Filter out top genes per cluster 
marker_genes4_19_df <- marker_genes4_19_df %>%
  arrange(cluster, desc(logfoldchanges))
write.csv(marker_genes4_19_df, "marker_genes4_19_df.csv" ,row.names = TRUE)

#Rename clusters according to literature
Idents(res_4) <- res_4$type
res_4 <- RenameIdents(res_4, `2` ="Migratory DC",`3` = "CD11b+ DC/Conventional DC2",
                    `19` ="Plasmacytoid DC") #Based on Single Cell Expression Atlas website
res_4[["Type_assigned"]] <- Idents(res_4)
Idents(res_4) <- res_4$seurat_clusters

DC_plot1 <- DimPlot(res_4, reduction = "wnn.umap", group.by = "Type_assigned", cols = col_con, repel = TRUE)
DC_plot1 <- LabelClusters(DC_plot1, id = "Type_assigned",  fontface = "bold", color = "black") +
  theme_classic() + NoLegend() + ggtitle("Automatic Cell Annotation - Clustifyr - Ref: Dendritic Cells (Kapoor et al., 2021)") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=20, face = "bold"))
ggsave("DC_plot1.pdf", width = 30, height = 20, units = "cm")

DC_plot2 <- FeaturePlot(res_4, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr - Ref: Dendritic Cells (Kapoor et al., 2021)") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
ggsave("DC_plot2.pdf", width = 30, height = 20, units = "cm")

## ##References:Seurat objects#####
##Mouse BoneMarrow
load(file.choose("NicheData10x.rda"), verbose = TRUE)
reference_5 <- NicheData10x
Idents(reference_5) <- levels(NicheData10x)
reference_5[["Clusters"]] <- Idents(reference_5)
Idents(reference_5) <- levels(NicheData10x)
head(reference_5[[]])

reference_5 <- NormalizeData(reference_5)
reference_5 <- FindVariableFeatures(reference_5)
reference_5 <- ScaleData(reference_5)
reference_5

#Make reference matrix 1.0
s_ref_5 <- seurat_ref(
  seurat_object = reference_5,
  cluster_col = "Clusters",
  assay_name = "RNA"
)
head(s_ref_5)


#Make reference matrix 2.0
s_ref_5.1 <- AverageExpression(object = reference_5, assay = "RNA")
head(s_ref_5.1)
s_ref_5.1 <- as.data.frame(s_ref_5.1)
head(s_ref_5.1)

res_5 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = s_ref_5.1,
  obj_out = TRUE,
  dr = "wnn.umap"
)
head(res_5)

#Visualise reults
DimPlot(res_5, reduction = "wnn.umap", group.by = "type")#Voelliger Quatsch!
FeaturePlot(res_3, reduction = "wnn.umap", features = "r", cols = viridis(10))


####B cells immunised with SRBCs####
expression_matrix_6 <- ReadMtx(
  mtx = "E-CURD-117.aggregated_filtered_counts.mtx", features = "E-CURD-117.aggregated_filtered_counts.mtx_rows",
  cells = "E-CURD-117.aggregated_filtered_counts.mtx_cols",
)

rownames(expression_matrix_6)
colnames(expression_matrix_6)

#Generate DF for changing row- and colnames
expression_df_6 <- as.data.frame(expression_matrix_6) #Convert to dataframe

#Convert ENSEMBL IDs to GENE Symbols
expression_df_6$Gene <- mapIds(org.Mm.eg.db, keys = rownames(expression_df_6), keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
expression_df_6_Gene <- expression_df_6[, "Gene", drop = FALSE] #Generate separate df
head(expression_df_6_Gene)

#Find NA in DF and replace with "Unknown"
is.na(expression_df_6_Gene)
colSums(is.na(expression_df_6_Gene)) #1473
expression_df_6_Gene[is.na(expression_df_6_Gene)] <- "Unknown"

#Replace rownmaes by GENE Symbols
rownames(expression_matrix_6) <- expression_df_6_Gene$Gene #Replace rownmaes by Genes
is.na(expression_matrix_6)

#Generate Seurat Object
seurat_object_6 <- CreateSeuratObject(counts = expression_matrix_6)
head(seurat_object_6[[]])

#Addd metadata to Seurat Object
meta.data_6 <- read_tsv(file.choose("E-CURD-117.clusters.tsv"))
meta.data_6 <- t(meta.data_6) #Columns to rows
meta.data_6 <- as.data.frame(meta.data_6[-(1:2), 5]) #as data frame + remove first two rows
colnames(meta.data_6) <- c("K24")

head(meta.data_6)
rownames(meta.data_6)
seurat_object_comp_6 <- AddMetaData(seurat_object_6, meta.data_6)
head(seurat_object_comp_6[[]])

#Normalise and scale seurat object
seurat_object_comp_6 <- NormalizeData(seurat_object_comp_6)
seurat_object_comp_6 <- FindVariableFeatures(seurat_object_comp_6)
seurat_object_comp_6 <- ScaleData(seurat_object_comp_6)
seurat_object_comp_6

#Generate Reference for clustify
s_ref_6 <- seurat_ref(
  seurat_object = seurat_object_comp_6,
  cluster_col = "K24"
)
head(s_ref_6)

#Run clustifyr
DefaultAssay(All_cells) <- "RNA"

res_6 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = s_ref_6,
  obj_out = TRUE,
  dr = "wnn.umap"
)
head(res_6[[]])

#Visualise matching cell type and correlation
DimPlot(res_6, reduction = "wnn.umap", group.by = "type")
FeaturePlot(res_6, reduction = "wnn.umap", features = "r", cols = viridis(10))

#Labeling of K24 clusters 
marker_genes6_24 <- read_tsv(file.choose("E-CURD-117.marker_genes_24.tsv"))
head(marker_genes6_24)

#Change ENSEMBL IDs of marker genes to GENE symbols
marker_genes6_24_df <- as.data.frame(marker_genes6_24) #Convert to dataframe
marker_genes6_24_df$Symbols <- mapIds(org.Mm.eg.db, keys = marker_genes6_24_df$genes, keytype = "ENSEMBL", column = "SYMBOL") #Replace ENSMBL IDs with GENE symbols
head(marker_genes6_24_df)

#Filter out top genes per cluster 
marker_genes6_24_df <- marker_genes6_24_df %>%
  arrange(cluster, desc(logfoldchanges))
write.csv(marker_genes6_24_df, "marker_genes6_24_df.csv" ,row.names = TRUE)
#!Not enough information on cluaster-associated markers given in publication, to deterimne cell types!
#Also, no detailed B cell characterisation!

Sple_plot1 <- DimPlot(res_6, reduction = "wnn.umap", group.by = "type", cols = col_con, repel = TRUE)
Sple_plot1 <- LabelClusters(Sple_plot1, id = "type",  fontface = "bold", color = "black") +
  theme_classic() + NoLegend() + ggtitle("Automatic Cell Annotation - Clustifyr - Ref: Splenocytes (Zhai et al., 2022)") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=20, face = "bold"))
Sple_plot1
ggsave("Sple_plot1.pdf", width = 30, height = 20, units = "cm")

Sple_plot2 <- FeaturePlot(res_6, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr - Ref: Splenocytes (Zhai et al., 2022)") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
ggsave("Sple_plot2.pdf", width = 30, height = 20, units = "cm")

## ##Clustifyr datahub####
eh <- ExperimentHub()

##query##
clus_refs <- query(eh, "clustifyrdatahub")
clus_refs

####Mouse Cell Atlas####
ref_MCA <- clus_refs[[1]] #Mouse cell atlas
ref_MCA <- clus_refs[["EH3444"]]  #Mouse cell atlas
head(ref_MCA)
str(ref_MCA)

##Run classification##
res_7 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = ref_MCA,
  obj_out = TRUE,
  dr = "wnn.umap"
)

MCA_plot1 <- DimPlot(res_7, reduction = "wnn.umap", group.by = "type", cols = col_con)
MCA_plot1 <- MCA_plot1 + theme_classic() + ggtitle("Automatic Cell Annotation - Clustifyr Ref: Mouse Cell Atlas") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
MCA_plot1
ggsave("MCA_plot1.pdf", width = 30, height =16, units = "cm")

MCA_plot2 <- FeaturePlot(res_7, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr Ref: Mouse Cell Atlas") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
MCA_plot2
ggsave("MCA_plot2.pdf", width = 30, height = 20, units = "cm")

####Tabula Muris Drop-Seq####
ref_TM <- clus_refs[["EH3445"]]

res_8 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = ref_TM,
  obj_out = TRUE,
  dr = "wnn.umap"
)

TM_plot1 <- DimPlot(res_8, reduction = "wnn.umap", group.by = "type", cols = col_con)
TM_plot1 <- TM_plot1 + theme_classic() + ggtitle("Automatic Cell Annotation - Clustifyr Ref: Tabula Muris Drop-Seq") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
TM_plot1
ggsave("TM_plot1.pdf", width = 30, height =16, units = "cm")

TM_plot2 <- FeaturePlot(res_8, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr Ref: Tabula Muris Drop-Seq") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
TM_plot2
ggsave("TM_plot2.pdf", width = 30, height = 20, units = "cm")

####Tabula Muris FACS####
ref_TM2 <- clus_refs[["EH3446"]]

res_9 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = ref_TM2,
  obj_out = TRUE,
  dr = "wnn.umap"
)

TM2_plot1 <- DimPlot(res_9, reduction = "wnn.umap", group.by = "type", cols = col_con)
TM2_plot1 <- TM2_plot1 + theme_classic() + ggtitle("Automatic Cell Annotation - Clustifyr Ref: Tabula Muris - FACS") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
TM2_plot1
ggsave("TM2_plot1.pdf", width = 30, height =16, units = "cm")

TM2_plot2 <- FeaturePlot(res_9, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr Ref: Mouse Cell Atlas") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
TM2_plot2
ggsave("TM2_plot2.pdf", width = 30, height = 20, units = "cm")

####Mouse cell types (mouse rna-seq)####
ref_celltypes <- clus_refs[["EH3447"]]

res_10 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = ref_celltypes,
  obj_out = TRUE,
  dr = "wnn.umap"
)

CT_plot1 <- DimPlot(res_10, reduction = "wnn.umap", group.by = "type", cols = col_con)
CT_plot1 <- CT_plot1 + theme_classic() + ggtitle("Automatic Cell Annotation - Clustifyr Ref: 28 Murine Cell Types - RNA-Seq") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
CT_plot1
ggsave("CT_plot1.pdf", width = 30, height =16, units = "cm")

CT_plot2 <- FeaturePlot(res_10, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr Ref: 28 Murine Cell Types - RNA-Seq") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
CT_plot2
ggsave("CT_plot2.pdf", width = 30, height = 20, units = "cm")

####Mouse organogenesis####
ref_org <- clus_refs[["EH3448"]]

res_11 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = ref_org,
  obj_out = TRUE,
  dr = "wnn.umap"
)

Org_plot1 <- DimPlot(res_11, reduction = "wnn.umap", group.by = "type", cols = col_con)
Org_plot1 <- Org_plot1 + theme_classic() + ggtitle("Automatic Cell Annotation - Clustifyr Ref: Mouse Organogenesis Cell Atlas") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
Org_plot1
ggsave("Org_plot1.pdf", width = 30, height =16, units = "cm")

Org_plot2 <- FeaturePlot(res_11, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr Ref: Mouse Organogenesis Cell Atlas") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
Org_plot2
ggsave("Org_plot2.pdf", width = 30, height = 20, units = "cm")

####Mouse sorted immune cells####
ref_imm <- clus_refs[["EH3449"]]

res_12 <- clustify(
  input = All_cells,
  cluster_col = "seurat_clusters",
  ref_mat = ref_imm,
  obj_out = TRUE,
  dr = "wnn.umap"
)

Imm_plot1 <- DimPlot(res_12, reduction = "wnn.umap", group.by = "type", cols = col_con)
Imm_plot1 <- Imm_plot1 + theme_classic() + ggtitle("Automatic Cell Annotation - Clustifyr Ref: Mouse Sorted Immune Cells") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
Imm_plot1
ggsave("Imm_plot1.pdf", width = 30, height =16, units = "cm")

Imm_plot2 <- FeaturePlot(res_12, reduction = "wnn.umap", features = "r", cols = viridis(10)) +
  theme_classic() + ggtitle("Automatic Cell Annotation Comfidence (r) Score - Clustifyr Ref: Mouse Sorted Immune Cells") + xlab("UMAP1") +
  ylab("UMAP2") + theme(plot.title = element_text(size=16, face = "bold"))
Imm_plot2
ggsave("Imm_plot2.pdf", width = 30, height = 20, units = "cm")


####Mouse cell atlas 2####
ref_atlas <- clus_refs[[11]] #File might be corrupted or non-existant anymore

## ##Cell Marker Identification - with addmodulescore#####
##Check and prepare Seurat object##
Idents(All_cells)
Test_Seurat <- All_cells

##Load cell marker data##
Cell_Marker_data <- read.csv(file.choose("Cell_Markers.csv"),header = T, sep = ',')
head(Cell_Marker_data)
str(Cell_Marker_data)

##Convert cell marker df into a list##
marker_list <- list() #Create an empty list

for(i in 1:ncol(Cell_Marker_data)) {
  marker_list[[i]] <- Cell_Marker_data[ , i]
} # Using for-loop to add columns to list

names(marker_list) <- colnames(Cell_Marker_data) #Add colnames to list elements

marker_list <- lapply(marker_list, str_to_title) #Capitalise gene symbols

head(marker_list) #Print marker_list
names(marker_list)

##Run addmodule score##
DefaultAssay(Test_Seurat) <- "RNA"

Test_Seurat <- AddModuleScore(Test_Seurat, features = marker_list[c(1,3,5)], name = "Couple", search = FALSE)
colnames(Test_Seurat[[]])
VlnPlot(Test_Seurat, c("Couple1", "Couple2", "Couple3"), pt.size  = 0.5) + NoLegend()
FeaturePlot(Test_Seurat, c("Couple1", "Couple2", "Couple3"), cols=viridis(10), reduction = "wnn.umap",pt.size  = 1.5, order = FALSE)

?FeaturePlot

## ##References: ArrayExpress or BioStudies data####
#Splenocytes
rawset <- ArrayExpress("E-MTAB-9769") #No raw files
processedset <- getAE("E-MTAB-9769") #No processed files

#Bone Marrow datset
mexp8630 <- getAE("E-MTAB-8630", type = "full")
MEXP8630raw <- ae2bioc(mageFiles = mexp8630) #No raw files

mexp8630 <- getAE("E-MTAB-8630", type = "full", extract = TRUE)
cn = getcolproc(mexp8630)
MEXP8630 <- procset(mexp8630, cn[2]) #which files to select?

## ##References: GEO data####











