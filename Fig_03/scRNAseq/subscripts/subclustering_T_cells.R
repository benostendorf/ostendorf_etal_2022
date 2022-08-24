## Re-cluster lymphocytes
T_cells <- subset(df, idents = c("T cells A", "T cells B", "T cells C"))
T_cells <- FindVariableFeatures(T_cells)
T_cells <- ScaleData(T_cells)
T_cells <- RunPCA(T_cells)
ElbowPlot(T_cells)
T_cells <- RunUMAP(T_cells, dims = 1:20)
T_cells <- FindNeighbors(T_cells, dims = 1:20)
T_cells <- FindClusters(T_cells, resolution = 0.5)


## Find DEG for lymphocytes
markers_T_cells <- FindAllMarkers(T_cells, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
 
markers_by_T_cells <- 
  markers_T_cells |>
  group_by(cluster) |>
  slice_max(n = 5, order_by = avg_log2FC)

T_cells.averages <- AverageExpression(T_cells, return.seurat = TRUE)

## ----------------------------------
## Rename clusters in T cells object
## ----------------------------------
rename_clusters_T_cells <- c(`0` = "misc_05", `1` = "T cells", 
                             `2` = "T cells naive", `3` = "Tregs", 
                             `4` = "misc_05", `5` = "misc_05", 
                             `6` = "misc_05", `7` = "T cells proliferating", 
                             `8` = "misc_05", `9` = "misc_05")

T_cells <- RenameIdents(object = T_cells, rename_clusters_T_cells)
Idents(T_cells) <- factor(Idents(T_cells), 
                     levels = levels(Idents(T_cells)), ordered = TRUE)
T_cells$celltype <- Idents(object = T_cells)

## ----------------------------------
## Rename clusters in main df object
## ----------------------------------
df$sub_cluster <- as.character(Idents(df))
df$sub_cluster[Cells(T_cells)] <- as.character(Idents(T_cells))
Idents(df) <- df$sub_cluster
df$celltype <- Idents(df)

## ----------------------------------
## Re-order clusters
## ----------------------------------
levels_clusters <- c("Alveolar mø A", "Alveolar mø B", 
                     "Alveolar mø prolif", "Interstitial mø", 
                     "Monocytes A", "Monocytes B", 
                     "Granulocytes", "DCs", 
                     "NK cells", 
                     "T cells naive", "T cells", 
                     "Tregs", "T cells proliferating", 
                     "B cells", 
                     "Myofibroblasts", "Lipofibroblasts" ,
                     "Col14a1pos fibroblasts", 
                     "Capillary ECs", "Vascular ECs A", 
                     "Vascular ECs B", "Other ECs", 
                     "Vcam1pos ECs A", "Vcam1pos ECs B", 
                     "Pericytes", 
                     "AT1", "AT2", 
                     "Ciliated cells", 
                     "Airway epithelial A", "Airway epithelial B", 
                     "Mesothelial cells", 
                     "Neuronal", 
                     "misc_01", "misc_02", "misc_03", 
                     "misc_04", "misc_05", "misc_06")

Idents(df) <- factor(Idents(df),
                     levels =  levels_clusters, 
                     ordered = TRUE)