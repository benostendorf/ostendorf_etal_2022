#### Subset and basic Seurat workflow

## Define thresholds
feature_min <- 150
feature_max <- 7500
count_max <- 40000
max_perc <- 15

# Subset
df <- subset(df_raw, 
             subset = 
               nFeature_RNA > feature_min & nFeature_RNA < feature_max & 
               nCount_RNA < count_max & 
               percent.mt < max_perc)

# Seurat workflow
df <- NormalizeData(df)
df <- FindVariableFeatures(df)
df <- ScaleData(df)
df <- RunPCA(df)
df <- RunUMAP(df, dims = 1:30)
ElbowPlot(df, ndims = 30)
df <- FindNeighbors(df, dims = 1:30)
df <- FindClusters(df, resolution = 1.4)

df <- BuildClusterTree(df, reorder = TRUE, reorder.numeric = TRUE)

saveRDS(df, "../../data/scRNAseq/seurat.RDS")
