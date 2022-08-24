## Load processed Seurat object
df <- readRDS("../../data/scRNAseq/seurat_processed.RDS")

## Define clusters to filter
clusters_to_filter <- grep("misc", levels(Idents(df)), value = TRUE)

## Remove ambiguous clusters
df_filt <- subset(df, idents = clusters_to_filter, 
                  invert = TRUE)
df_filt$celltype <- droplevels(df_filt$celltype)
df_filt.averages <- AverageExpression(df_filt, return.seurat = TRUE)
celltypes <- as.character(Idents(df_filt))

## Group cell types
df_filt$celltype_grouped <- case_when(grepl("Monocyte", celltypes) ~ "Monocytes",
                                      grepl("Alveolar", celltypes) ~ "Alveolar Mø",
                                      grepl("T cells", celltypes) ~ "T cells",
                                      grepl("Fibro", celltypes, ignore.case = T) ~ "Fibroblasts",
                                      grepl("ECs", celltypes) ~ "Endothelial",
                                      grepl("AT1|AT2|epith|Meso|Cilia", celltypes) ~ "Epithelium", 
                                      TRUE ~ celltypes)
df_filt$celltype_grouped <- 
  factor(df_filt$celltype_grouped, 
         levels = c("Alveolar Mø", "Interstitial mø", "Monocytes", "Granulocytes", 
                    "DCs", "NK cells", "T cells", "Tregs", 
                    "B cells", "Endothelial", "Pericytes", "Fibroblasts", 
                    "Epithelium", "Neuronal"), 
         ordered = TRUE)
saveRDS(df_filt, "../../data/scRNAseq/df_filt.RDS")