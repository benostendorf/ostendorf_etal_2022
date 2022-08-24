df_inf <- subset(df_filt, condition == "inf")
df_ctrl <- subset(df_filt, condition == "ctrl")

## df_immune
immune_subsets <- names(table(df_filt$celltype_grouped))[1:9]
df_immune <- subset(df_filt, celltype_grouped %in% immune_subsets)
df_immune$celltype_grouped <- droplevels(df_immune$celltype_grouped)
df_immune$celltype <- droplevels(df_immune$celltype)

## df_immune_inf
df_immune_inf <- subset(df_filt, celltype_grouped %in% immune_subsets & condition == "inf")
df_immune_inf$celltype_grouped <- droplevels(df_immune_inf$celltype_grouped)
df_immune_inf$celltype <- droplevels(df_immune_inf$celltype)

## df_nonimmune
nonimmune_subsets <- names(table(df_filt$celltype_grouped))[10:14]
df_nonimmune <- subset(df_filt, celltype_grouped %in% nonimmune_subsets)
df_nonimmune$celltype_grouped <- droplevels(df_nonimmune$celltype_grouped)
df_nonimmune$celltype <- droplevels(df_nonimmune$celltype)

## df_nonimmune_inf
df_nonimmune_inf <- subset(df_filt, celltype_grouped %in% nonimmune_subsets & condition == "inf")
df_nonimmune_inf$celltype_grouped <- droplevels(df_nonimmune_inf$celltype_grouped)
df_nonimmune_inf$celltype <- droplevels(df_nonimmune_inf$celltype)

## Extract UMAP data for df_filt, df_inf, and df_ctrl
extract_UMAP <- function(seurat_object){
  UMAP <- as.data.frame(Embeddings(seurat_object, reduction = "umap"))
  UMAP$genotype <- seurat_object$genotype
  UMAP$celltype <- Idents(seurat_object)
  UMAP$celltype_grouped <- seurat_object$celltype_grouped
  UMAP$condition <- seurat_object$condition
  return(UMAP)
}
UMAP_df_filt <- extract_UMAP(df_filt)
UMAP_inf <- extract_UMAP(df_inf)
UMAP_ctrl <- extract_UMAP(df_ctrl)