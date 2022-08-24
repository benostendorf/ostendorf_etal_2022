clusters <- levels(df_inf)

clusters_immune <- clusters[c(1:14)]
clusters_fibroblasts <- clusters[c(15:17)]
clusters_endothelial <- clusters[c(18:23)]
clusters_epithelial <- clusters[c(25:30)]
clusters_other <- clusters[c(24, 31)]


## =============================================================
## GSEA using hallmark pathways
## =============================================================
H_pathways_mm <- 
  msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  mutate(gs_name = gsub("HALLMARK_", "", gs_name)) %>%
  mutate(gs_name = gsub("_", " ", gs_name)) %>%
  mutate(gs_name = str_to_sentence(gs_name))

immune_hallmark <- levels(as.factor(H_pathways_mm$gs_name))[c(2, 4, 7, 10:11, 22:27, 42, 44, 45, 46, 49)]
H_pathways_mm <- filter(H_pathways_mm, gs_name %in% immune_hallmark)

for (genotype in c("E2", "E4")) {
  
  ## -----------------------
  ## Run GSEA
  ## -----------------------
  if (!(file.exists(paste0("../../data/scRNAseq/gsea_hallmark_", genotype, ".rds")))) {
    
    DGEA_all_clusters_ranked <- readRDS(paste0("../../data/scRNAseq/DGEA_geneLists_", genotype, ".rds"))
    
    gsea_all_clusters_H <- vector("list", length = length(DGEA_all_clusters_ranked))
    names(gsea_all_clusters_H) <- names(DGEA_all_clusters_ranked)
  
    for (cluster in names(DGEA_all_clusters_ranked)) {
      print(cluster)
      gsea_cluster_H <- clusterProfiler::GSEA(DGEA_all_clusters_ranked[[cluster]], 
                                              TERM2GENE = H_pathways_mm, 
                                              nPerm = 10000,
                                              pvalueCutoff = 0.1,
                                              pAdjustMethod = "BH")
  
      gsea_cluster_H <- DOSE::setReadable(gsea_cluster_H, org.Mm.eg.db, keyType = "ENTREZID")
      res_gsea_cluster_H <- as.data.frame(gsea_cluster_H)
      gsea_all_clusters_H[[cluster]] <- res_gsea_cluster_H
    }  
    
    saveRDS(gsea_all_clusters_H, paste0("../../data/scRNAseq/gsea_hallmark_", genotype, ".rds"))
  }
  gsea_all_clusters_H <- readRDS(paste0("../../data/scRNAseq/gsea_hallmark_", genotype, ".rds"))
  length(gsea_all_clusters_H)
  
  gsea_all_clusters_H_red <- gsea_all_clusters_H[sapply(gsea_all_clusters_H, nrow) > 0]
  
  ## -----------------------
  ## Wrangle for plotting
  ## -----------------------

  gsea_all_clusters_H_mut <- vector("list", length(gsea_all_clusters_H_red))
  names(gsea_all_clusters_H_mut) <- names(gsea_all_clusters_H_red)
  for (cluster in names(gsea_all_clusters_H_red)) {
    
    ## Wrangle output of GSEA for nicer plotting
    gsea_cluster_count_H <- 
      gsea_all_clusters_H_red[[cluster]] %>% 
      group_by(ID) %>% 
      summarise(count = sum(str_count(core_enrichment, "/")) + 1)
    
    gsea_cluster_H_df<- 
      left_join(gsea_all_clusters_H_red[[cluster]], gsea_cluster_count_H, by = "ID") %>% 
      mutate(GeneRatio = count/setSize ,
             type = case_when(NES < 0 ~ "Suppressed", 
                              NES > 0 ~ "Activated")) %>%
      mutate(type = ordered(type))
    
    gsea_all_clusters_H_mut[[cluster]] <- gsea_cluster_H_df
  }
  
  saveRDS(gsea_all_clusters_H_mut, paste0("../../data/scRNAseq/GSEA_H_mut_", genotype, ".RDS"))
}


clusters <- levels(df_inf$celltype_grouped)

for (genotype in c("E2", "E4")) {
  
  ls_gsea_H <- vector("list")
  
  gsea_all_clusters_H_mut <- readRDS(paste0("../../data/scRNAseq/GSEA_H_mut_", genotype, ".RDS"))
  
  for (celltype in names(gsea_all_clusters_H_mut)){
    df_for_each_cluster <- dplyr::select(gsea_all_clusters_H_mut[[celltype]], 
                                  Description, enrichmentScore, p.adjust)
    ls_gsea_H[[celltype]] <- df_for_each_cluster
  }

  df_gsea_H <- 
    bind_rows(ls_gsea_H, .id = "cluster") %>%
    pivot_wider(names_from = cluster, values_from = enrichmentScore, id_cols = Description)
  
  ## Add no hits to clusters with no significant pathway enrichment
  missing_clusters <- clusters[!clusters %in% names(ls_gsea_H)]
  missing_clusters_df <- data.frame(matrix(ncol = length(missing_clusters), 
                                           nrow = nrow(df_gsea_H)))
  rownames(missing_clusters_df) <- df_gsea_H$Description
  colnames(missing_clusters_df) <- missing_clusters
  df_gsea_h_combined <- 
    as_tibble(missing_clusters_df, rownames = "Description") |>
    full_join(df_gsea_H)
  
  df_gsea_h_combined <- as.data.frame(df_gsea_h_combined)
  rownames(df_gsea_h_combined) <- df_gsea_h_combined$Description
  
  mtx_gsea_H <- as.matrix(df_gsea_h_combined[, clusters])
  rownames(mtx_gsea_H) <- df_gsea_h_combined$Description
  mtx_gsea_H[is.na(mtx_gsea_H)] <- 0

  print(
    Heatmap(mtx_gsea_H, border = TRUE, 
            show_row_dend = FALSE, show_column_dend = FALSE, cluster_columns = FALSE, 
            row_names_gp = gpar(fontsize = 5), 
            column_names_gp = gpar(fontsize = 5), 
            col = circlize::colorRamp2(c(-1, 0, 1), 
                                           c("blue", "grey 95", "red")), 
            rect_gp = gpar(col = "white", lwd = 1), row_order = immune_hallmark, 
            heatmap_legend_param = list(grid_width = unit(0.2, "cm"), 
                                        legend_height = unit(0.4, "cm"))))
}
