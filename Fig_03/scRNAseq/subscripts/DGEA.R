## Calculate DGEA for infected cells across genotypes

## Include all genes by setting parameters to `-Inf` for gene ranking
Idents(df_inf) <- df_inf$celltype_grouped

for (genotype in c("E2", "E4")) {

  if (!(file.exists(paste0("../../data/scRNAseq/DGEA_geneLists_", genotype, ".rds")))) {

    if (!(file.exists(paste0("../../data/scRNAseq/DGEA_", genotype, ".rds")))) {

      DGEA_all_clusters <- vector("list", length = nlevels(df_inf))
      names(DGEA_all_clusters) <- levels(df_inf)

      for (cluster in levels(df_inf)) {
        print(cluster)
        DGEA_all_clusters[[cluster]] <-
          FindMarkers(df_inf,
                      group.by = "genotype",
                      subset.ident = cluster,
                      ident.1 = genotype,
                      ident.2 = "E3",
                      min.pct = -Inf,
                      logfc.threshold = -Inf,
                      min.diff.pct = -Inf
          )
        saveRDS(DGEA_all_clusters, paste0("../../data/scRNAseq/DGEA_", genotype, ".rds"))
      }
    }

    DGEA_all_clusters <- readRDS(paste0("../../data/scRNAseq/DGEA_", genotype, ".rds"))
    length(DGEA_all_clusters)

    ## ----------------------------
    ## Data wrangling
    ## ----------------------------
    DGEA_all_clusters_ranked <- vector("list", length = length(DGEA_all_clusters))
    names(DGEA_all_clusters_ranked) <- names(DGEA_all_clusters)

    for (cluster in seq_along(names(DGEA_all_clusters))) {

      DGEA_cluster_subset <- DGEA_all_clusters[[cluster]]

      ## Calculate ranking metric for GSEA
      DGEA_cluster_subset$ranking_metric <-
        -log10(DGEA_cluster_subset$p_val) / sign(DGEA_cluster_subset$avg_log2FC)
      DGEA_cluster_subset <-
        DGEA_cluster_subset[order(DGEA_cluster_subset$ranking_metric, decreasing = TRUE), ]

      ## Add Entrez ID
      DGEA_cluster_subset$entrez <- AnnotationDbi::mapIds(org.Mm.eg.db,
                                                          keys = rownames(DGEA_cluster_subset),
                                                          column = "ENTREZID",
                                                          keytype = "ALIAS",
                                                          multiVals = "first")

      geneList_cluster_subset <- DGEA_cluster_subset$ranking_metric
      names(geneList_cluster_subset) <- DGEA_cluster_subset$entrez
      geneList_cluster_subset <- geneList_cluster_subset[!(is.na(geneList_cluster_subset))]
      geneList_cluster_subset <- geneList_cluster_subset[!duplicated(names(geneList_cluster_subset))]

      DGEA_all_clusters_ranked[[cluster]] <- geneList_cluster_subset
    }
    saveRDS(DGEA_all_clusters_ranked, paste0("../../data/scRNAseq/DGEA_geneLists_", genotype, ".rds"))
  }
}

