APOE_contrasts_E3 <- list(c("APOE2", "APOE3"), c("APOE4", "APOE3"))

## --------------------------------------
## Plot cell proportions
## --------------------------------------
plot_cell_prop <- function(df, x, y, fill, facet_var, 
                           comparisons, x_axis_face = "italic", 
                           ylab = "Cell fraction (%)", 
                           fill_pal = pal_APOE, nrow = NULL){
  
  facet_formula <- as.formula(paste0("~ ", facet_var))
  
  ## Create blank data frame to force individual y-limits
  blank_data <- data.frame(celltype = levels(df[["celltype"]]),
                           condition = "inf",
                           genotype = "APOE2", 
                           y_limit_indiv = df |>
                             group_by(celltype) |>
                             summarize(y_max = max(Freq) * 1.35) |> pull(y_max)) |>
    mutate(celltype = ordered(celltype, levels = names(table(df_filt$celltype_grouped)))) |>
    mutate(celltype = droplevels(celltype))
  
  ggplot(df, aes(x = {{x}}, y = {{y}}, fill = {{fill}})) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.75, alpha = 0.6, shape = 21, 
                                 color = "white", fill = "black", stroke = 0.2) +
    stat_compare_means(method = "t.test", comparisons = comparisons,
                       label = "p.format", size = 1.75) + 
    geom_blank(data = blank_data, aes(x = {{x}}, y = y_limit_indiv)) +
    scale_fill_manual(values = fill_pal) +
    expand_limits(y = 0) +
    guides(x = guide_axis(angle = 60)) +
    labs(y = ylab) +
    theme_custom2 +
    theme(legend.position = c(0.9, 0.1), 
          axis.text.x = element_text(face = x_axis_face),
          axis.title.x = element_blank(), 
          strip.text = element_text(size = 5, margin = margin(0.5,0,0.5,0, "mm")))  +
    facet_wrap(facet_formula, scales = "free_y", nrow = nrow) +
    lemon::coord_capped_cart(left = "both", bottom = "both")
}

## --------------------------------------
## Annotate subset proportion dfs
## --------------------------------------
annotate_prop_df <-
  function(df){
    
    df_mut <- 
      df |>
      as_tibble() |>
      dplyr::rename(celltype = Var1, sample = Var2) |>
      mutate(condition = case_when(grepl("ctrl", sample) ~ "ctrl", 
                                   grepl("inf", sample) ~ "inf"), 
             genotype = case_when(grepl("E2", sample) ~ "APOE2", 
                                  grepl("E3", sample) ~ "APOE3", 
                                  grepl("E4", sample) ~ "APOE4"), 
             celltype = ordered(celltype, levels = names(table(df_filt$celltype_grouped)))) |>
      mutate(celltype = droplevels(celltype))
    return(df_mut)
  }

## --------------------------------------
## Density plot
## --------------------------------------
plot_density <- 
  function(UMAP_df, facet_var){
    
    facet_formula <- as.formula(paste0(". ~ ", facet_var))
    
    density_plot <- 
      ggplot(data = UMAP_df, aes(x = UMAP_1, y = UMAP_2)) +
      stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 1000) + 
      ggrastr::geom_point_rast(color = "white", size = 0.02, 
                               stroke = 0, shape = 16, raster.dpi = 600) +
      scale_colour_viridis_c(alpha = 1, option = "A", aesthetics = "fill") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "black", size = 0.01), 
            panel.grid = element_blank(), 
            axis.title = element_text(size = 6), 
            axis.text = element_blank(), 
            plot.title = element_text(size = 6), 
            legend.title = element_text(size = 6), 
            legend.text = element_text(size = 6), 
            legend.key.size = unit(0.5, "line"), 
            strip.text.x = element_text(size = 7)) +
      xlab("UMAP dim. 1") +
      ylab("UMAP dim. 2") +
      facet_grid(facet_formula)
    return(density_plot)
  }

## --------------------------------------
## Canonical markers as combination of DEG results and Angelidis et al, Nature Communications, 2019
## --------------------------------------
lineage_markers_lung <- 
  list("alveolar_macs" = c("Ear2", "Mrc1"), 
       "interstitial_macs" = c("Adgre1", "C1qa"), 
       "monocytes" = c("Vcan", "Plac8"), 
       "granulocytes" = c("Ly6g", "Retnlg"), 
       "DCs" = c("Flt3", "Xcr1"),
       "NK_cells" = c("Ncr1", "Klra7"), 
       "T_cells" = c("Cd3g", "Tox", "Cd8a", "Cd4", "Foxp3"), 
       "proliferation" = c("Top2a", "Mki67"), 
       "B_cells" = c("Bank1", "Cd79b"), 
       "myofibroblasts" = c("Pdgfra", "Col3a1", "Eln"), 
       "other_fibroblasts" = c("Gpx3", "Inmt", "Col14a1"), 
       "capillary_ECs" = c("Kdr", "Ednrb"), 
       "vascular_ECs" = c("Gpihbp1", "Hpgd", "Tspan7"), 
       "vcam1_ECs" = c("Vcam1", "Vwf"), 
       "pericytes" = c("Trpc6", "Pdgfrb"), 
       "AT1" = c("Rtkn2", "Col4a3"), 
       "AT2" = c("Sftpc", "Sftpa1"), 
       "ciliated_cells" = c("Ccdc153", "Dynlrb2"), 
       "airway_epithelial" = c("Chst9", "Ntm"), 
       "mesothelial_cells" = c("Igfbp5", "Rarres2"), 
       "neuronal" = c("Reln", "Sema3a"))

## --------------------------------------
## Custom ggplot theme
## --------------------------------------
theme_custom <- theme(axis.line = element_line(size = custom_linewidth), 
                      axis.text = element_text(size = 6),
                      axis.title = element_text(size = 6),
                      axis.title.x = element_blank(),
                      axis.ticks = element_line(size = custom_linewidth), 
                      axis.ticks.length = unit(0.075, "cm"),
                      plot.title = element_text(size = 7, face = "plain", hjust = 0.5), 
                      legend.title = element_blank(),
                      legend.key = element_rect(size = 0), 
                      legend.text = element_text(size = 6),
                      legend.key.size = unit(0.4, "line")
)

## --------------------------------------
## Custom ComplexHeatmap options
## --------------------------------------
ht_opt(heatmap_column_title_gp = gpar(fontsize = 6),
       heatmap_border = TRUE,
       annotation_border = TRUE, 
       legend_labels_gp = gpar(fontsize = 6),
       legend_title_gp = gpar(fontsize = 6))
