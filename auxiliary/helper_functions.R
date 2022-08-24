## Define color palettes
pal_sex <- c("#E56C0F", "#0D775F")
pal_APOE <- c("#ca0020", "gray60",  "#0571b0")

## Define contrasts for statistical tests
APOE_groups <- c("APOE2", "APOE3", "APOE4")
APOE_contrasts <- list(c("APOE2", "APOE3"), c("APOE4", "APOE3"), c("APOE2", "APOE4"))
APOE_contrasts_small <- list(c("APOE2", "APOE3"), c("APOE4", "APOE3"))

## Custom aesthetics for forest plots
custom_fm_panels <- list(
  list(width = 0.02),
  list(width = 0.03, display = ~variable, heading = "Variable"),
  list(width = 0.19, display = ~level, fontface = "plain"), 
  list(width = 0.07, display = ~n, hjust = 1, heading = "n"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.24, item = "forest", hjust = 0.5, heading = "HR", linetype = "dashed",
    line_x = 0
  ),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.2, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(
    width = 0.35,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.03)
)

## Define custom GGplot theme
custom_linewidth <-  5/8 * 72.27 / 96 * 0.5
theme_custom2 <-
  ggplot2::theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = custom_linewidth),
    axis.ticks = element_line(size = custom_linewidth),
    axis.ticks.length = unit(0.075, "cm"),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5, colour = "black"),
    plot.title = element_text(size = 7, hjust = 0.5),
    panel.background = element_blank(),
    legend.text = element_text(size = 5),
    legend.title = element_blank(),
    legend.key.size = unit(0.25, "line")
  )

## Plotting functions
IF_plot <- function(df, x, y, fill, 
                    ylab = "ylab", title = "Plot title", label.y = NULL, 
                    expand_limits = NULL) {
  ggplot(df, aes_string(x = x, y = y, fill = fill)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.75,
                                 alpha = 0.6, 
                                 shape = 21,
                                 color = "white", 
                                 fill = "black", 
                                 stroke = 0.2) +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = APOE_contrasts, 
                       size = 1.75, 
                       label.y = label.y) +
    scale_fill_manual(values = pal_APOE) +
    expand_limits(y = expand_limits) +
    labs(title = title, x = "Genotype", y = ylab) +
    guides(x = guide_axis(angle = 60)) +
    theme_custom2 +
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "italic")) +
    lemon::coord_capped_cart(left = "both", bottom = "both")
}