custom_fm_panels <- list(
  list(width = 0.02),
  list(width = 0.03, display = ~variable, heading = "Variable"),
  list(width = 0.19, display = ~level, fontface = "plain"), #ifelse(grepl("APOE", ~level), "italic", "plain")),
  list(width = 0.07, display = ~n, hjust = 1, heading = "n"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.24, item = "forest", hjust = 0.5, heading = "HR", linetype = "dotted",
    line_x = 0
  ),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.2, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.3f (%0.3f, %0.3f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(
    width = 0.35,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.03)
)