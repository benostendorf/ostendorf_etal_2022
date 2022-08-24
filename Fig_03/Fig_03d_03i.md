Flow cytometry Fig 3, ED Fig 4
================
Benjamin Ostendorf
2022-08-08

## Preamble

``` r
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

source("../auxiliary/helper_functions.R")
knitr::opts_chunk$set(fig.retina = 3)
```

## Fig 3d, ED Fig 4b

``` r
## Import data
df_lungs <- read_tsv("../data/FC_lungs.tsv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   CD11b = col_double(),
    ##   Ly6C = col_double(),
    ##   Ly6G = col_double(),
    ##   CD4 = col_double(),
    ##   CD8 = col_double(),
    ##   B_cells = col_double(),
    ##   NK_cells = col_double(),
    ##   Genotype = col_character()
    ## )

``` r
## Select variables
pops_to_plot <- colnames(df_lungs)[1:7]

## Plot
for (v in pops_to_plot) {
  
  label.y_coordinates <- c(rep(max(as.numeric(df_lungs[[v]])) * 1.05, 2), 
                           max(as.numeric(df_lungs[[v]])) * 1.2)
  
  alternative_t_test <- if  (v %in% c("CD11b", "Ly6C", "Ly6G")) "greater" else "less"
  
  p_variable <-
    df_lungs |>
    ggplot(aes_string(x = "Genotype", fill = "Genotype", y = v)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.75, alpha = 0.6, shape = 21, 
                               color = "white", fill = "black", stroke = 0.2) +
    stat_compare_means(method = "t.test", method.args = list(alternative = alternative_t_test),
                       comparisons = APOE_contrasts_small,
                       label = "p.format", size = 1.75, label.y = label.y_coordinates) +
    scale_fill_manual(values = pal_APOE) +
    guides(x = guide_axis(angle = 60)) +
    expand_limits(y = c(0, max(as.numeric(df_lungs[[v]])) * 1.3)) +
    labs(title = gsub("_", " ", v), 
         y = bquote('Fraction of '~CD45^{textstyle("+")}), 
         x = "Genotype") +
    theme_custom2  +
    theme(legend.position = "none", 
          title = element_text(size = 7), 
          axis.text.x = element_text(face = "italic")) +
    lemon::coord_capped_cart(left = "both", bottom = "both")
  
  print(p_variable)
}
```

<img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-2-1.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-2-2.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-2-3.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-2-4.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-2-5.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-2-6.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-2-7.png" width="672" />

## Fig 3i

``` r
df_tetramer <- read_tsv("../data/FC_tetramer.tsv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   genotype = col_character(),
    ##   tetramer = col_double()
    ## )

``` r
df_tetramer |>
  ggplot(aes(x = genotype, fill = genotype, y = tetramer)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(size = 0.75, alpha = 0.6, shape = 21, 
                             color = "white", fill = "black", stroke = 0.2) +
  stat_compare_means(method = "t.test", method.args = list(var.equal = FALSE), 
                     comparisons = APOE_contrasts,
                     label = "p.format", size = 1.75, label.y = c(20, 20, 22)) +
  scale_fill_manual(values = pal_APOE) +
  guides(x = guide_axis(angle = 60)) +
  expand_limits(y = c(0, max(as.numeric(df_tetramer$tetramer)) * 1.4)) +
  labs(title = "Tetramer+", y = "Fraction of CD8+ (%)", x = "Genotype") +
  theme_custom2  +
  theme(legend.position = "none", 
        title = element_text(size = 7), 
        axis.text.x = element_text(face = "italic")) +
  lemon::coord_capped_cart(left = "both", bottom = "both")
```

<img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-3-1.png" width="672" />

## ED Fig 4d-f

``` r
## ---------------------------------------------
## Import data
## ---------------------------------------------
df_PB_abs <- read_tsv("../data/FC_PB_abs.tsv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   genotype = col_character(),
    ##   CD45 = col_double()
    ## )

``` r
df_PB <- read_tsv("../data/FC_PB.tsv")
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   genotype = col_character(),
    ##   CD11b = col_double(),
    ##   Ly6C = col_double(),
    ##   Ly6G = col_double(),
    ##   CD4 = col_double(),
    ##   CD8 = col_double(),
    ##   B_cells = col_double(),
    ##   NK_cells = col_double()
    ## )

``` r
## ---------------------------------------------
## ED Fig 4d
## ---------------------------------------------
df_PB_abs |>
  ggplot(aes(x = genotype, fill = genotype, y = CD45)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(size = 0.75, alpha = 0.6, shape = 21, 
                             color = "white", fill = "black", stroke = 0.2) +
  stat_compare_means(method = "t.test", comparisons = APOE_contrasts,
                     label = "p.format", size = 1.75, label.y = c(3600, 3600, 4000)) +
  scale_fill_manual(values = pal_APOE) +
  guides(x = guide_axis(angle = 60)) +
  expand_limits(y = c(0, 4000)) +
  labs(title = "CD45+", 
       y = "Cells per µl", 
       x = "Genotype") +
  theme_custom2  +
  theme(legend.position = "none", 
        title = element_text(size = 7), 
        axis.text.x = element_text(face = "italic")) +
  lemon::coord_capped_cart(left = "both", bottom = "both")
```

<img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-4-1.png" width="672" />

``` r
## ---------------------------------------------------------------
## ED Fig 4e-f
## ---------------------------------------------------------------
cell_populations <- colnames(df_PB)[2:8]

for (v in cell_populations) {
  
  label.y_coordinates <- c(rep(max(as.numeric(df_PB[[v]])) * 1.05, 2), 
                           max(as.numeric(df_PB[[v]])) * 1.15)
  
  p_variable <-
    df_PB |>
    ggplot(aes_string(x = "genotype", fill = "genotype", y = v)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.75, alpha = 0.6, shape = 21, 
                               color = "white", fill = "black", stroke = 0.2) +
    stat_compare_means(method = "t.test", comparisons = APOE_contrasts,
                       label = "p.format", size = 1.75, label.y = label.y_coordinates) +
    scale_fill_manual(values = pal_APOE) +
    guides(x = guide_axis(angle = 60)) +
    expand_limits(y = c(0)) +
    labs(title = gsub("_", " ", v), 
         y = "Fraction of CD45+ (%)", 
         x = "Genotype") +
    theme_custom2  +
    theme(legend.position = "none", 
          title = element_text(size = 7), 
          axis.text.x = element_text(face = "italic")) +
    lemon::coord_capped_cart(left = "both", bottom = "both")
  print(p_variable)
}
```

<img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-4-2.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-4-3.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-4-4.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-4-5.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-4-6.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-4-7.png" width="672" /><img src="Fig_03d_03i_files/figure-gfm/unnamed-chunk-4-8.png" width="672" />

## Session info

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 4.1.0 (2021-05-18)
    ##  os       macOS Big Sur 10.16         
    ##  system   x86_64, darwin17.0          
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  ctype    en_US.UTF-8                 
    ##  tz       Europe/Berlin               
    ##  date     2022-08-24                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package      * version date       lib source        
    ##  abind          1.4-5   2016-07-21 [1] CRAN (R 4.1.0)
    ##  assertthat     0.2.1   2019-03-21 [1] CRAN (R 4.1.0)
    ##  backports      1.2.1   2020-12-09 [1] CRAN (R 4.1.0)
    ##  beeswarm       0.3.1   2021-03-07 [1] CRAN (R 4.1.0)
    ##  broom          0.7.6   2021-04-05 [1] CRAN (R 4.1.0)
    ##  cachem         1.0.5   2021-05-15 [1] CRAN (R 4.1.0)
    ##  callr          3.7.0   2021-04-20 [1] CRAN (R 4.1.0)
    ##  car            3.0-10  2020-09-29 [1] CRAN (R 4.1.0)
    ##  carData        3.0-4   2020-05-22 [1] CRAN (R 4.1.0)
    ##  cellranger     1.1.0   2016-07-27 [1] CRAN (R 4.1.0)
    ##  cli            2.5.0   2021-04-26 [1] CRAN (R 4.1.0)
    ##  colorspace     2.0-1   2021-05-04 [1] CRAN (R 4.1.0)
    ##  crayon         1.4.1   2021-02-08 [1] CRAN (R 4.1.0)
    ##  curl           4.3.1   2021-04-30 [1] CRAN (R 4.1.0)
    ##  data.table     1.14.0  2021-02-21 [1] CRAN (R 4.1.0)
    ##  DBI            1.1.1   2021-01-15 [1] CRAN (R 4.1.0)
    ##  dbplyr         2.1.1   2021-04-06 [1] CRAN (R 4.1.0)
    ##  desc           1.3.0   2021-03-05 [1] CRAN (R 4.1.0)
    ##  devtools       2.4.1   2021-05-05 [1] CRAN (R 4.1.0)
    ##  digest         0.6.29  2021-12-01 [1] CRAN (R 4.1.0)
    ##  dplyr        * 1.0.6   2021-05-05 [1] CRAN (R 4.1.0)
    ##  ellipsis       0.3.2   2021-04-29 [1] CRAN (R 4.1.0)
    ##  evaluate       0.14    2019-05-28 [1] CRAN (R 4.1.0)
    ##  fansi          0.5.0   2021-05-25 [1] CRAN (R 4.1.0)
    ##  farver         2.1.0   2021-02-28 [1] CRAN (R 4.1.0)
    ##  fastmap        1.1.0   2021-01-25 [1] CRAN (R 4.1.0)
    ##  forcats      * 0.5.1   2021-01-27 [1] CRAN (R 4.1.0)
    ##  foreign        0.8-81  2020-12-22 [1] CRAN (R 4.1.0)
    ##  fs             1.5.0   2020-07-31 [1] CRAN (R 4.1.0)
    ##  generics       0.1.0   2020-10-31 [1] CRAN (R 4.1.0)
    ##  ggbeeswarm     0.6.0   2017-08-07 [1] CRAN (R 4.1.0)
    ##  ggplot2      * 3.3.5   2021-06-25 [1] CRAN (R 4.1.0)
    ##  ggpubr       * 0.4.0   2020-06-27 [1] CRAN (R 4.1.0)
    ##  ggsignif       0.6.1   2021-02-23 [1] CRAN (R 4.1.0)
    ##  glue           1.6.0   2021-12-17 [1] CRAN (R 4.1.0)
    ##  gridExtra      2.3     2017-09-09 [1] CRAN (R 4.1.0)
    ##  gtable         0.3.0   2019-03-25 [1] CRAN (R 4.1.0)
    ##  haven          2.4.1   2021-04-23 [1] CRAN (R 4.1.0)
    ##  highr          0.9     2021-04-16 [1] CRAN (R 4.1.0)
    ##  hms            1.1.0   2021-05-17 [1] CRAN (R 4.1.0)
    ##  htmltools      0.5.2   2021-08-25 [1] CRAN (R 4.1.0)
    ##  httr           1.4.2   2020-07-20 [1] CRAN (R 4.1.0)
    ##  jsonlite       1.7.2   2020-12-09 [1] CRAN (R 4.1.0)
    ##  knitr          1.37    2021-12-16 [1] CRAN (R 4.1.0)
    ##  labeling       0.4.2   2020-10-20 [1] CRAN (R 4.1.0)
    ##  lattice        0.20-44 2021-05-02 [1] CRAN (R 4.1.0)
    ##  lemon          0.4.5   2020-06-08 [1] CRAN (R 4.1.0)
    ##  lifecycle      1.0.0   2021-02-15 [1] CRAN (R 4.1.0)
    ##  lubridate      1.7.10  2021-02-26 [1] CRAN (R 4.1.0)
    ##  magrittr       2.0.1   2020-11-17 [1] CRAN (R 4.1.0)
    ##  memoise        2.0.0   2021-01-26 [1] CRAN (R 4.1.0)
    ##  modelr         0.1.8   2020-05-19 [1] CRAN (R 4.1.0)
    ##  munsell        0.5.0   2018-06-12 [1] CRAN (R 4.1.0)
    ##  openxlsx       4.2.3   2020-10-27 [1] CRAN (R 4.1.0)
    ##  patchwork    * 1.1.1   2020-12-17 [1] CRAN (R 4.1.0)
    ##  pillar         1.6.1   2021-05-16 [1] CRAN (R 4.1.0)
    ##  pkgbuild       1.2.0   2020-12-15 [1] CRAN (R 4.1.0)
    ##  pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.1.0)
    ##  pkgload        1.2.1   2021-04-06 [1] CRAN (R 4.1.0)
    ##  plyr           1.8.6   2020-03-03 [1] CRAN (R 4.1.0)
    ##  prettyunits    1.1.1   2020-01-24 [1] CRAN (R 4.1.0)
    ##  processx       3.5.2   2021-04-30 [1] CRAN (R 4.1.0)
    ##  ps             1.6.0   2021-02-28 [1] CRAN (R 4.1.0)
    ##  purrr        * 0.3.4   2020-04-17 [1] CRAN (R 4.1.0)
    ##  R6             2.5.1   2021-08-19 [1] CRAN (R 4.1.0)
    ##  RColorBrewer * 1.1-2   2014-12-07 [1] CRAN (R 4.1.0)
    ##  Rcpp           1.0.7   2021-07-07 [1] CRAN (R 4.1.0)
    ##  readr        * 1.4.0   2020-10-05 [1] CRAN (R 4.1.0)
    ##  readxl         1.3.1   2019-03-13 [1] CRAN (R 4.1.0)
    ##  remotes        2.3.0   2021-04-01 [1] CRAN (R 4.1.0)
    ##  reprex         2.0.0   2021-04-02 [1] CRAN (R 4.1.0)
    ##  rio            0.5.26  2021-03-01 [1] CRAN (R 4.1.0)
    ##  rlang          0.4.12  2021-10-18 [1] CRAN (R 4.1.0)
    ##  rmarkdown      2.11    2021-09-14 [1] CRAN (R 4.1.0)
    ##  rprojroot      2.0.2   2020-11-15 [1] CRAN (R 4.1.0)
    ##  rstatix        0.7.0   2021-02-13 [1] CRAN (R 4.1.0)
    ##  rstudioapi     0.13    2020-11-12 [1] CRAN (R 4.1.0)
    ##  rvest          1.0.0   2021-03-09 [1] CRAN (R 4.1.0)
    ##  scales         1.1.1   2020-05-11 [1] CRAN (R 4.1.0)
    ##  sessioninfo    1.1.1   2018-11-05 [1] CRAN (R 4.1.0)
    ##  stringi        1.7.6   2021-11-29 [1] CRAN (R 4.1.0)
    ##  stringr      * 1.4.0   2019-02-10 [1] CRAN (R 4.1.0)
    ##  testthat       3.0.2   2021-02-14 [1] CRAN (R 4.1.0)
    ##  tibble       * 3.1.2   2021-05-16 [1] CRAN (R 4.1.0)
    ##  tidyr        * 1.1.3   2021-03-03 [1] CRAN (R 4.1.0)
    ##  tidyselect     1.1.1   2021-04-30 [1] CRAN (R 4.1.0)
    ##  tidyverse    * 1.3.1   2021-04-15 [1] CRAN (R 4.1.0)
    ##  usethis        2.0.1   2021-02-10 [1] CRAN (R 4.1.0)
    ##  utf8           1.2.1   2021-03-12 [1] CRAN (R 4.1.0)
    ##  vctrs          0.3.8   2021-04-29 [1] CRAN (R 4.1.0)
    ##  vipor          0.4.5   2017-03-22 [1] CRAN (R 4.1.0)
    ##  withr          2.4.2   2021-04-18 [1] CRAN (R 4.1.0)
    ##  xfun           0.29    2021-12-14 [1] CRAN (R 4.1.0)
    ##  xml2           1.3.2   2020-04-23 [1] CRAN (R 4.1.0)
    ##  yaml           2.2.1   2020-02-01 [1] CRAN (R 4.1.0)
    ##  zip            2.1.1   2020-08-27 [1] CRAN (R 4.1.0)
    ## 
    ## [1] /Library/Frameworks/R.framework/Versions/4.1/Resources/library
