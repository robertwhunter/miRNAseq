#### NOTES ----

# Input dataframe must have rows = variables (e.g. gene expression) & columns = samples (e.g. sequencing libraries)
# This is then transposed in the functions below
#
# Meta-data needs to have "sample_name" and "group" columns

# See https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/

# When scaling, transform as [xi - mean(x)] / sd(x) - so will not work if sd(x) is zero...
# ...so first filter these data out


#### setup ----
library(ggdendro)
library(broom)

#### clustering functions ----
reshape_for_clustering <- function(df) df %>% select(-1) %>% as.matrix %>% t() %>% scale()
remove_rows_with_stdev_zero <- function(df) {
  df %>% 
    rowwise() %>% 
    mutate(stdev_row = sd(c_across(cols = -1))) %>% 
    dplyr::filter(stdev_row != 0) %>% select(-stdev_row)
  }

dendrogram_RWH <- function(df, df_meta = NULL) {
  df %>% 
    reshape_for_clustering() %>% 
    dist() %>% 
    hclust() %>% 
    dendro_data() -> df_dendro
  
    df_dendro %>% 
    ggdendrogram() -> p
  
  if (is.null(df_meta)) {
    p %>% return()
    } else {
    p +
    geom_point(
      data = df_meta, 
      aes(x = match(df_meta$sample_name, df_dendro$label$label), y = -0.7, fill = group),
      size = 5,
      shape = 21) %>% return()
  }
}

pca_plot_RWH_prcomp <- function(df, df_meta = NULL) {
  df %>% 
    remove_rows_with_stdev_zero() %>% 
    reshape_for_clustering() %>% 
    prcomp(scale. = FALSE) %>% 
    augment() -> df_pca
    
    if (is.null(df_meta)) {
      df_pca %>%
        ggplot(aes(x = .fittedPC1, y = .fittedPC2, label = .rownames)) +
        geom_point(alpha = 1, size = 3) +
        geom_text_repel(size = 3, segment.size = 0.1, nudge_y = 0.05, nudge_x = 0.05, point.padding = unit(0.5, 'lines'), force = 2, max.iter = 10000) +
        xlab("PC1") +
        ylab("PC2") +
        theme_miR()
    } else {
      df_pca <- df_pca %>% left_join(df_meta, by = c(".rownames" = "sample_name"))
      df_pca %>%
        ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = group, label = .rownames)) +
        geom_point(alpha = 1, size = 3) +
        geom_text_repel(size = 3, segment.size = 0.1, nudge_y = 0.05, nudge_x = 0.05, point.padding = unit(0.5, 'lines'), force = 2, max.iter = 10000) +
        xlab("PC1") +
        ylab("PC2") +
        theme_miR() 
    }
}
