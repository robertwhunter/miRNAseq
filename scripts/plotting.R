# pca_plot() ----

library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggridges)

library_plot <- function(cpm_plot) {
  
  cpm_plot %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") -> df_plot
  
  df_plot %>% 
    pivot_longer(cols = -gene, names_to = "library", values_to = "count") %>% 
    mutate(count_log2 = log2(count)) -> df_plotL
  
  df_plotL %>% 
    ggplot(aes(x = library, y = count_log2)) +
    geom_boxplot() + 
    ylab("read cpms \n(log2 transformed)") +
    theme_miR_horizontal() +
    theme(axis.text.x = element_text(angle = 90))
  
}

library_plot_density <- function(cpm_raw, cpm_filtered) {
  
  cpm_raw %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") -> df_plot_raw
  
  df_plot_raw %>% 
    pivot_longer(cols = -gene, names_to = "library", values_to = "count") %>% 
    mutate(count_log2 = log2(count+1), state = "raw") -> df_plot_raw_L
  
  cpm_filtered %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") -> df_plot_filtered
  
  df_plot_filtered %>% 
    pivot_longer(cols = -gene, names_to = "library", values_to = "count") %>% 
    mutate(count_log2 = log2(count+1), state = "filtered") -> df_plot_filtered_L
  
  df_plot <- rbind(df_plot_raw_L, df_plot_filtered_L)
  
  df_plot$state <- df_plot$state %>% as.factor() %>% fct_relevel("raw")
  
  df_plot %>% 
    ggplot(aes(x = count_log2, y = 1, group = library)) +
    geom_density_ridges(alpha = 0.1, fill = NA, color = "#0000001A") +
    xlab("read cpms \n(log2 transformed)") +
    ylab("density") +
    facet_wrap(~state) +
    theme_miR_horizontal()
  
}


pca_plot_RWH <- function(dgl, labels, color, PC1, PC2) { 
  
  ## amended to extract NORMALISED counts from dgl
  
  ## PCA 
  normalised_counts <- cpm(dgl)                                                   # extract normalised counts from dgl
  pcavals <- log2(normalised_counts + 1)                                          # log transform (add 1 to avoid zeros)
  pcavals_var <- pcavals[apply(pcavals, 1, function(x) length(unique(x))) > 1, ]  # remove genes that have the same expression in every sample
  pcavals_m <- pcavals %>% t() %>% as.matrix(scale = T)
  
  pca <- prcomp(pcavals_m)                                                        # do pca
  fraction_explained <- round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100              # get raction of variance explained by each pc 
  
  ## pca plot
  guides <- guides(shape=NULL)
  if(PC2=="PCA2"){guides <- guides(shape=NULL,color=NULL)}
  
  plotdata<-as.data.frame(pca$x)
  plotdata$label = dgl$samples[[labels]]  
  plotdata$color <- as.factor(dgl$samples[[color]])     
  # plotdata$shape <- as.factor(dgl$samples[[shape]])  
  
  ggplot(plotdata, aes_string(PC1, PC2, color = 'color', label = 'label')) + # shape = 'shape'
    geom_point(size = 3) + 
    geom_text_repel(size = 3, segment.size = 0.1, nudge_y = 0.05, nudge_x = 0.05, point.padding = unit(0.5, 'lines'), force = 2, max.iter = 10000) +
    theme_miR() -> pca_plot
  
  
  ## scree plot
  variance <- pca$sdev^2                          # get the variance for each pc
  variance <- variance[1:(length(variance)-1)]    # miss the last one because it's 0 and makes the plot look weird
  pev<-(cumsum(variance)/sum(variance))*100       # calculate the culmative variance explained by increasing numbers of principle components
  plotdata <- data.frame(
    PC = c(1:length(pev)),
    pev = pev,
    type = 'scree'
  )
  
  origin_line <- data.frame(PC = c(0, 1), pev = c(0, pev[1]), type = 'orign')
  
  ggplot(rbind(plotdata, origin_line), aes(PC, pev,color = type, linetype = type)) +
    geom_point(size = 3) + 
    geom_path() + 
    scale_color_manual(values = c('black', 'grey')) +
    ylab('cumulative variance \nexplained (%)') + 
    xlab('principal component') + 
    theme_miR() + theme(legend.position="none") -> scree_plot
  
  list(
    pca = pca,
    pca_plot = pca_plot,
    scree_plot  = scree_plot
  ) %>% return()
  
}


## volcano plot ----

df_prep_volcano <- function(df_results, fc_threshold = 2, fdr_threshold = 0.05) {

  df_results %>% 
    add_rownames("gene") %>% 
    mutate(
      change = case_when(
        # abs(logFC) > log2(fc_threshold) & PValue < fdr_threshold & FDR >= fdr_threshold ~ "p sig (but not fdr)",
        # abs(logFC) > log2(fc_threshold) & FDR < fdr_threshold ~ "fdr sig",
        # abs(logFC) <= log2(fc_threshold) | PValue >= fdr_threshold ~ "not significant"
        logFC > log2(fc_threshold) & FDR < fdr_threshold ~ "up",
        logFC < -log2(fc_threshold) & FDR < fdr_threshold ~ "down",
        abs(logFC) <= log2(fc_threshold) | FDR >= fdr_threshold ~ "not significant"
      ),
      label = case_when(
  #      change %in% c("p sig (but not fdr)", "fdr sig") ~ gene,
        change %in% c("up", "down") ~ gene,
        change == "not significant" ~ "NA"
      )) -> df_dea_toplot
  
  df_dea_toplot$label <- df_dea_toplot$label %>% na_if("NA")
  df_dea_toplot %>% return()
  
}

volcano_plot_RWH <- function(df_results){
  
  df_results %>% 
    ggplot(aes(logFC, -log10(PValue), color=change )) + 
    geom_point(alpha = 0.5) + 
    geom_text_repel(aes(label = label), show.legend = FALSE) +
    # geom_vline(xintercept = log2(fc_threshold)) +
    # geom_vline(xintercept = -log2(fc_threshold)) +
    # geom_hline(yintercept = -log10(0.05)) +
    xlab("log fold-change") +
    ylab("-log<sub>10</sub>(p-value)") +
    theme_miR() +
    theme(
      axis.title.y = element_markdown(),
      plot.title = element_markdown())
}


## gene plot ----

plot_gene <- function(gene, dgl){
  
  expr <- cpm(dgl)
  
  df_plot <- cbind(dgl$samples, expression = expr[gene,])
# plot_data <- plot_data[order(plot_data$group),] 
  df_plot$sample <- factor(rownames(df_plot), levels=rownames(df_plot))
  
  df_plot %>% 
  ggplot(aes(x=sample, y=expression, fill=group)) + 
    geom_bar(stat = "identity") +
    xlab("library") +
    ylab("expression (cpm)") +
    ggtitle(gene) +
    theme_miR_horizontal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) %>% 
    return()

}


plot_gene_dotplot_RWH <- function(gene, dgl){
  
  expr <- cpm(dgl)
  
  df_plot <- cbind(dgl$samples, expression = expr[gene,])
  # plot_data <- plot_data[order(plot_data$group),] 
  df_plot$sample <- factor(rownames(df_plot), levels=rownames(df_plot))
  
  df_plot %>% 
    ggplot(aes(x=group, y=expression)) + 
    geom_boxplot(width = 0.2, fill = NA, alpha = 0.5, colour = "Grey", outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, colour = "Blue") +
    geom_text_repel(aes(label = sample), show.legend = FALSE) +
    xlab("group") +
    ylab("expression (cpm)") +
    ggtitle(gene) +
    theme_miR_vertical() +
    coord_flip()
  # coord_flip() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) %>% 
  # return()
  
}

plot_gene_dotplot_RWH_many <- function(gene_to_plot, dgl){
  
  df_plot <- cpm(dgl) %>% as.data.frame() %>% rownames_to_column("gene")
  df_plot <- df_plot %>% filter(gene %in% gene_to_plot)
  
  df_plot %>% 
    pivot_longer(cols = -gene, names_to = "sample", values_to = "cpm") %>% 
    left_join(dgl$samples)-> df_plot_L
  
  df_plot_L %>% 
    ggplot(aes(x=group, y=cpm)) + 
    geom_boxplot(width = 0.2, fill = NA, alpha = 0.5, colour = "Grey", outlier.shape = NA) +
    geom_jitter(width = 0.05, alpha = 0.3, colour = "Blue") +
#    geom_text_repel(aes(label = sample), show.legend = FALSE) +
    xlab("group") +
    ylab("expression (cpm)") +
    facet_wrap(~gene, scales = "free", ncol = 2) +
    theme_miR_vertical() +
    theme(strip.text = element_text(size = rel(0.5)))
#    coord_flip()

}

# pca_plot_archived <- function(dgl, labels, color, PC1, PC2) { # shape
#   
#   ## PCA 
#   pcavals <- log2(dgl$counts + 1)                                                 # log transform (add 1 to avoid zeros)
#   pcavals_var <- pcavals[apply(pcavals, 1, function(x) length(unique(x))) > 1, ]  # remove genes that have the same expression in every sample
#   pcavals_m <- pcavals %>% t() %>% as.matrix(scale = T)
#   
#   pca <- prcomp(pcavals_m)                                                        # do pca
#   fraction_explained <- round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100              # get raction of variance explained by each pc 
#   
#   ## pca plot
#   guides <- guides(shape=NULL)
#   if(PC2=="PCA2"){guides <- guides(shape=NULL,color=NULL)}
#   
#   plotdata<-as.data.frame(pca$x)
#   plotdata$label = dgl$samples[[labels]]  
#   plotdata$color <- as.factor(dgl$samples[[color]])     
#   # plotdata$shape <- as.factor(dgl$samples[[shape]])  
#   
#   ggplot(plotdata, aes_string(PC1, PC2, color = 'color', label = 'label')) + # shape = 'shape'
#     geom_point(size = 3) + 
#     geom_text_repel(size = 3, segment.size = 0.1, nudge_y = 0.05, nudge_x = 0.05, point.padding = unit(0.5, 'lines'), force = 2, max.iter = 10000) +
#     theme_miR() -> pca_plot
#   
#   
#   ## scree plot
#   variance <- pca$sdev^2                          # get the variance for each pc
#   variance <- variance[1:(length(variance)-1)]    # miss the last one because it's 0 and makes the plot look weird
#   pev<-(cumsum(variance)/sum(variance))*100       # calculate the culmative variance explained by increasing numbers of principle components
#   plotdata <- data.frame(
#     PC = c(1:length(pev)),
#     pev = pev,
#     type = 'scree'
#   )
#   
#   origin_line <- data.frame(PC = c(0, 1), pev = c(0, pev[1]), type = 'orign')
#   
#   ggplot(rbind(plotdata, origin_line), aes(PC, pev,color = type, linetype = type)) +
#     geom_point(size = 3) + 
#     geom_path() + 
#     scale_color_manual(values = c('black', 'grey')) +
#     ylab('cumulative variance \nexplained (%)') + 
#     xlab('principal component') + 
#     theme_miR() + theme(legend.position="none") -> scree_plot
#   
#   list(
#     pca = pca,
#     pca_plot = pca_plot,
#     scree_plot  = scree_plot
#   ) %>% return()
#   
# }

# plot_gene_dotplot_archived <- function(gene, dgl){
#   
#   expr <- cpm(dgl)
#   
#   df_plot <- cbind(dgl$samples, expression = expr[gene,])
#   # plot_data <- plot_data[order(plot_data$group),] 
#   df_plot$sample <- factor(rownames(df_plot), levels=rownames(df_plot))
#   
#   df_plot %>% 
#     ggplot(aes(x=group, y=expression)) + 
#     geom_boxplot(width = 0.2, fill = NA, alpha = 0.5, colour = "Grey", outlier.shape = NA) +
#     geom_jitter(width = 0.1, alpha = 0.3, colour = "Blue") +
#     xlab("group") +
#     ylab("expression (cpm)") +
#     ggtitle(gene) +
#     theme_miR_vertical() +
#     coord_flip()
#     # coord_flip() +
#     # theme(axis.text.x = element_text(angle = 90, hjust = 1)) %>% 
#     # return()
#   
# }

# volcano_plot_archived <- function(df_results, fc_threshold = 2, fdr_threshold = 0.05, log=TRUE){
#   
#   df_results$significant <- 'no'
#   df_results$significant[ abs(df_results$logFC) >= log2(fc_threshold) &
#                               df_results$FDR <= fdr_threshold  ] <- 'yes'
#   
#   if(! log){
#     df_results$logFC <- sign(df_results$logFC)*2**abs(df_results$logFC) 
#   }
#   
#   df_results %>% 
#     ggplot(aes(logFC, -log10(FDR), color=significant )) + 
#     geom_point(alpha = 0.5) + 
#     scale_colour_manual(name = 'significant', 
#                         values = setNames(c('red','grey'),c('yes', 'no'))) +
#     theme_miR()
# }

