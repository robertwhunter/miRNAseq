# pca_plot() ----

library(ggplot2)
library(ggrepel)

pca_plot <- function(dgl, labels, color, PC1, PC2) { # shape
  
  ## PCA 
  pcavals <- log2(dgl$counts + 1)                                                 # log transform (add 1 to avoid zeros)
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

volcano_plot <- function(df_results, fc_threshold = 2, fdr_threshold = 0.05, log=TRUE){
  
  df_results$significant <- 'no'
  df_results$significant[ abs(df_results$logFC) >= log2(fc_threshold) &
                              df_results$FDR <= fdr_threshold  ] <- 'yes'
  
  if(! log){
    df_results$logFC <- sign(df_results$logFC)*2**abs(df_results$logFC) 
  }
  
  df_results %>% 
    ggplot(aes(logFC, -log10(FDR), color=significant )) + 
    geom_point(alpha = 0.5) + 
    scale_colour_manual(name = 'significant', 
                        values = setNames(c('red','grey'),c('yes', 'no'))) +
    theme_miR()
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

