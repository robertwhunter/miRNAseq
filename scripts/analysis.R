## pca_pvals ----

pca_pvals <- function(pca, variables, dgl) {

  pcameta <- dgl$samples[,variables]  # variables to test for association with pcs 
  pcameta <- data.frame(pcameta)       # RWH added this line 
  
  last_pc <- 10                       # number of pcs to include
  
  ## create an empty matrix (cols = variables; rows = pc)
  pvals <- matrix(
    data = NA, 
    ncol = ncol(pcameta), 
    nrow = last_pc, 
    dimnames = list(as.character(1:last_pc), colnames(pcameta))
  )
  
  # fill the matrix with anova p values
  for (i in 1:ncol(pcameta)) {
    for (j in 1:last_pc) {
      fit <- aov(pca$x[, j] ~ as.factor(pcameta[, i]))
      if ("Pr(>F)" %in% names(summary(fit)[[1]])) {
        pvals[j, i] <- summary(fit)[[1]][["Pr(>F)"]][[1]]
      }
    }
  }
  
  # calculate the percent variance explained by each component
  fraction_explained <- round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100
  names(fraction_explained) <- colnames(pca$x)
  rownames(pvals) <- paste(paste('PC', 1:last_pc, sep = ""), " (", fraction_explained[1:last_pc],  "%)", sep = "" )
  
  pvals %>% return()

}