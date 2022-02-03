TSNE_genelist <- function(beta_mat, ncomp = 50) {
  ## Do PCA
  beta_mat_svd <- RSpectra::svds(beta_mat, k = ncomp)
  beta_mat_pcscore <- t(t(beta_mat_svd$u) * beta_mat_svd$d)

  ## Do TSNE
  TSNE <- Rtsne::Rtsne(beta_mat_pcscore, pca = FALSE, check_duplicates = TRUE)$Y

  ## Calculate the Mahalanobis distance
  gene_list <- mahalanobis(TSNE, center = colMeans(TSNE), cov = cov(TSNE))
  names(gene_list) <- toupper(rownames(beta_mat))
  gene_list <- sort(gene_list, decreasing = TRUE)

  ## return results
  res <- list()
  res$TSNE <- TSNE
  res$gene_list <- gene_list
  return(res)
}