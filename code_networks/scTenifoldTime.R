source("code_networks/tensorDecomposition.R")
source("code_networks/pcNet.R")
source("code_networks/regression_coefficienets.R")
source("code_networks/TSNE_genelist.R")
scTenifoldTime <- function(dta_list, method = "pcnet", time_vec, nComp = 5, q = 0,
                           K = 10, maxIter = 10000, maxError = 1e-5) {
  res <- list()
  n_net <- length(dta_list)
  nGenes <- nrow(dta_list[[1]])
  gene_name <- rownames(dta_list[[1]])
  ## generate network
  network_list <- list()
  for (i in seq_len(n_net)) {
    if (method == "pcnet") {
      network_list[[i]] <- pcNet(dta_list[[i]], nComp = nComp, scaleScores = TRUE, symmetric = FALSE, q = q, verbose = TRUE)
      network_list[[i]] <- round(network_list[[i]], 3)
      diag(network_list[[i]]) <- 0
      network_list[[i]] <- as(network_list[[i]], "dgCMatrix")
    } else {
      network_list[[i]] <- cor(t(dta_list[[i]]), method = "spearman")
      network_list[[i]] <- round(network_list[[i]], 3)
      diag(network_list[[i]]) <- 0
      network_list[[i]] <- as(network_list[[i]], "dgCMatrix")
    }
    rownames(network_list[[i]]) <- colnames(network_list[[i]]) <- gene_name
    print(paste0("Finish network ", i))
  }
  res$network_list <- network_list
  ## tensor decomposition
  set.seed(1)
  network_tensor <- tensorDecomposition(network_list, K = K, maxIter = maxIter, maxError = maxError)
  for (i in seq_len(length(network_tensor))) {
    network_tensor[[i]] <- round(network_tensor[[i]], 3)
    diag(network_tensor[[i]]) <- 0
    network_tensor[[i]] <- as(network_tensor[[i]], "dgCMatrix")
  }
  res$network_tensor <- network_tensor
  print("Finish tensor decomposition part.")

  ## Do regression without tensor decomposition
  beta_mat <- regression_coefficienets(network_list = network_list, time_vec = time_vec)
  diag(beta_mat) <- 0
  rownames(beta_mat) <- colnames(beta_mat) <- gene_name
  res$beta_mat <- beta_mat

  ## Get low dimensional representation and rank gene list
  set.seed(1)
  res$gene_list <- TSNE_genelist(beta_mat)


  ## Do regression with tensor decomposition
  beta_mat_tensor <- regression_coefficienets(network_list = network_tensor, time_vec = time_vec)
  diag(beta_mat_tensor) <- 0
  rownames(beta_mat_tensor) <- colnames(beta_mat_tensor) <- gene_name
  res$beta_mat_tensor <- beta_mat_tensor

  ## Get low dimensional representation and rank gene list
  set.seed(1)
  res$gene_list_tensor <- TSNE_genelist(beta_mat_tensor)

  ## return results
  return(res)
}