regression_coefficienets <- function(network_list, time_vec) {
  ## Getting beta mat
  n_net <- length(network_list)
  X_mat <- cbind(1, time_vec)
  tmp <- as.numeric(network_list[[1]])
  nGenes <- nrow(network_list[[1]])
  Y <- tmp
  if (n_net > 1) {
    for (i in 2:n_net) {
      tmp <- as.numeric(network_list[[i]])
      Y <- rbind(Y, tmp)
    }
  }
  ## Do regress
  beta_coef <- solve(crossprod(X_mat), crossprod(X_mat, Y))
  beta_mat <- matrix(beta_coef[2, ], nrow = nGenes)
  # return results
  return(beta_mat)
}