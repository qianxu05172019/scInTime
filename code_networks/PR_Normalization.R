PR_Normalization <- function(X) {
  X_sum <- sum(X)
  X_rowsum <- rowSums(X)
  X_colSum <- colSums(X)
  u_hat <- X_rowsum %*% t(X_colSum) / X_sum
  Z <- (X - u_hat) / sqrt(u_hat + u_hat^2 / 100)
  return(Z)
}