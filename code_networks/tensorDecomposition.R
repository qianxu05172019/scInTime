source("code_networks/cpDecomposition.R")
tensorDecomposition <- function(xList, K = 5, maxError = 1e-5, maxIter = 1e3) {
  xNets <- length(xList)
  nNet <- xNets
  xGenes <- unique(unlist(lapply(xList, rownames)))
  sGenes <- xGenes
  nGenes <- length(sGenes)

  tensorX <- array(data = 0, dim = c(nGenes, nGenes, 1, nNet))

  for (i in seq_len(nNet)) {
    tempX <- matrix(0, nGenes, nGenes)
    rownames(tempX) <- colnames(tempX) <- sGenes
    temp <- as.matrix(xList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempX[tGenes, tGenes] <- temp[tGenes, tGenes]
    tensorX[, , , i] <- tempX
  }

  set.seed(1)
  tensorX <- as.tensor(tensorX)
  tensorX <- cpDecomposition(tnsr = tensorX, num_components = K, max_iter = maxIter, tol = maxError)
  tensorOutput <- list()
  for (i in seq_len(nNet)) {
    tensorOutput[[i]] <- tensorX$est$data[, , , i]
  }
  return(tensorOutput)
}