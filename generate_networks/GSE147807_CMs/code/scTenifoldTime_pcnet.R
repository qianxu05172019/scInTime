## load data and packages
source("code_networks/PR_Normalization.R")
source("code_networks/scTenifoldTime.R")
dta <- readRDS("generate_networks/GSE147807_CMs/data/dta_raw.rds")
dta_sudotime <- readRDS("generate_networks/GSE147807_CMs/data/sudotime.rds")
dta <- scTenifoldNet::scQC(as.matrix(dta),  maxMTratio = 0.5, minPCT = 0.2)
dta <- PR_Normalization(as.matrix(dta))
n_cell <- ncol(dta)
n_gene <- nrow(dta)
gene_name <- rownames(dta)

#### Do scTenifoldTime
dta_list <- list()
n_group <- 15
time_vec <- rep(NA, n_group)
cell_diff <- round(n_cell / (n_group + 1))
if (cell_diff > 250) {
  for (i in seq_len(n_group - 1)) {
    dta_list[[i]] <- as.matrix(dta[, (cell_diff * (i - 1) + 1):(cell_diff * (i + 1))])
    time_vec[i] <- mean(dta_sudotime[colnames(dta_list[[i]]), 1])
  }
  dta_list[[n_group]] <- as.matrix(dta[,  (cell_diff * (n_group - 1) + 1):n_cell])
  time_vec[n_group] <- mean(dta_sudotime[colnames(dta_list[[n_group]]), 1])
  time_vec <- time_vec / max(time_vec)
} else{
  cell_diff <- round((n_cell - 500) / (n_group - 1))
  if (cell_diff < 5) stop("The number of cell is too small.")
  for (i in seq_len(n_group - 1)) {
    dta_list[[i]] <- as.matrix(dta[, (cell_diff * (i - 1) + 1):(cell_diff * (i - 1) + 500)])
    time_vec[i] <- mean(dta_sudotime[colnames(dta_list[[i]]), 1])
  }
  dta_list[[n_group]] <- as.matrix(dta[, (n_cell - 500 + 1):n_cell])
  time_vec[n_group] <- mean(dta_sudotime[colnames(dta_list[[n_group]]), 1])
  time_vec <- time_vec / max(time_vec)
}

#### save results
res <- scTenifoldTime(dta_list, method = "pcnet", time_vec, nComp = 5, q = 0,
                      K = 5, maxIter = 10000, maxError = 1e-5)
saveRDS(res, paste0("generate_networks/GSE147807_CMs/results/networks/res_pcnet_ngroup", n_group, ".rds"))