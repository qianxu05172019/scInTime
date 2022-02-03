# ## load data and packages
# setwd("/data/victorxun/scTenifoldTime/inst/GSE137524_HNSCC_cells/")
# library(scTenifoldTime)
#
# dta <- readRDS("data/dta_rep2.rds")
# sudotime <- readRDS("data/sudotime2.rds")
# # dta <- scTenifoldNet::scQC(as.matrix(dta), minPCT = 0.25)
# dta <- PR_Normalization(as.matrix(dta))
# n_cell <- ncol(dta)
# n_gene <- nrow(dta)
# gene_name <- rownames(dta)
#
# ## Do scTenifoldTime
# dta_list <- list()
# time_vec <- rep(NA, 10)
# cell_diff <- 346
# for (i in 1:9) {
#   dta_list[[i]] <- as.matrix(dta[, (346 * (i - 1) + 1):(346 * (i + 1))])
#   time_vec[i] <- mean(sudotime[colnames(dta_list[[i]]), 3])
# }
# dta_list[[10]] <- as.matrix(dta[, 3115:n_cell])
# time_vec[10] <- mean(sudotime[colnames(dta_list[[10]]), 3])
# time_vec <- time_vec / max(time_vec)
#
# ## return results
# res <- scTenifoldTime(dta_list = dta_list, method = "cor", time_vec = time_vec,
#                       nComp = 5, q = 0, K = 5)
# saveRDS(res, "results/networks/res_cor2.rds")
