# ## load data and packages
# setwd("/data/victorxun/scTenifoldTime/inst/GSE137524_HNSCC_cells/")
# library(scTenifoldTime)
#
# dta <- readRDS("data/dta_rep1.rds")
# sudotime <- readRDS("data/sudotime1.rds")
# # dta <- scTenifoldNet::scQC(as.matrix(dta), minPCT = 0.25)
# dta <- PR_Normalization(as.matrix(dta))
# n_cell <- ncol(dta)
# n_gene <- nrow(dta)
# gene_name <- rownames(dta)
#
# ## Do scTenifoldTime
# dta_list <- list()
# time_vec <- rep(NA, 10)
# cell_diff <- 488
# for (i in seq_len(9)) {
#   dta_list[[i]] <- as.matrix(dta[, (cell_diff * (i - 1) + 1):(cell_diff * (i + 1))])
#   time_vec[i] <- mean(sudotime[colnames(dta_list[[i]]), 9])
# }
# dta_list[[10]] <- as.matrix(dta[, 4393:n_cell])
# time_vec[10] <- mean(sudotime[colnames(dta_list[[10]]), 9])
# time_vec <- time_vec / max(time_vec)
#
#
# ## return results
# res <- scTenifoldTime(dta_list = dta_list, method = "cor", time_vec = time_vec,
#                       nComp = 5, q = 0, K = 5)
# saveRDS(res, "results/networks/res_cor1.rds")
