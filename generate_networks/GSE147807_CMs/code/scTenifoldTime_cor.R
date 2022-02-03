# #### load data and packages
# setwd("/data/victorxun/scTenifoldTime/inst/GSE147807_CMs/")
# library(scTenifoldTime)
#
# dta <- readRDS("data/dta_raw.rds")
# dta_sudotime <- readRDS("data/sudotime.rds")
# dta <- scTenifoldNet::scQC(as.matrix(dta),  maxMTratio = 0.5, minPCT = 0.2)
# dta <- PR_Normalization(as.matrix(dta))
# n_cell <- ncol(dta)
# n_gene <- nrow(dta)
# gene_name <- rownames(dta)
#
# #### Do scTenifoldTime
# dta_list <- list()
# time_vec <- rep(NA, 10)
# for (i in 1:9) {
#   dta_list[[i]] <- as.matrix(dta[, (100 * (i - 1) + 1):(100 * (i - 1) + 500)])
#   time_vec[i] <- mean(dta_sudotime[colnames(dta_list[[i]]), 1])
# }
# dta_list[[10]] <- as.matrix(dta[, (n_cell - 500 + 1):n_cell])
# time_vec[10] <- mean(dta_sudotime[colnames(dta_list[[10]]), 1])
# time_vec <- time_vec / max(time_vec)
#
# #### correlation method
# res <- scTenifoldTime(dta_list, method = "cor", time_vec, nComp = 5, q = 0,
#                       K = 5, maxIter = 10000, maxError = 1e-5)
#
# saveRDS(res, "results/networks/res_cor.rds")
