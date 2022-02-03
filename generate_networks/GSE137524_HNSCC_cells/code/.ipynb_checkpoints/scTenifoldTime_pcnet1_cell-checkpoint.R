## load data and packages
setwd("/data/victorxun/scTenifoldTime/inst/GSE137524_HNSCC_cells/")
library(scTenifoldTime)
library(Seurat)

dta <- readRDS("data/dta_rep1.rds")
sudotime <- readRDS("data/sudotime1.rds")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps
dta_seurat <- CreateSeuratObject(counts = dta)
dta_seurat <- NormalizeData(dta_seurat)
dta_seurat <- FindVariableFeatures(dta_seurat, selection.method = "vst")
dta_seurat <- ScaleData(dta_seurat, features = rownames(dta_seurat))

dta_seurat <- CellCycleScoring(dta_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
dta_seurat <- ScaleData(dta_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(dta_seurat))

# test 1
dta <- as.matrix(dta_seurat@assays$RNA@scale.data)

n_cell <- ncol(dta)
n_gene <- nrow(dta)
gene_name <- rownames(dta)

## Do scTenifoldTime
dta_list <- list()
time_vec <- rep(NA, 10)
cell_diff <- 467
for (i in seq_len(9)) {
  dta_list[[i]] <- as.matrix(dta[, (cell_diff * (i - 1) + 1):(cell_diff * (i + 1))])
  time_vec[i] <- mean(sudotime[colnames(dta_list[[i]]), 3])
}
dta_list[[10]] <- as.matrix(dta[, 4204:n_cell])
time_vec[10] <- mean(sudotime[colnames(dta_list[[10]]), 3])
time_vec <- time_vec / max(time_vec)

## return results
res <- scTenifoldTime(dta_list = dta_list, method = "pcnet", time_vec = time_vec,
                      nComp = 5, q = 0, K = 5)
saveRDS(res, "results/networks/res_pcnet1_v1.rds")

# test 2
dta_seurat$CC.Difference <- dta_seurat$S.Score - dta_seurat$G2M.Score
dta_seurat <- ScaleData(dta_seurat, vars.to.regress = "CC.Difference", features = rownames(dta_seurat))

dta <- as.matrix(dta_seurat@assays$RNA@scale.data)

n_cell <- ncol(dta)
n_gene <- nrow(dta)
gene_name <- rownames(dta)

## Do scTenifoldTime
dta_list <- list()
time_vec <- rep(NA, 10)
cell_diff <- 467
for (i in seq_len(9)) {
  dta_list[[i]] <- as.matrix(dta[, (cell_diff * (i - 1) + 1):(cell_diff * (i + 1))])
  time_vec[i] <- mean(sudotime[colnames(dta_list[[i]]), 3])
}
dta_list[[10]] <- as.matrix(dta[, 4204:n_cell])
time_vec[10] <- mean(sudotime[colnames(dta_list[[10]]), 3])
time_vec <- time_vec / max(time_vec)

## return results
res <- scTenifoldTime(dta_list = dta_list, method = "pcnet", time_vec = time_vec,
                      nComp = 5, q = 0, K = 5)
saveRDS(res, "results/networks/res_pcnet1_v2.rds")
