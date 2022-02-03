## load data
dta <- read.csv("data/data_original/Cardiomyocytes.csv")
sudotime <- read.csv("data/data_original/Barcode_pseudotime_r4.csv")

## remove replicated genes
gene_name <- dta$X
dta <- dta[, -1]
rownames(dta) <- gene_name

## save data set
rownames(sudotime) <- sudotime$V1
sudotime <- sudotime[, -1, drop = FALSE]
sudotime <- sudotime[order(sudotime$PT, decreasing = FALSE), , drop = FALSE]
dta <- dta[, rownames(sudotime)]

saveRDS(dta, "data/dta_raw.rds")
saveRDS(sudotime, "data/sudotime.rds")
