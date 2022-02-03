library(Seurat)
library(Matrix)

dta <- Read10X("data/data_original/")
sudotime <- read.csv("data/data_original/Zebrafish_Pseudotime_r4.0.csv")
rownames(sudotime) <- sudotime$barcode
sudotime <- sudotime[, -1]
sudotime <- sudotime[order(sudotime$Pseudotime, decreasing = FALSE), ]
dta <- dta[, rownames(sudotime)]

saveRDS(dta, "data/dta.rds")
saveRDS(sudotime, "data/sudotime.rds")
