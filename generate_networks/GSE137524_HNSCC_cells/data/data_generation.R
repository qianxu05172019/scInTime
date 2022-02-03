## load data
dta1 <- read.csv("data/data_original/SCC6_repeat1.csv")
dta2 <- read.csv("data/data_original/SCC6_repeat2.csv")
sudotime1 <- read.csv("data/data_original/Ordered_SCC6_repeat1_r4.csv")
sudotime2 <- read.csv("data/data_original/Ordered_SCC6_repeat2_r4.csv")

## remove replicated genes
gene_name <- dta1$X
dta1 <- dta1[, -1]
index_dup <- which(duplicated(gene_name)==TRUE)
rownames(dta1) <- gene_name

gene_name <- dta2$X
dta2 <- dta2[, -1]
index_dup <- which(duplicated(gene_name)==TRUE)
rownames(dta2) <- gene_name

## ordered the cells
colnames(dta1) <- gsub(pattern = "[.]", replacement = "_", x = colnames(dta1))
rownames(sudotime1) <- gsub(pattern = "[.]", replacement = "_", x = sudotime1$X)
sudotime1 <- sudotime1[, -c(1,2)]
sudotime1 <- sudotime1[order(sudotime1$Pseudotime, decreasing = FALSE), ]
dta1 <- dta1[, rownames(sudotime1)]


colnames(dta2) <- gsub(pattern = "[.]", replacement = "_", x = colnames(dta2))
rownames(sudotime2) <- gsub(pattern = "[.]", replacement = "_", x = sudotime2$X)
sudotime2 <- sudotime2[, -c(1,2)]
sudotime2 <- sudotime2[order(sudotime2$Pseudotime, decreasing = FALSE), ]
dta2 <- dta2[, rownames(sudotime2)]

saveRDS(dta1, file = paste0("data/dta_rep1.rds"))
saveRDS(dta2, file = paste0("data/dta_rep2.rds"))
saveRDS(sudotime1, file = paste0("data/sudotime1.rds"))
saveRDS(sudotime2, file = paste0("data/sudotime2.rds"))
