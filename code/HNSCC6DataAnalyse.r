library(Seurat)
library(monocle3)
library(scTenifoldNet)
library(ggpubr)
library(statsExpressions)
library(tidyverse)
library(ggplot2)
library(viridis)
library(GSVA)
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)

setwd("Z:/Cailab/scTenifoldTime/")
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/filterGenes.R')
SCC6<- read.csv('Data/HNSCC/GSE137524_exprsSCC6Matrix.csv')
rownames(SCC6) <- make.unique(SCC6$X)
SCC6 <- SCC6[,-c(1)]
SCC6 <- scQC(as.matrix(SCC6))
SCC6 <- filterGenes(SCC6)
cellnumber <- colnames(SCC6)
SCC6_CTX_R1ID <- which(grepl('SCC6_CTX_R1',cellnumber))
SCC6_PBS_R1ID <- which(grepl('SCC6_PBS_R1',cellnumber))

CTX <- SCC6[,c(SCC6_CTX_R1ID)]
PBS <- SCC6[,c(SCC6_PBS_R1ID)]
rm(SCC6_CTX_R1ID,SCC6_PBS_R1ID)
CTX <- CreateSeuratObject(count=CTX,project = "CTX")
PBS <- CreateSeuratObject(count=PBS,project = "PBS")
A <- merge(CTX,PBS)
A <- NormalizeData(A)
A <- FindVariableFeatures(A)
A <- ScaleData(A)
A <- RunPCA(A)
DimPlot(A, reduction = "pca")
ElbowPlot(A)
A <- RunUMAP(A, dims = 1:10)
Pumap <- UMAPPlot(A) +
  xlab("UMAP 1") + ylab("UMAP 2")
Pumap
A <- RunTSNE(A, dims = 1:10)
ptsne <- DimPlot(A,reduction = 'tsne') +
  xlab("tSNE 1") + ylab("tSNE 2")
ptsne

FeaturePlot(A,features = 'ERBB2',pt.size = 1,reduction = 'tsne')+
  xlab("tSNE 1") + ylab("tSNE 2")+ scale_colour_viridis()
FeaturePlot(A,features = 'AREG',pt.size = 1,reduction = 'tsne')+
  xlab("tSNE 1") + ylab("tSNE 2")+ scale_colour_viridis()

FeaturePlot(A,features = 'CXCL8',pt.size = 1,reduction = 'tsne')+
  xlab("tSNE 1") + ylab("tSNE 2")+ scale_colour_viridis()
FeaturePlot(A,features = 'VIM',pt.size = 1,reduction = 'tsne')+
  xlab("tSNE 1") + ylab("tSNE 2")+ scale_colour_viridis()

FeaturePlot(A,features = 'MMP2',pt.size = 1,reduction = 'tsne')+
  xlab("UMAP 1") + ylab("UMAP 2")+ scale_colour_viridis()

FeaturePlot(A,features = 'KI67',pt.size = 1,reduction = 'tsne')+
  xlab("UMAP 1") + ylab("UMAP 2")+ scale_colour_viridis()
FeaturePlot(A,features = c('PIK3CA','KRAS','HRAS','PTEN'),pt.size = 1,reduction = 'tsne')
FeaturePlot(A,features = c('RANKL'),pt.size = 1,reduction = 'tsne')
HVG <- VariableFeatures(A)

###########  Trajectory analysis  ##################
rm(CTX,PBS)
All <- A@assays$RNA@counts
cell_metadata <- as.data.frame(cbind(colnames(All)))
cell_metadata$Treatment <- "NA"
cell_metadata$Treatment[which(grepl(".*SCC6_CTX_R1",colnames(All)))] <- 'CTX'
cell_metadata$Treatment[which(grepl(".*SCC6_PBS_R1",colnames(All)))] <- 'PBS'
colnames(cell_metadata) <- c('barcode','Treatment')
rownames(cell_metadata) <- cell_metadata$barcode
gene_metadata <- as.data.frame(cbind(rownames(All)))
rownames(gene_metadata) <- gene_metadata$V1
colnames(gene_metadata) <- 'gene_short_name'
cds <- new_cell_data_set(expression_data = All,gene_metadata = gene_metadata,cell_metadata=cell_metadata)
cds <- preprocess_cds(cds, num_dim = 50, method ="PCA")
cds <- reduce_dimension(cds,reduction_method='UMAP',umap.min_dist = 0.2)
cds <- cluster_cells(cds, k=30)
cds <- learn_graph(cds, verbose = FALSE,use_partition = F)
P6 <- plot_cells(cds, label_cell_groups = F,label_branch_points = F,label_leaves = F,color_cells_by = 'Treatment')
P6

VIMID <- which(rownames(assays(cds)[[1]])=='VIM')
IVLID <- which(rownames(assays(cds)[[1]])=='IVL')
CDH1ID <- which(rownames(assays(cds)[[1]])=='CDH1')
CDH2ID <- which(rownames(assays(cds)[[1]])=='CDH2')
colData(cds)$root.cells <- grepl(".*SCC6_PBS_R1",colnames(assays(cds)[[1]])) & assays(cds)[[1]][VIMID,] ==0 & assays(cds)[[1]][IVLID,]==0 & assays(cds)[[1]][CDH1ID,]>3 & assays(cds)[[1]][CDH2ID,]==0 & cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_dist"]][1,] < -5 & cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_dist"]][2,] < 0
length(which(colData(cds)$root.cells))
which(colData(cds)$root.cells)
get_earliest_principal_node <- function(cds) {
  cell_ids <- which(cds$root.cells)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, 
  ]))))]
  root_pr_nodes
}
root.nodes <- get_earliest_principal_node(cds)
cds <- order_cells(cds, root_pr_nodes = root.nodes)
P8 <- plot_cells(cds,reduction_method = 'UMAP', color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 0.2,label_roots =F,cell_size = 0.7)
P8

A$pseudotime <- pseudotime(cds)[match(colnames(A), colnames(cds))]
umap.pt <- FeaturePlot(A, feature = "pseudotime")
umap.pt <- umap.pt + scale_colour_viridis(option = "plasma")+labs(title = 'Pseudotime')
umap.pt

tsne.pt <- FeaturePlot(A, feature = "pseudotime",reduction = 'tsne')
tsne.pt <- tsne.pt + scale_colour_viridis(option = "plasma")+labs(title = 'Pseudotime')
tsne.pt

SCC6_repeat1 <- as.data.frame(pData(cds))
SCC6_repeat1$Pseudotime <- pseudotime(cds)[match(SCC6_repeat1$barcode,colnames(cds))]
#write.csv(SCC6_repeat1,file = 'Data/HNSCC/Ordered_SCC6_repeat1_r4.csv')

violin <- ggplot(SCC6_repeat1, aes(x=Treatment, y=Pseudotime, fill = Treatment)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.03, fill="white")+
  labs(y='Pseudotime') + theme_minimal() + stat_compare_means(method = "t.test",label.x = 1.2, label.y = 40)+scale_fill_manual(values=c('#F8766D','#00BFC4'))
violin


genes<- c('VIM','MT2A','PFN1','CAV1','FN1')
cds_subset <- cds[rownames(cds) %in% genes,]
png('figures/HNSCC/plot_sudotime_rep1.png',res = 1000,height = 5000,width = 4000)
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="pseudotime",
                         min_expr=0)
dev.off()


#####################   find genes changing across time using monocle ##############
gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime",expression_family="negbinomial")
fit_coefs <- coefficient_table(gene_fits)
pseudotime_terms <- fit_coefs %>% filter(term %in% c("pseudotime"))
pseudotime_terms <- pseudotime_terms %>%
  select(gene_short_name, term, q_value, estimate)
#write.csv(pseudotime_terms,file = 'figures/HNSCC/DE_by_monocle3.csv',  quote = F,row.names = F,)

genes<- c('ERBB2' )
genes<- c('VIM','MT2A','CLTB','CAV1' )
cds_subset <- cds[rownames(cds) %in% genes]
png('figures/HNSCC/plot_sudotime_rep2.png',res = 800,height = 5000,width = 4000)
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="pseudotime",
                         min_expr=0)
dev.off()

#Construct networks and regression to get the beta matrix
#SCC1
library(scTenifoldTime)
dta <- readRDS("data/dta_rep1.rds")
sudotime <- readRDS("data/sudotime1.rds")
# dta <- scTenifoldNet::scQC(as.matrix(dta), minPCT = 0.25)
dta <- PR_Normalization(as.matrix(dta))
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
# saveRDS(res, "results/networks/res_pcnet1.rds")


################################### Downstream analysis starting from beta-matrix ###################################
setwd("Z:/Cailab/scTenifoldTime/")
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
res_pcnet1 <- readRDS("Z:/Cailab/scTenifoldTime/network/HNSCC6/newest/res_pcnet1.rds")
beta <- res_pcnet1$beta_mat
set.seed(123)
beta_mat_svd <- RSpectra::svds(beta, k = 50)
beta_mat_pcscore <-  t(t(beta_mat_svd$u) * beta_mat_svd$d)
TSNE <- Rtsne::Rtsne(beta_mat_pcscore, pca = FALSE, check_duplicates = TRUE)$Y
gene_list <- mahalanobis(TSNE, center = colMeans(TSNE), cov = cov(TSNE))
names(gene_list) <- toupper(rownames(beta))
gene_list <- sort(gene_list, decreasing = TRUE)
U <-gene_list
U <- U[!grepl('RPL|RPS|RP[[:digit:]]+|MT-', names(U))]
eD <- fgseaMultilevel(BIOP, U)
eD$leadingEdge <- unlist(lapply(eD$leadingEdge, function(X){paste0(X,collapse = ';')}))
eD <- eD[eD$padj < 0.05,,drop=FALSE]
#write.csv(eD, file = paste0('figures/HNSCC/SCC6_rep1_gsea_BIOP.csv')) 

gSet <- 'Beta-3 integrin cell surface interactions'
pTitle <- 'Beta-3 integrin cell\nsurface interactions'
plotEnrichment(BIOP[[gSet]], U,ticksSize = 0.001) + 
  labs(title = pTitle) +
  theme_bw() + 
  ylab('Enrichment score') +
  xlab('Gene rank') +
  theme(plot.title = element_text(face = 2, size = 20))

# plot the example networks that change across pseudotime
df <- eD[eD$pathway=='Beta-3 integrin cell surface interactions',]
genes <- unlist(strsplit(df$leadingEdge,split = ';'))
Networks <- list()
for (i in 1:10){
  Net <- res_pcnet1$network_list[[i]]
  Net <- Net[colnames(Net) %in% genes,rownames(Net) %in% genes]
  Networks[[i]] <- Net
}
Network <- as.matrix(cbind(Networks[[1]],Networks[[3]],Networks[[5]],
                           Networks[[7]],Networks[[9]]))

col_fun = colorRamp2(c(-0.2,0.2), c('green',  "#e34a33"))
split = rep(c('Time 1','Time 3',"Time 5", "Time 7",'Time 9'), each=length(genes))
Heatmap(Network, name = "Beta", col = col_fun,cluster_rows = F,cluster_columns = FALSE,
        column_split = split,column_title_gp = gpar(fill = c("#fff7ec", "#fee8c8", "#fdd49e","#fdbb84","#fc8d59")),
        column_gap  = unit(3, "mm") )

# plot the gene set in the tsne plot
TSNE.orig <- TSNE
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!TSNE$X %in% genes] <- NA
TSNE$col <- ifelse(is.na(TSNE$X), '#2ca25f', 'red')
TSNE$size <- ifelse(is.na(TSNE$X),1,2.5)
ggplot(TSNE, aes(`TSNE  1`, `TSNE 2`, label = X)) + geom_point(color = TSNE$col,alpha=0.5,size=TSNE$size) + 
  geom_label_repel(aes(label = X),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',max.overlaps = 20) +
  theme_classic()

TSNE <- as.data.frame(TSNE.orig)
rownames(TSNE) <- rownames(beta)
colnames(TSNE) <- c('TSNE1','TSNE2')
set.seed(123)
fit_cluster_kmeans=kmeans(TSNE, 50,iter.max = 1000) 
TSNE$Group = factor(fit_cluster_kmeans$cluster)
ggplot() +
  geom_point(data = TSNE, 
             mapping = aes(x = TSNE$TSNE1, 
                           y = TSNE$TSNE2, 
                           colour = TSNE$Group)) +
  geom_point(mapping = aes_string(x = fit_cluster_kmeans$centers[, "TSNE1"], 
                                  y = fit_cluster_kmeans$centers[, "TSNE2"]),
             color = "black", size = 0.02)+
  geom_text(mapping = aes_string(x =fit_cluster_kmeans$centers[, "TSNE1"], 
                                 y =fit_cluster_kmeans$centers[, "TSNE2"],
                                 label = 1:50),
            color = "black", size = 4) +
  theme_light()+
  xlab("TSNE 1") + ylab("TSNE 2")+ labs("")+ NoLegend()

TSNE$gene <- rownames(beta)
TSNE2 <- TSNE

############## MYC ##############
group = '37'
geneid <- which(TSNE2$Group==group,) 
genes <- rownames(TSNE2)[geneid]
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!TSNE$X %in% genes] <- NA
TSNE$col <- ifelse(is.na(TSNE$X), '#2ca25f', 'red')
TSNE$size <- ifelse(is.na(TSNE$X),0.5,2)
ggplot(TSNE, aes(`TSNE  1`, `TSNE 2`, label = X)) + geom_point(color = TSNE$col,alpha=0.5,size=TSNE$size) +
  theme_classic()+ggtitle(paste0('Cluster ',group))

writeLines(genes)

############## VIM ##############
gene <- 'VIM'
vim <- beta[rownames(beta)==gene,]
vim <- vim[!grepl(pattern = 'RPS|RPL|RP[[:digit:]]+',names(vim))]
vim <- vim[order(abs(vim),decreasing = T)]
genes <- names(vim)[0:30]
Networks <- list()
for (i in 1:10){
  Net <- res_pcnet1$network_list[[i]]
  Net <- Net[rownames(Net) %in% genes,colnames(Net) == gene]
  Networks[[i]] <- Net
}
Network <- as.matrix(rbind(Networks[[1]],Networks[[2]],Networks[[3]],Networks[[4]],Networks[[5]],
                           Networks[[6]],Networks[[7]],Networks[[8]],Networks[[9]],Networks[[10]]))

col_fun = colorRamp2(range(Network), c('green',  "red"))
Heatmap(Network, name = "Beta", col = col_fun,cluster_rows = F,cluster_columns = T)


genes <- c(genes,'VIM')
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!TSNE$X %in% genes] <- NA
TSNE$col <- ifelse(is.na(TSNE$X), '#2ca25f', 'red')
TSNE$size <- ifelse(is.na(TSNE$X),1,2)
ggplot(TSNE, aes(`TSNE  1`, `TSNE 2`, label = X)) + geom_point(color = TSNE$col,alpha=0.5,size=TSNE$size)+ 
  geom_label_repel(aes(label = X),
                   box.padding   = 0.05, 
                   point.padding = 0.5,
                   segment.color = 'grey50',max.iter = 5000,max.overlaps = 300)+theme_classic()
