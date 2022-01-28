library(Seurat)
library(ggplot2)
library(viridis)
library(monocle3)
library(statsExpressions)
library(tidyverse)
library(GSVA)
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(Rtsne)
library(phateR)

set.seed(1)
setwd("Z:/Cailab/scTenifoldTime/")
hbp <- Read10X(data.dir = "Data/Zebrafish/", gene.column=1)
annotation <- read.csv(file = 'Data/Zebrafish/Zebrafish_Pseudotime.csv')
hbp <- hbp[,colnames(hbp) %in% annotation$barcode]
hbp <- CreateSeuratObject(counts = hbp)
hbp$stage <- annotation$stage[match(colnames(hbp),annotation$barcode)]
hbp$cell <- annotation$celltype[match(colnames(hbp),annotation$barcode)]
hbp <- NormalizeData(hbp, normalization.method = "LogNormalize", scale.factor = 10000)
hbp <- FindVariableFeatures(hbp, selection.method = "vst", nfeatures = 2000)
hbp <- ScaleData(hbp, features = rownames(hbp))
hbp <- RunPCA(hbp, features = VariableFeatures(object = hbp), npcs = 50, verbose = FALSE, 
              ndims.print = NULL, nfeatures.print = NULL, reduction.key = "PC_", seed.use = 1, 
              approx = FALSE, reduction.name = "pca")
hbp <- FindNeighbors(hbp, dims = 1:50, verbose = FALSE, k.param = 20, compute.SNN = TRUE, 
                     prune.SNN = 1/15, nn.eps = 0)
hbp <- FindClusters(hbp, verbose = FALSE, resolution = 0.8)
hbp <- RunUMAP(hbp, dims = 1:50, verbose = FALSE,seed.use = 1)
hbp <- RunTSNE(hbp,dims = 1:50,seed.use = 1)

umap_stage <- DimPlot(hbp, reduction = "umap", group.by = "stage", label = TRUE, 
              label.size = 4) + NoLegend()+xlab("UMAP 1") + ylab("UMAP 2")
umap_stage
tsne_stage <- DimPlot(hbp,reduction = 'tsne',group.by = 'stage',label = TRUE, 
                label.size = 4)+xlab("tSNE 1") + ylab("tSNE 2")+ NoLegend()
tsne_stage
umap_celltype <- DimPlot(hbp, reduction = "umap", group.by = "cell", label = TRUE, 
              label.size = 4) + NoLegend()+ xlab("UMAP 1") + ylab("UMAP 2")
umap_celltype
tsne_celltype <- DimPlot(hbp, reduction = "tsne", group.by = "cell", label = TRUE, 
              label.size = 4) + NoLegend()+ xlab("tSNE 1") + ylab("tSNE 2")
tsne_celltype
FeaturePlot(hbp, features = c("mki67", 'sox3','elavl3'), 
                   cols = viridis(10, direction = 1), ncol = 3,reduction = 'tsne')

HVG <- VariableFeatures(hbp)

##Trajectory analysis
A <- hbp@assays$RNA@counts
cell_metadata <- as.data.frame(cbind(colnames(A)))
cell_metadata$stage <- annotation$stage[match(cell_metadata$V1,annotation$barcode)]
cell_metadata$cellType <- annotation$celltype[match(cell_metadata$V1,annotation$barcode)]
colnames(cell_metadata) <- c('barcode','stage','cellType')
rownames(cell_metadata) <- cell_metadata$barcode
rownames(annotation) <- annotation$barcode
gene_metadata <- as.data.frame(cbind(rownames(A)))
rownames(gene_metadata) <- gene_metadata$V1
colnames(gene_metadata) <- 'gene_short_name'
cds <- new_cell_data_set(expression_data = A,gene_metadata = gene_metadata,cell_metadata=cell_metadata)
cds <- preprocess_cds(cds, num_dim = 50, method ="PCA")
cds <- align_cds(cds, num_dim = 50, alignment_group = "stage")
cds <- reduce_dimension(cds,reduction_method='UMAP',umap.min_dist = 0.2)
cds <- cluster_cells(cds, k=30)
cds <- learn_graph(cds,  verbose = FALSE)
monocle_stage <- plot_cells(cds, color_cells_by="stage", label_cell_groups = F,label_branch_points = F,label_leaves = F)
monocle_stage
monocle_celltype <- plot_cells(cds,
                  color_cells_by = "cellType",
                  label_groups_by_cluster=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE)
monocle_celltype

sox2ID <- which(rownames(assays(cds)[[1]])=='sox2')
egr2bID <- which(rownames(assays(cds)[[1]])=='egr2b')
mafbaID <- which(rownames(assays(cds)[[1]])=='mafba')
mki67ID <- which(rownames(assays(cds)[[1]])=='mki67')
colData(cds)$root.cells <- grepl(".*-1",colnames(assays(cds)[[1]])) & assays(cds)[[1]][mki67ID,] >0 & assays(cds)[[1]][sox2ID,]>0 & assays(cds)[[1]][egr2bID,]>0 & assays(cds)[[1]][mafbaID,]>0 & cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_dist"]][1,] < -5 & cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_dist"]][2,] < -2
length(which(colData(cds)$root.cells))
which(colData(cds)$root.cells)

monocle_rootcell <- plot_cells(cds, color_cells_by = "root.cells",cell_size = 0.7,label_cell_groups = F,label_branch_points = F,label_leaves = F,label_roots = F)
monocle_rootcell
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

monocle_PT <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2,label_roots =F,cell_size = 0.7)
monocle_PT
genes <- c("dla",'dlb' ,'her4.4' )
cds_subset <- cds[rowData(cds)$gene_short_name %in% genes,]
#png('figures/Zebrafish/plot_sudotime.png',res = 800,height = 2500,width = 4000)
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="pseudotime",
                         min_expr=0)
#dev.off() 
#gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime+stage",expression_family="negbinomial")
gene_fits <- fit_models(cds_subset, model_formula_str = "~pseudotime",expression_family="negbinomial")
fit_coefs <- coefficient_table(gene_fits)
pseudotime_terms <- fit_coefs %>% filter(term %in% c("pseudotime",'stagehpf24','stagehpf44'))
pseudotime_terms <- pseudotime_terms %>%
  select(gene_short_name, term, q_value, estimate)

hbp$pseudotime <- pseudotime(cds)[match(colnames(hbp), colnames(cds))]
umap.pt <- FeaturePlot(hbp, feature = "pseudotime")
umap.pt <- umap.pt + scale_colour_viridis(option = "plasma")+labs(title = 'Pseudotime')
umap.pt

tsne.pt <- FeaturePlot(hbp, feature = "pseudotime",reduction = 'tsne')
tsne.pt <- tsne.pt + scale_colour_viridis(option = "plasma")+labs(title = 'Pseudotime')+xlab("tSNE 1") + ylab("tSNE 2")
tsne.pt

#CCAT
Ref <- read.csv("Data/Zebrafish/Human_entrz_esembl.Ref.csv")
Ref <- Ref[Ref$entrezIDs!='nooutput',]
Aname <- rownames(hbp)
Aname <- toupper(Aname)
AID <- c()
for (p in 1 : length(Aname)){
  GENENAME <- Ref$entrezIDs[which(Ref$SYMBOL==Aname[p]) ]
  if (length(GENENAME)!=0){
    AID[p] <-GENENAME
  }else{
    AID[p] <- 'X'             
  }
}
hbp@assays$RNA@counts@Dimnames[[1]] <- AID
Amtrix <- hbp@assays$RNA@counts

load("Data/Zebrafish/net17Jan16.rda")
ppiA.m <- net17Jan16.m
exp.m <- log2(Amtrix+1)
exp.m <- as.matrix(exp.m)
commonEID.v <- intersect(rownames(ppiA.m),rownames(exp.m))
k.v <- rowSums(ppiA.m[match(commonEID.v, rownames(ppiA.m)),])
ccat.v <- as.vector(cor(exp.m[match(commonEID.v,rownames(exp.m)),],k.v))

cell <- colnames(exp.m)
CCAT.celltype <-as.data.frame(cbind(cell,ccat.v))
CCAT.celltype$celltype <- annotation$cell[match(CCAT.celltype$cell, annotation$barcode)]
hbp$CCAT <-rank( as.numeric(CCAT.celltype$ccat.v[match(colnames(hbp), CCAT.celltype$cell)]))
CCAT_tsne <- FeaturePlot(hbp, feature = "CCAT",reduction = 'tsne')+ scale_colour_viridis(option = "plasma")
CCAT_tsne

CCAT.celltype$rank <- rank(CCAT.celltype$ccat.v)
celltype <- as.character(CCAT.celltype$celltype)
CCAT.celltype$celltype <- with(CCAT.celltype,reorder(celltype,rank,median))
ordercell <- as.data.frame(pseudotime(cds))
ordercell$barcode <- rownames(ordercell)
ordercell <- cbind(ordercell,colData(cds)$stage)

ordercell  <- merge(ordercell,CCAT.celltype,by.x='barcode',by.y='cell')
colnames(ordercell) <- c("barcode" ,"Pseudotime","stage","ccat.v","celltype","CCAT_rank")
ordercell <- ordercell[!ordercell$Pseudotime=='Inf',]
levels(CCAT.celltype$celltype)

violin_PT <- ggplot(ordercell, aes(y=celltype, x=rank(Pseudotime), fill=celltype)) + 
  geom_violin()+ theme_bw()+ 
  theme(legend.position = 'none')+
  geom_jitter(position = position_jitter(width = 0.5, height = 0.3),size=0.7,alpha=0.20)+
  scale_y_discrete(limits = c('Progenitors (16hpf)', "DP (24hpf)","VMP (24hpf)","Differentiating Progenitors (44hpf)",
                              "Immature Neurons", "Hindbrain Neurons" ))+ NoLegend()+xlab("")+ylab("")+theme(axis.text=element_text(size=12))
violin_PT

violin_CCAT <- ggplot(ordercell, aes(y=celltype, x=CCAT_rank,fill=celltype)) + 
  geom_violin()+ theme_bw()+ 
  theme(legend.position = 'none')+
  geom_jitter(position = position_jitter(width = 0.5, height = 0.3),size=0.7,alpha=0.20)+
  scale_y_discrete(limits = c('Progenitors (16hpf)', "DP (24hpf)","VMP (24hpf)","Differentiating Progenitors (44hpf)",
                              "Immature Neurons", "Hindbrain Neurons" ))+ NoLegend()+xlab("")+ylab("")+theme(axis.text=element_text(size=12)) + scale_color_brewer(palette = "Dark2")

violin_CCAT

####################### Downstream analysis starting from beta-matrix #######################
setwd("Z:/Cailab/scTenifoldTime/")
STRINGDB <- read.table('figures/Zebrafish/enrichment.PMID.tsv',sep = '\t',header = T)
gSets <- strsplit(STRINGDB$matching.proteins.in.your.input..labels., ',')
names(gSets) <- STRINGDB$term.description
DE <- read.csv('figures/Zebrafish/DE_by_monocle3.csv')
DE <- DE$x
#res_pcnet <- readRDS("Z:/Cailab/scTenifoldTime/network/Zebrafish/res_pcnet.rds")
# r4
res_pcnet <- readRDS("Z:/Cailab/scTenifoldTime/network/Zebrafish/r4/results/networks/res_pcnet.rds")
beta <- res_pcnet$beta_mat
set.seed(123)
beta_mat_svd <- RSpectra::svds(beta, k = 50)
beta_mat_pcscore <-  t(t(beta_mat_svd$u) * beta_mat_svd$d)
TSNE <- Rtsne::Rtsne(beta_mat_pcscore, pca = FALSE, check_duplicates = TRUE)$Y
gene_list <- mahalanobis(TSNE, center = colMeans(TSNE), cov = cov(TSNE))
names(gene_list) <- toupper(rownames(beta))
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!grepl('RPL|RPS|RP[[:digit:]]+|MT-', names(gene_list))]
U <- gene_list
U <- U[!grepl('RPL|RPS|RP[[:digit:]]+|MT-', names(U))]
names(U) <- tolower(names(U))
eD <- fgseaMultilevel(gSets, U)
eD$leadingEdge <- unlist(lapply(eD$leadingEdge, function(X){paste0(X,collapse = ';')}))
eD <- eD[eD$padj < 0.05,,drop=FALSE]
gSet <- '(2010) Zebrafish Pou5f1-dependent transcriptional networks in temporal control of early development.'
pTitle <- 'Pou5f1-dependent\ntranscriptional network'
plotEnrichment(gSets[[gSet]], U,ticksSize = 0.01) + 
  labs(title = pTitle) +
  theme_bw() + 
  ylab('Enrichment score') +
  xlab('Gene rank') +
  theme(plot.title = element_text(face = 2, size = 20))

# plot the example networks that change across pseudotime
df <- eD[eD$pathway==gSet,]
genes <- unlist(strsplit(df$leadingEdge,split = ';'))
Networks <- list()
for (i in 1:10){
  Net <- res_pcnet$network_list[[i]]
  Net <- Net[colnames(Net) %in% genes,rownames(Net) %in% genes]
  Networks[[i]] <- Net
}
Network <- as.matrix(cbind(Networks[[1]],Networks[[3]],Networks[[5]],
                           Networks[[7]],Networks[[9]]))
range(Network)
col_fun = colorRamp2(c(-0.2,0.2), c('green',  "red"))
split = rep(c('Time 1','Time 3',"Time 5", "Time 7",'Time 9'), each=length(genes))
Heatmap(Network, name = "Beta", col = col_fun,cluster_rows = T,cluster_columns = FALSE,
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
                   segment.color = 'grey50',segment.size = 0.7) +
  theme_classic() + theme(
    plot.title = element_text(color="black", size=16)
  )

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
             color = "black", size = 0.2)+
  geom_text(mapping = aes_string(x =fit_cluster_kmeans$centers[, "TSNE1"], 
                                 y =fit_cluster_kmeans$centers[, "TSNE2"],
                                 label = 1:50),
            color = "black", size = 4) +
  theme_light()+
  xlab("TSNE 1") + ylab("TSNE 2")+ labs("")+ NoLegend()

TSNE$gene <- rownames(beta)
TSNE2 <- TSNE

############## MYC ##############
group = '4'
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
genes <- gsub(pattern = '\\(1 of many)',replacement = '',x = genes)
genes <- gsub(pattern = '\\:',replacement = '.',x = genes)
genes <- gsub(pattern = '\\ ',replacement = '',x = genes)
writeLines(genes)

############## gene regulators ##############
gene <- 'dlb'
vim <- beta[rownames(beta)==gene,]
vim <- vim[!grepl(pattern = 'rps|rpl|rp[[:digit:]]+',names(vim))]
vim <- vim[order(abs(vim),decreasing = T)]
genes <- names(vim)[0:30]
Networks <- list()
for (i in 1:10){
  Net <- res_pcnet$network_list[[i]]
  Net <- Net[rownames(Net) %in% genes,colnames(Net) == gene]
  Networks[[i]] <- Net
}
Network <- as.matrix(rbind(Networks[[1]],Networks[[2]],Networks[[3]],Networks[[4]],Networks[[5]],
                           Networks[[6]],Networks[[7]],Networks[[8]],Networks[[9]],Networks[[10]]))

col_fun = colorRamp2(range(Network), c('green',  "red"))
Heatmap(Network, name = "Beta", col = col_fun,cluster_rows = F,cluster_columns = T)

genes <- c(genes,gene)
writeLines(genes)
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!TSNE$X %in% genes] <- NA
TSNE$col <- ifelse(is.na(TSNE$X), '#2ca25f', 'red')
TSNE$size <- ifelse(is.na(TSNE$X),1,2)
ggplot(TSNE, aes(`TSNE  1`, `TSNE 2`, label = X)) + geom_point(color = TSNE$col,alpha=0.5,size=TSNE$size)+ 
  geom_label_repel(aes(label = X),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',max.iter = 50000)+theme_classic()

genes <- c("dla",'dlb' ,'her4.4' )
cds_subset <- cds[rowData(cds)$gene_short_name %in% genes,]
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="pseudotime",
                         min_expr=0)

######################## cluster 31 #############################
geneid <- which(TSNE2$Group==31,) 
genes <- rownames(TSNE2)[geneid]
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!TSNE$X %in% genes] <- NA
TSNE$col <- ifelse(is.na(TSNE$X), '#2ca25f', 'red')
TSNE$size <- ifelse(is.na(TSNE$X),1,2)
ggplot(TSNE, aes(`TSNE  1`, `TSNE 2`, label = X)) + geom_point(color = TSNE$col,alpha=0.5,size=TSNE$size)+ 
  geom_label_repel(aes(label = X),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',max.iter = 7000,segment.alpha = 0.5,label.size = 0.5)+theme_classic()

df <- eD[eD$pathway=='(2008) Genomewide expression analysis in zebrafish mind bomb alleles with pancreas defects of different severity identifies putative Notch responsive genes.',]
genes <- unlist(strsplit(df$leadingEdge,split = ';'))
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!TSNE$X %in% genes] <- NA
TSNE$col <- ifelse(is.na(TSNE$X), '#2ca25f', 'red')
TSNE$size <- ifelse(is.na(TSNE$X),1,2)
ggplot(TSNE, aes(`TSNE  1`, `TSNE 2`, label = X)) + geom_point(color = TSNE$col,alpha=0.5,size=TSNE$size)+ 
  geom_label_repel(aes(label = X),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',max.iter = 7000,segment.alpha = 0.5,label.size = 0.6)+theme_classic()+ggtitle('Genomewide expression analysis in zebrafish mind bomb alleles')

data_phate <- phate(beta_mat_pcscore)
plot(data_phate$embedding)
pa <- data_phate$embedding
pa <- as.data.frame(pa)
rownames(pa) <- rownames(beta)
geneid <- which(TSNE2$gene %in% HVG,) 
genes <- rownames(TSNE2)[geneid]
pa$X <- rownames(beta)
pa$X[!pa$X %in% genes] <- NA
pa$col <- ifelse(is.na(pa$X), '#2ca25f', 'red')
pa$size <- ifelse(is.na(pa$X),1,2)
ggplot(pa, aes(`PHATE1`, `PHATE2`, label = X)) + geom_point(color = pa$col,alpha=0.5,size=pa$size)+theme_classic()+ggtitle('HVG')

group = '2'
geneid <- which(TSNE2$Group==group,) 
genes <- rownames(TSNE2)[geneid]
pa$X <- rownames(beta)
pa$X[!pa$X %in% genes] <- NA
pa$col <- ifelse(is.na(pa$X), '#2ca25f', 'red')
pa$size <- ifelse(is.na(pa$X),0.1,0.5)
ggplot(pa, aes(`PHATE1`, `PHATE2`, label = X)) + geom_point(color = pa$col,alpha=0.5,size=pa$size)+theme_classic()+ggtitle(paste0('Cluster ',group))

geneid <- which(TSNE2$gene %in% HVG,) 
genes <- rownames(TSNE2)[geneid]
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!TSNE$X %in% genes] <- NA
TSNE$col <- ifelse(is.na(TSNE$X), '#2ca25f', 'red')
TSNE$size <- ifelse(is.na(TSNE$X),1,2)
ggplot(TSNE, aes(`TSNE  1`, `TSNE 2`)) + geom_point(color = TSNE$col,alpha=0.5,size=TSNE$size) +
  theme_classic()+ggtitle('HVG')

