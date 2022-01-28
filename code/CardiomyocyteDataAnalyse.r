library(ggplot2)
library(Matrix)
library(monocle3)
library(Seurat)
library(viridis)
library(statsExpressions)
library(tidyverse)
library(GSVA)
library(fgsea)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(Rtsne)

setwd("Z:/Cailab/scTenifoldTime/")
data <- read.csv('Data/Cardiomyocytes/Cardiomyocytes.csv')
rownames(data) <- make.unique(data$X)
data <- data[,-(1)]
pheno_data <- read.csv('Z:/Cailab/scTenifoldTime/Data/Cardiomyocytes/phenodata.csv')
rownames(pheno_data) <- pheno_data$X
timepoint_levels = c("e14", "e18", "p0", "p1", "p4", "p8", "p11", "p14", "p15", "p18", "p22", "p28", "p35", "p56", "p84")
pheno_data$Timepoint = factor(pheno_data$Timepoint, levels = timepoint_levels)
data = data[, rownames(pheno_data)]

S <- CreateSeuratObject(counts = data)
S$Timepoint <- pheno_data$Timepoint[match(colnames(S),pheno_data$X)]
S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst",nfeatures = 2000)
S <- ScaleData(S, features = rownames(S))
S <-RunPCA(S, features = VariableFeatures(object = S), verbose = FALSE, 
       ndims.print = NULL, nfeatures.print = NULL, reduction.key = "PC_")

S <- FindNeighbors(S, dims = 1:20, verbose = FALSE, k.param = 20, compute.SNN = TRUE, 
                     prune.SNN = 1/15, nn.eps = 0)


S <- FindClusters(S, verbose = FALSE, resolution = 2)
S <- RunUMAP(S, dims = 1:20, verbose = FALSE)
S <- RunTSNE(S,dims = 1:20)
HVG <- VariableFeatures(S)
umap_stage <- DimPlot(S, reduction = "umap", group.by = "Timepoint") +xlab("UMAP 1") + ylab("UMAP 2")
umap_stage

tsne_stage <- DimPlot(S,reduction = 'tsne',group.by = "Timepoint")+xlab("tSNE 1") + ylab("tSNE 2")
tsne_stage 

#early
FeaturePlot(object = S,features = c('Tnni1'),reduction = 'tsne')+scale_color_viridis()
FeaturePlot(object = S,features = c('Myh7'),reduction = 'tsne')+scale_color_viridis()
#late
FeaturePlot(object = S,features = c('Tnni3'),reduction = 'tsne')+scale_color_viridis()
FeaturePlot(object = S,features = c('Myh6'),reduction = 'tsne')+scale_color_viridis()

pheno_data$group = factor(pheno_data$group, levels = c("in vivo"))
G_list <- cbind(rownames(data))
rownames(G_list) <- rownames(data)
colnames(G_list) <- 'gene_short_name'
cds_invivo <- new_cell_data_set(as.sparse(data), 
                                cell_metadata = pheno_data
                                , gene_metadata = G_list)
fData(cds_invivo)$num_cells_expressed = rowSums(exprs(cds_invivo) >= 1) #I turned on multimappers in Kallisto, but I don't feel too bad about considering as "not expressed" if the final count is under 1
cds_invivo <- preprocess_cds(cds_invivo, num_dim = 5)
cds_invivo <- align_cds(cds_invivo, alignment_group = "library")
cds_invivo <- reduce_dimension(cds_invivo, umap.min_dist = 0.2) #Set umap min_dist here just to get a *slightly* neater trajectory (there was an exceedingly mild discontinuity with the default; I don't think setting it makes a huge difference here other than aesthetics)
cds_invivo <- cluster_cells(cds_invivo, k = 100) #The k is set to aggressively ensure there is only one trajectory - sometimes, I think Monocle gets unduly confused because it *expects* multiple trajectories
cds_invivo <- learn_graph(cds_invivo, learn_graph_control = list(ncenter = 300)) #ncenter set to 300 so that there's more range captured in the pseudotime scores

colData(cds_invivo)$root.cells <- cds_invivo@colData$Timepoint %in% c('e14')  & cds_invivo@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_dist"]][1,]> 10
length(which(colData(cds_invivo)$root.cells))
UMAP_root <- plot_cells(cds_invivo, color_cells_by = "root.cells",label_groups_by_cluster = FALSE, 
                 label_cell_groups = FALSE,label_leaves = FALSE, cell_size = 1)
Monocle_stage <- plot_cells(cds_invivo, label_groups_by_cluster = FALSE, 
                 label_cell_groups = FALSE, color_cells_by = "Timepoint",
                 label_branch_points = FALSE,label_leaves = FALSE, cell_size = 1)


get_earliest_principal_node <- function(cds) {
  cell_ids <- which(cds$root.cells)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, 
  ]))))]
  root_pr_nodes
}

root.nodes <- get_earliest_principal_node(cds_invivo)
cds_invivo <- order_cells(cds_invivo, root_pr_nodes = root.nodes)

cds_invivo@colData$pseudotime = cds_invivo@principal_graph_aux$UMAP$pseudotime
UMAP_PT <- plot_cells(cds_invivo,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5, cell_size = 1)+labs(color = "legend title")

markers <- plot_cells(cds_invivo, genes=c("Yy1", "Nfe2l2", "Sox9", "Nrf1",'Mef2a'),show_trajectory_graph = F,cell_size = 1)
UMAP_PT
markers

PT <- cds_invivo@colData$pseudotime
df <- as.data.frame(cbind(names(PT),PT))
# write.csv(df,file = 'Data/Cardiomyocytes/Barcode_pseudotime_r4.csv',row.names = F)
# pheno_data$pseudotime <- df$PT[match(pheno_data$X,df$V1)]
# write.csv(pheno_data,file='Data/Cardiomyocytes/phenodata.csv',row.names = F)

S$pseudotime <- pseudotime(cds_invivo)[match(colnames(S), colnames(cds_invivo))]
umap.pt <- FeaturePlot(S, feature = "pseudotime")
umap.pt <- umap.pt + scale_colour_viridis(option = "plasma")+labs(title = 'Pseudotime')+xlab("UMAP 1") + ylab("UMAP 2")
umap.pt

tsne.pt <- FeaturePlot(S, feature = "pseudotime",reduction = 'tsne')
tsne.pt <- tsne.pt + scale_colour_viridis(option = "plasma")+labs(title = 'Pseudotime')+xlab("tSNE 1") + ylab("tSNE 2")
tsne.pt

# Run CCAT to calculate cell entropy
Ref <- read.csv("Z:/Cailab/scTenifoldTime/Data/Zebrafish/Human_entrz_esembl.Ref.csv")
Ref <- Ref[Ref$entrezIDs!='nooutput',]
Aname <- rownames(cds_invivo)
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
rownames(cds_invivo) <- AID
Amtrix <- cds_invivo@assays@data$counts
rownames(Amtrix) <- rownames(cds_invivo)

load("Z:/Cailab/scTenifoldTime/Data/Zebrafish/net17Jan16.rda")
ppiA.m <- net17Jan16.m
exp.m <- log2(Amtrix+1)
exp.m <- as.matrix(exp.m)
commonEID.v <- intersect(rownames(ppiA.m),rownames(exp.m))
k.v <- rowSums(ppiA.m[match(commonEID.v, rownames(ppiA.m)),])
ccat.v <- as.vector(cor(exp.m[match(commonEID.v,rownames(exp.m)),],k.v))
cell <- colnames(exp.m)
CCAT.celltype <-as.data.frame(cbind(cell,ccat.v))
CCAT.celltype$stage <- pheno_data$timepoint[match(CCAT.celltype$cell,pheno_data$X)]
pheno_data$CCAT <-rank(as.numeric(CCAT.celltype$ccat.v[match(pheno_data$X, CCAT.celltype$cell)]))
S$CCAT <-rank( as.numeric(CCAT.celltype$ccat.v[match(colnames(S), CCAT.celltype$cell)]))

CCAT_tsne <- FeaturePlot(S, feature = "CCAT",reduction = 'tsne')+ scale_colour_viridis(option = "plasma")+xlab("tSNE 1") + ylab("tSNE 2")+labs(title='Differentiation potency')
CCAT_tsne

CCAT_umap <- FeaturePlot(S, feature = "CCAT")+ scale_colour_viridis(option = "plasma")+xlab("UMAP 1") + ylab("UMAP 2")+labs(title='Differentiation potency')
CCAT_umap

pseudotimePlot <- ggplot(pheno_data, aes(x = Timepoint, y = as.numeric(pseudotime), fill =Timepoint)) + 
  geom_violin() +
  xlab("Time Points") +ylab("Pseudotime")+
  theme(legend.position="none") + 
  xlab("")+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_jitter(position = position_jitter(width = 0.15, height = 0.1),size=0.5,alpha=0.2)+NoLegend()
pseudotimePlot
CCATPlot <- ggplot(pheno_data, aes(x = Timepoint, y = CCAT, fill =Timepoint)) + 
  geom_violin() +
  xlab("Time Points") +ylab("CCAT")+
  theme(legend.position="none") + 
  xlab("")+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_jitter(position = position_jitter(width = 0.15, height = 0.1),size=0.5,alpha=0.2)+NoLegend()

CCATPlot


################################# DE by monocle3 ######################################################
#takes 20 mins to run, so I saved the file here
#library(tidyverse)
# gene_fits <- fit_models(cds_invivo, model_formula_str = "~Timepoint")
# fit_coefs <- coefficient_table(gene_fits)
# time_terms <- fit_coefs %>% filter(q_value < 0.05)
# time_terms <- time_terms[,c('gene_short_name','term','p_value','q_value')]
# DE <- unique(time_terms$gene_short_name)
# DE <- DE[DE %in% colnames(beta)]
#write.table(DE,file = 'figures/cardiomyocytes/DE_by_monocle3.csv',  quote = F,row.names = F)
DE <- read.csv('figures/cardiomyocytes/DE_by_monocle3.csv')
DE <- DE$x

################# Downstream analysis starting from beta-matrix #################
setwd("Z:/Cailab/scTenifoldTime/")
library(GSVA)
library(fgsea)
library(ggplot2)
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
res_pcnet <- readRDS("network/Cardiomyocyte/res_pcnet_7000genes.rds")

#res_pcnet <- readRDS("network/Cardiomyocyte/r4/results/networks/res_pcnet.rds")
beta <- res_pcnet$beta_mat
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

# Plot the enrichment plot of the specific pathway
gSet <- 'Pyruvate metabolism and citric acid (TCA) cycle'
pTitle <- 'Pyruvate metabolism'
plotEnrichment(BIOP[[gSet]], U,ticksSize = 0.001) + 
  labs(title = pTitle) +
  theme_bw() + 
  ylab('Enrichment score') +
  xlab('Gene rank') +
  theme(plot.title = element_text(face = 2, size = 20))

# plot the example networks that change across pseudotime
df <- eD[eD$pathway=='Pyruvate metabolism and citric acid (TCA) cycle',]
genes <- unlist(strsplit(df$leadingEdge,split = ';'))
Networks <- list()
for (i in 1:10){
  Net <- res_pcnet$network_list[[i]]
  Net <- Net[toupper(colnames(Net)) %in% genes,toupper(rownames(Net)) %in% genes]
  Networks[[i]] <- Net
}
Network <- as.matrix(cbind(Networks[[1]],Networks[[3]],Networks[[5]],
                           Networks[[7]],Networks[[9]]))
range(Network)
col_fun = colorRamp2(c(-0.01,0.02), c('green',  "red"))
split = rep(c('Time 1','Time 3',"Time 5", "Time 7",'Time 9'), each=length(genes))
Heatmap(Network, name = "Beta", col = col_fun,cluster_rows = T,cluster_columns = FALSE,
        column_split = split,column_title_gp = gpar(fill = c("#fff7ec", "#fee8c8", "#fdd49e","#fdbb84","#fc8d59")),
        column_gap  = unit(3, "mm") )

# plot the gene set in the tsne plot
TSNE.orig <- TSNE
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!toupper(TSNE$X) %in% genes] <- NA
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

############## Group ##############
geneid <- which(TSNE2$Group==42,) 
genes <- rownames(TSNE2)[geneid]
TSNE <- as.data.frame(TSNE.orig)
TSNE$X <- rownames(beta)
colnames(TSNE) <- c('TSNE  1','TSNE 2','X')
TSNE$X[!TSNE$X %in% genes] <- NA
TSNE$col <- ifelse(is.na(TSNE$X), '#2166ac', 'red')
TSNE$size <- ifelse(is.na(TSNE$X),1,1)
ggplot(TSNE, aes(`TSNE  1`, `TSNE 2`, label = X)) + geom_point(color = TSNE$col,alpha=0.3,size=TSNE$size) +
 theme_classic()+ggtitle('Cluster 42')
genes <- genes[!grepl('Gm|Rik', genes)]

############## Ppara ##############
gene <- 'Ppara'
plot_cells(cds_invivo, genes=gene,show_trajectory_graph = F,cell_size = 1)

vim <- beta[rownames(beta)==gene,]
vim <- vim[!grepl(pattern = 'Rps|Rpl|Rp[[:digit:]]+|Hb|Gm|mt|Atp',names(vim))]
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

range(Network)

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
                   segment.color = 'grey50',max.iter = 50000,max.overlaps = 100)+theme_classic()

genes_PT <- c('Acadvl','Acaa2','Ndufa10','Ppara','Cox6c')
cds_subset <- cds_invivo[rowData(cds_invivo)$gene_short_name %in% genes_PT,]
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="pseudotime",
                         min_expr=0)

cds_subset <- cds_invivo[rowData(cds_invivo)$gene_short_name %in% c('Acadm','Aco2','Eef1a1','Dpysl3','Gja1',genename),]
gene_fits <- fit_models(cds_subset,model_formula_str = '~Timepoint')
fit_coefs <- coefficient_table(gene_fits)

plot_genes_violin(cds_subset, group_cells_by="Timepoint", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="Timepoint")



