# load packages
library(SPATA2,lib.loc = '/cluster/home/yzy_jh/sbin/R/library/R_TMP/')
library(ggplot2)
library(monocle3)
library(magrittr)
library(tidyverse)
library(patchwork)

project_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/SPATA2/'
setwd(project_dir)
CNV_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/SPATA2/CNV'
checkdir(CNV_dir)
monocle_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/SPATA2/monocle'
checkdir(monocle_dir)
phenotype_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/SPATA2/phenotype'
checkdir(phenotype_dir)
fig_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/SPATA2/fig'
checkdir(fig_dir)

#load data
samples <- c('shcc1C','shcc1T','shcc2T','shcc3C','shcc3T','shcc4C','shcc4T','shcc5T')

spata_obj <- parallel::mclapply(samples,function(x){
  spata_obj <-
    initiateSpataObject_10X(
      directory_10X = as.character(glue("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/counts/{x}/")), # the directory from which to load the data
      sample_name = x)
  return(spata_obj)
  
},mc.cores = 8)

names(spata_obj) <- samples
saveRDS(spata_obj,file = "spata_obj.rds")

spata_obj <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/SPATA2/spata_obj.rds')
#CNV
spata_cnv_obj <- parallel::mclapply(samples,function(x){
                  sample_CNV_diR <- glue("{CNV_dir}/{x}")
                  checkdir(sample_CNV_diR)
                  spata_obj[[x]] <-
                    runCnvAnalysis(
                      object = spata_obj[[x]],
                      directory_cnv_folder = sample_CNV_diR ,
                      cnv_prefix = "Chr",
                      gene_pos_df = SPATA2::gene_pos_df
                    )
                  return(spata_obj[[x]])
},mc.cores = 8)
names(spata_cnv_obj) <- samples
write_rds(spata_cnv_obj,file = 'CNV/spata_cnv_obj.rds')
#monocle
spata_monocle_obj <- parallel::mclapply(samples,function(x){
                                    sample_monocle_diR <- glue("{monocle_dir}/{x}") %>% as.character()
                                    checkdir(sample_monocle_diR)
                                    spata_obj[[x]] <-
                                      runCnvAnalysis(
                                        object = spata_obj[[x]],
                                        directory_cnv_folder = sample_CNV_diR ,
                                        cnv_prefix = "Chr",
                                        gene_pos_df = SPATA2::gene_pos_df
                                      )
                                    return(spata_obj[[x]])
},mc.cores = 8)

p0 <- plotSurface(object = spata_cnv_obj$shcc1T, color_by = "seurat_clusters")

shcc1C_cnv <- plotCnvResults(object = spata_cnv_obj$shcc1T)
shcc1T_cnv <- plotCnvResults(object = spata_cnv_obj$shcc1C)
shcc1C_cnv/shcc1T_cnv

shcc3C_cnv <- plotCnvResults(object = spata_cnv_obj$shcc3C)
shcc3T_cnv <- plotCnvResults(object = spata_cnv_obj$shcc3T)
shcc3C_cnv/shcc3T_cnv

shcc4C_cnv <- plotCnvResults(object = spata_cnv_obj$shcc4C)
shcc4T_cnv <- plotCnvResults(object = spata_cnv_obj$shcc4T)
shcc4C_cnv/shcc4T_cnv
plotSurface(object = spata_cnv_obj$shcc4C, color_by = "Chr17", pt_clrsp = "Reds 3", c1 = 1)
plotSurface(object = spata_cnv_obj$shcc4T, color_by = "Chr17", pt_clrsp = "Reds 3", c1 = 1)

shcc5T_cnv <- plotCnvResults(object = spata_cnv_obj$shcc5T)

plotSurface(object = spata_obj$shcc5T, color_by = "seurat_clusters")
plotCnvResults(object = spata_cnv_obj$shcc5T, across = "seurat_clusters")
plotSurface(object = spata_cnv_obj$shcc5T, color_by = "Chr7", pt_clrsp = "Reds 3", c1 = 1)

# compile a cell-data-set
cortex_cds <- transformSpataToCDS(object = spata_obj$shcc5T)
# the pseudotime values for all cells/barcode-spots are obtained via
cortex_cds <- order_cells(cortex_cds)
pseudotime_vec <- monocle3::pseudotime(cortex_cds)

# subset output
pseudotime_vec[1:10]
# 1. convert to data.frame 
pseudotime_df <- 
  base::as.data.frame(pseudotime_vec)

# output
head(pseudotime_df)
# 2. create a barcodes-key variable and rename the pseudotime variable 
feature_df <- 
  magrittr::set_colnames(x = pseudotime_df, value = "Pseudotime") %>% 
  tibble::rownames_to_column(var = "barcodes")

# 3. add to spata-object via 
spata_obj$shcc5T <- addFeatures(object = spata_obj$shcc5T,
                         feature_names = "Pseudotime", 
                         feature_df = feature_df, 
                         overwrite = TRUE)
# visualize pseudotime on the surface
plotSurface(spata_obj$shcc5T, color_by = "Pseudotime", smooth = TRUE, smooth_span = 0.2, pt_size = 1.8) # output
# open interactive application
spata_obj$shcc5T <- createTrajectories(object = spata_obj$shcc5T)
# get trajectory names 
getTrajectoryNames(object = spata_obj$shcc5T)
plotTrajectory(object = spata_obj$shcc5T, 
               trajectory_name = "tumor_areas",
               color_by = "seurat_clusters",
               pt_clrp = "npg",
               pt_alpha = 0.25, # reduce alpha to highlight the trajectory's course
               display_image = FALSE) +
  legendTop()

plotTrajectory(object = spata_obj$shcc5T, 
               trajectory_name = "tumor_areas",
               color_by = "nCount_Spatial",
               smooth = TRUE, 
               pt_alpha = 0.25, 
               display_image = FALSE) +
  legendTop()
plotTrajectoryFeatures(object = spata_obj$shcc5T,
                       trajectory_name = "tumor_areas",
                       features = "nCount_Spatial", 
                       smooth_method = "loess", 
                       smooth_span = 0.2, 
                       smooth_se = TRUE) 

plotTrajectoryFeaturesDiscrete(object = spata_obj$shcc5T,
                               trajectory_name = "tumor_areas",
                               discrete_feature = "seurat_clusters", 
                               clrp = "npg",
                               display_trajectory_parts = FALSE) 
all_genes <- getGenes(spata_obj$shcc5T)
all_gene_sets <- getGeneSets(spata_obj$shcc5T)

# obtain an assessed trajectory data.frame for all genes
atdf_genes <- assessTrajectoryTrends(object = spata_obj$shcc5T, 
                                     trajectory_name = "tumor_areas", 
                                     variables = all_genes)

# obtain an assessed trajectory data.frame for all gene-sets
atdf_gene_sets <- assessTrajectoryTrends(object = spata_obj$shcc5T, 
                                         trajectory_name = "tumor_areas", 
                                         variables = '')

# output example
atdf_genes

spata_obj <- 
  runDeAnalysis(object = spata_obj$shcc5T,
                across = "seurat_clusters", # across which identity groups
                method_de = "wilcox" # with which methods
  )

spata_obj_gsea <- 
  runGsea(
    object = spata_obj$shcc5T, 
    across = "seurat_clusters",
    methods_de = "wilcox" 
  )

# compare the trend of a variable to different models
plotTrajectoryFit(object = spata_obj$shcc5T,
                  trajectory_name = "tumor_areas", 
                  variable = c("APOC4",'FGG','FGB','FTH1','SERPINC1','SAA2','GPX2'), 
                  display_residuals = TRUE) + 
  legendTop()

# example 1: extract all variables that follow the linear descending trend while moving towards the hypoxic area (See Figure 2.1)
# with an auc-evaluation equal to or lower than 2
descending_genes <-
  filterTrajectoryTrends(atdf = atdf_genes,
                         limit = 3,
                         trends = "Gradient descending", 
                         variables_only = FALSE) # return a data.frame

descending_genes
descending_genes_vec <- descending_genes$variables
hm_colors <- viridis::inferno(n = 100)
plotTrajectoryHeatmap(object = spata_obj$shcc5T, 
                      trajectory_name = "tumor_areas",
                      variables = descending_genes_vec,
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_rownames = TRUE,
                      split_columns = FALSE, 
                      smooth_span = 0.5)

peaking_genes <-
  filterTrajectoryTrends(atdf = atdf_genes,
                         limit = 6,
                         trends = c("Early peak", "Late peak", "One peak"),
                         variables_only = TRUE) # return a vector of variables
plotTrajectoryHeatmap(object = spata_obj$shcc5T, 
                      trajectory_name = "tumor_areas",
                      variables = peaking_genes, 
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_row_names = FALSE, 
                      split_columns = TRUE, # splits the heatmap to highlight the trajectory parts
                      smooth_span = 0.5)

# example 2: extract variables that featured a peak along the trajectory 
# with auc-evaluation equal to or lower than 4


head(peaking_genes)

trajectory_length <- getTrajectoryLength(spata_obj$shcc5T, trajectory_name = "tumor_areas", binwidth = 5)

trajectory_length

trajectory_direction <- 1:trajectory_length

linear_ascending <- scales::rescale(1:trajectory_length, to = c(0,1))

plotTrajectoryGenes(object = spata_obj$shcc5T, trajectory_name = "tumor_areas", genes = "FTH1") + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = trajectory_direction, y = linear_ascending),
                     color = "blue",
                     size = 1)

traj_df <- getTrajectoryDf(object = spata_obj$shcc5T,
                           trajectory_name = "tumor_areas",
                           variables = c("HM_HYPOXIA","",# hypoxia gene set
                                         "FTH1"), # METRN gene 
                           binwidth = 5, 
                           shift_wider = TRUE)

dplyr::select(traj_df, trajectory_part, trajectory_order, HM_HYPOXIA) 


similar_to_HM_HYPOXIA <- filterTrajectoryTrends(atdf_cust, limit = 0.6, trends = "HM_HYPOXIA")

# print names of similar genes
similar_to_HM_HYPOXIA
# plot both in comparison 
plotTrajectoryGeneSets(object = spata_obj$shcc5T, 
                       trajectory_name = "tumor_areas", 
                       gene_sets = c("HM_HYPOXIA"), 
                       display_facets = TRUE, 
                       clrp = "default")

plotTrajectoryGenes(object = spata_obj$shcc5T, 
                    trajectory_name = "tumor_areas", 
                    genes = similar_to_HM_HYPOXIA[1:4], 
                    display_facets = TRUE)


# gene-set names
genes_of_interest <- c("CALM1", "VIM", "GFAP", "METRN")

# plot lineplot
plotTrajectoryGenes(object = spata_obj,
                    trajectory_name = "tumor-areas", 
                    genes = genes_of_interest,
                    smooth_span = 0.2,
                    smooth_se = TRUE, 
                    display_facets = TRUE, # use facet_wrap() to split the plot in four parts
                    nrow = 2 # align the sub plots in two rows 
)

plotTrajectoryGeneSets(
  object = spata_obj,
  trajectory_name = "tumor-areas",
  gene_sets = "HM_HYPOXIA",
  display_trajectory_parts = FALSE) + # results in missing vertical lines 
  legendTop()

all_genes <- getGenes(spata_obj$shcc5T)
all_gene_sets <- getGeneSets(spata_obj$shcc5T)

# obtain an assessed trajectory data.frame for all genes
atdf_genes <- assessTrajectoryTrends(object = spata_obj$shcc5T, 
                                     trajectory_name = "tumor_areas", 
                                     variables = all_genes)

# obtain an assessed trajectory data.frame for all gene-sets
atdf_gene_sets <- assessTrajectoryTrends(object = spata_obj$shcc5T, 
                                         trajectory_name = "tumor_areas", 
                                         variables = all_gene_sets)

# output example
atdf_genes
# compare the trend of a variable to different models
plotTrajectoryFit(object = spata_obj$shcc5T,
                  trajectory_name = "tumor_areas", 
                  variable = "GFAP", 
                  display_residuals = TRUE) + 
  legendTop()

# example 2: extract variables that featured a peak along the trajectory 
# with auc-evaluation equal to or lower than 4
peaking_genes <-
  filterTrajectoryTrends(atdf = atdf_genes,
                         limit = 4,
                         trends = c("Early peak", "Late peak", "One peak"),
                         variables_only = TRUE) # return a vector of variables

head(peaking_genes)
trajectory_length <- getTrajectoryLength(spata_obj$shcc5T, trajectory_name = "tumor_areas", binwidth = 5)

trajectory_length
trajectory_direction <- 1:trajectory_length

linear_ascending <- scales::rescale(1:trajectory_length, to = c(0,1))

plotTrajectoryGenes(object = spata_obj$shcc5T, trajectory_name = "tumor_areas", genes = "APOA5") + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = trajectory_direction, y = linear_ascending),
                     color = "blue",
                     size = 1)
traj_df <- getTrajectoryDf(object = spata_obj$shcc5T,
                           trajectory_name = "tumor_areas",
                           variables = c("HM_HYPOXIA", # hypoxia gene set
                                         "METRN"), # METRN gene 
                           binwidth = 5, 
                           shift_wider = TRUE)

dplyr::select(traj_df, trajectory_part, trajectory_order, HM_HYPOXIA) 
atdf_cust <- assessTrajectoryTrendsCustomized(object = spata_obj$shcc5T,
                                              trajectory_name = "tumor_areas",
                                              customized_trends = dplyr::select(traj_df, HM_HYPOXIA, METRN), 
                                              variables = all_genes
)

atdf_cust
similar_to_HM_HYPOXIA <- filterTrajectoryTrends(atdf_cust, limit = 0.8, trends = "HM_HYPOXIA")

# print names of similar genes
similar_to_HM_HYPOXIA
# plot both in comparison 
plotTrajectoryGeneSets(object = spata_obj$shcc5T, 
                       trajectory_name = "tumor_areas", 
                       gene_sets = c("HM_HYPOXIA"), 
                       display_facets = TRUE, 
                       clrp = "default")

plotTrajectoryGenes(object = spata_obj$shcc5T, 
                    trajectory_name = "tumor_areas", 
                    genes = similar_to_HM_HYPOXIA[1:4], 
                    display_facets = TRUE)

similar_to_METRN <- filterTrajectoryTrends(atdf_cust, limit = 2, trends = "METRN")

# print names of similar gene sets
similar_to_METRN

# plot both in comparison 
plotTrajectoryGenes(object = spata_obj$shcc5T, 
                    trajectory_name = "tumor_areas", 
                    genes = c("METRN"), 
                    display_facets = TRUE, 
                    clrp = "default")

plotTrajectoryGenes(object = spata_obj$shcc5T, 
                    trajectory_name = "tumor_areas", 
                    genes = similar_to_METRN, 
                    display_facets = TRUE)

descending_genes_vec <- descending_genes$variables

hm_colors <- viridis::inferno(n = 100)

plotTrajectoryHeatmap(object = spata_obj$shcc5T, 
                      trajectory_name = "tumor_areas",
                      variables = descending_genes_vec,
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_rownames = TRUE,
                      split_columns = FALSE, 
                      smooth_span = 0.5)

plotTrajectoryHeatmap(object = spata_obj$shcc5T, 
                      trajectory_name = "tumor_areas",
                      variables = peaking_genes, 
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_row_names = FALSE, 
                      split_columns = TRUE, # splits the heatmap to highlight the trajectory parts
                      smooth_span = 0.5)



spata_obj$shcc5TrunDea(spata_obj$shcc5T)

getFeatureNames(spata_obj$shcc5T)
monocle_clusters <- findMonocleClusters(object = spata_obj$shcc5T, 
                                        preprocess_method = "PCA", 
                                        reduction_method = c("UMAP", "PCA", "tSNE"), 
                                        cluster_method = c("leiden", "louvain"), 
                                        k = 5, 
                                        num_iter = 5)

# output
monocle_clusters
# add the cluster results
spata_obj$shcc5T <- 
  addFeatures(object = spata_obj$shcc5T, 
              feature_names = c("cluster_leiden_UMAP_k5", "cluster_leiden_tSNE_k5","cluster_louvain_PCA_k5"), 
              feature_df = monocle_clusters,
              overwrite = TRUE,
              key = "barcodes")
spata_obj$shcc5T <- renameFeatures(object = spata_obj$shcc5T, "Leiden_UMAP" = "cluster_leiden_UMAP_k5")
# feature names afterwards
getFeatureNames(spata_obj$shcc5T)
spata_obj$shcc5T <- 
  runDeAnalysis(object = spata_obj$shcc5T,
                across = "Leiden_UMAP", # across which identity groups
                method_de = c("wilcox", "bimod") # with which methods
  )

# get an overview about the de-analysis results stored in your object
printDeaOverview(spata_obj$shcc5T)
de <- getDeaResultsDf(object = spata_obj$shcc5T, 
                across = "Leiden_UMAP", 
                method_de = "wilcox")

spata_obj$shcc5T <- 
  runGsea(
    object = spata_obj$shcc5T, 
    across = "seurat_clusters",
    methods_de = "wilcox" 
  )
getGseaDf(
  object = spata_obj$shcc5T, 
  across = "seurat_clusters",
  method_de = "wilcox", 
  n_gsets = 20 # extract top 20 most significant gene sets
) 
plotGseaDotPlot(
  object = spata_obj$shcc5T,
  across = "seurat_clusters",
  across_subset = "3",
  n_gsets = 15,
  transform_with = list("fdr" = "log10"), 
  by_group = TRUE
)
plotGseaDotPlot(
  object = spata_obj$shcc5T,
  across = "seurat_clusters",
  across_subset = c("1", "2", "4", "5", "6"),
  n_gsets = 7,
  pt_alpha = 0.8,
  transform_with = list("fdr" = c("log10")),
  by_group = FALSE # merge in one plot
) 




head(feature_df)
monocle_plot <- 
  plot_cells(cortex_cds,
             reduction_method = "UMAP",
             color_cells_by = "seurat_clusters" # visualize a transferred variable
  ) +
  scale_color_add_on(aes = "color", variable = "discrete", clrp = "npg") +
  legendNone() 
monocle_plot
# 3. add to spata-object via 
spata_obj$shcc5T <- addFeatures(object = spata_obj$shcc5T,
                         feature_names = "Pseudotime", 
                         feature_df = feature_df, 
                         overwrite = TRUE)


plotSurface(object = spata_obj$shcc1C, 
            color_by = "Leiden_UMAP", 
            pt_clrp = "jama", 
            pt_size = 1.9
)
getFeatureNames(spata_obj$shcc5T)
monocle_clusters <- findMonocleClusters(object = spata_obj$shcc1C, 
                                        preprocess_method = "PCA", 
                                        reduction_method = c("UMAP", "PCA", "tSNE"), 
                                        cluster_method = c("leiden", "louvain"), 
                                        k = 5, 
                                        num_iter = 5)

# get an overview about the de-analysis results stored in your object
printDeaOverview(spata_obj$shcc5T)
# open interactive application
seurat_obj_shcc5T <- transformSpataToSeurat(object = spata_obj$shcc5T)
DefaultAssay(seurat_obj_shcc5T)
p0 <- Seurat::SpatialDimPlot(object = seurat_obj_shcc5T,alpha = 0) + NoLegend()
p1 <- Seurat::SpatialDimPlot(object = seurat_obj_shcc5T,label = T)
p2 <- Seurat::SpatialFeaturePlot(seurat_obj_shcc5T,features = 'HYPOXIA')
p3 <- Seurat::SpatialFeaturePlot(seurat_obj_shcc5T,features = 'Pseudotime')
p_m <- p0|p1|p2|p3


P3 <- VlnPlot(seurat_obj_shcc5T,features = 'HYPOXIA',pt.size = 0)
P4 <- VlnPlot(seurat_obj_shcc5T,features = 'Pseudotime',pt.size = 0)
p_n <- P3|P4
fig_shcc5 <- p_m/p_n
ggsave(fig_shcc5,filename = glue("{fig_dir}/fig_shcc5.png"),width = 8,height = 6)


seurat_obj_shcc5T_cnv <- transformSpataToSeurat(object = spata_cnv_obj$shcc5T)
p1 <- Seurat::SpatialFeaturePlot(object = seurat_obj_shcc5T_cnv, features = "Chr4")
p2 <-Seurat::SpatialFeaturePlot(object = seurat_obj_shcc5T_cnv, features = "Chr15")
p3 <- Seurat::SpatialFeaturePlot(object = seurat_obj_shcc5T_cnv, features = "Chr17")
p4 <- Seurat::SpatialFeaturePlot(object = seurat_obj_shcc5T_cnv, features = "Chr7")
p5 <- Seurat::SpatialFeaturePlot(object = seurat_obj_shcc5T_cnv, features = "Chr23")
p_right <- p1|p2|p3
p_left <- p4|p5
ggsave(p_right/p_left,filename = '../../visium/Integrative/SPATA2/CNV/shcc5T/CNV.png',width = 10,height = 8)
p_cnv_cluster <- plotCnvResults(object = spata_cnv_obj$shcc5T, across = "seurat_clusters") + 
  theme(axis.title.x = element_text(size = 5))
ggsave(p_cnv_cluster,filename = '../../visium/Integrative/SPATA2/CNV/shcc5T/p_cnv_cluster.png',width = 12,height = 8)

p_cnv <- p1|p2
p_cnv

spata_obj$shcc5T <- 
  runDeAnalysis(object = spata_obj$shcc5T,
                across = "seurat_clusters", # across which identity groups
                method_de ="wilcox"# with which methods
  )
de <- getDeaResultsDf(object = spata_obj$shcc5T, 
                across = "seurat_clusters", 
                method_de = "wilcox", 
                max_adj_pval = 0.025, # of every cluster take genes with an adj. p-value of 0.025 or lower
                n_lowest_pval = 20)

cluster_of_interest <- c("3", "4", "5","7")

plotDeaHeatmap(object = spata_obj$shcc5T, 
               across = "seurat_clusters", # the grouping variable
               across_subset = cluster_of_interest, # the identity groups of interest 
               method_de = "wilcox", # the method with which the results were computed
               max_adj_pval = 0.025, # the adjusted p-value threshold
               n_lowest_pval = 20, 
               n_highest_lfc = 20, 
               clrp = "jama", 
               fontsize = 8) 

marker_sp5 <- FindAllMarkers(seurat_obj_shcc5T)
marker_sp5_top50 <- de %>% group_by(seurat_clusters) %>% top_n(30,wt = avg_logFC)
ids = bitr(marker_sp5_top50$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
markers_bitr = merge(marker_sp5_top50,ids,by.x='gene',by.y='SYMBOL')
markers_bitr <- split(markers_bitr$ENTREZID,marker_sp5_top50$seurat_clusters)

KEGG_cluster <- compareCluster(markers_bitr, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
p_KEGG <- dotplot(KEGG_cluster) + theme(axis.text.x = element_text(angle = 45, 
                                        vjust = 0.5, hjust=0.5,size=15),
                                        axis.text.y = element_text(size=15)) + 
                                  ggtitle("KEGG")

GO_cluster_BP <- compareCluster(markers_bitr, fun="enrichGO",
                     OrgDb =org.Hs.eg.db, ont = "BP",pvalueCutoff=0.05) 
GO_cluster_CC <- compareCluster(markers_bitr, fun="enrichGO",
                             OrgDb =org.Hs.eg.db, ont = "CC",pvalueCutoff=0.05)
GO_cluster_MF <- compareCluster(markers_bitr, fun="enrichGO",
                             OrgDb =org.Hs.eg.db, ont = "MF",pvalueCutoff=0.05)

p_BP <- dotplot(GO_cluster_BP) + theme(axis.text.x = element_text(angle = 45, 
                                       vjust = 0.5, hjust=0.5,size=15),
                                       axis.text.y = element_text(size=15)) + ggtitle("GO-BP")
p_CC <- dotplot(GO_cluster_CC) + theme(axis.text.x = element_text(angle = 45, 
                                       vjust = 0.5, hjust=0.5,size=15),
                                       axis.text.y = element_text(size=15)) + ggtitle("GO-CC")
p_MF <- dotplot(GO_cluster_MF) + theme(axis.text.x = element_text(angle = 45, 
                                       vjust = 0.5, hjust=0.5,size=15),
                                       axis.text.y = element_text(size=15)) + ggtitle("GO-MF")
p_enrich <- p_KEGG|p_BP|p_MF
p_enrich
ggsave(p_KEGG,filename = glue("{fig_dir}/p_KEGG.png"),width = 10,height = 15)
ggsave(p_BP,filename = glue("{fig_dir}/p_BP.png"),width = 10,height = 15)
ggsave(p_CC,filename = glue("{fig_dir}/p_CC.png"),width = 10,height = 15)
ggsave(p_MF,filename = glue("{fig_dir}/p_MF.png"),width = 15,height = 25)
ggsave(p_enrich,filename = glue("{fig_dir}/p_enrich.png"),width = 35,height = 30)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

p1 <- SpatialFeaturePlot(seurat_obj,features ='G2M.Score')
p2 <- SpatialFeaturePlot(tmp,features ='G2M.Score')
p1|p2

#fig1 cluster 
parallel::mclapply(samples,function(x){
  p <- plotSurface(spata_obj[[x]],
                   color_by = "HYPOXIA",
                   pt_size = 1.9,
                   smooth = TRUE, 
                   smooth_span = 0.5,
                   pt_clrsp = "magma") 
  p1 <- plotSurface(spata_obj[[x]],
                    color_by = "nCount_Spatial",
                    pt_size = 1.9,
                    smooth = TRUE, 
                    smooth_span = 0.5,
                    pt_clrsp = "BuPu") 
  p2 <- plotSurface(spata_obj[[x]],
                    color_by = "nFeature_Spatial",
                    pt_size = 1.9,
                    smooth = TRUE, 
                    smooth_span = 0.5,
                    pt_clrsp = "viridis") 
  p3 <-   plotSurface(object = spata_obj[[x]],
                      color_by = "seurat_clusters",
                      pt_size = 1.5) +
    ggplot2::labs(color = "Clusters")
  ggsave(p|p1|p2,filename = glue("{phenotype_dir}/{x}_HYPOXIA.png"),width = 15,height = 8)
  ggsave(p3,filename = glue("{phenotype_dir}/{x}_cluster.png"),width = 8,height = 8)
},mc.cores = 8)

#



# the default
plotSurface(object = spata_obj$shcc5T, color_by = "seurat_clusters") 

# overwrite default by specifying the argument within the function call
plotSurface(object = spata_obj$shcc5T, color_by = "seurat_clusters", pt_clrp = "jco")

################################### CNV ########################################
#monocle3

# compile a cell-data-set
sp_cds <- transformSpataToCDS(object = spata_obj)

pseudotime_vec <- monocle3::pseudotime(sp_cds)
# 1. convert to data.frame 
pseudotime_df <- 
  base::as.data.frame(pseudotime_vec)

# output
head(pseudotime_df)
# 2. create a barcodes-key variable and rename the pseudotime variable 
feature_df <- 
  magrittr::set_colnames(x = pseudotime_df, value = "Pseudotime") %>% 
  tibble::rownames_to_column(var = "barcodes")

# output
head(feature_df)

# 3. add to spata-object via 
spata_obj <- addFeatures(object = spata_obj,
                         feature_names = "Pseudotime", 
                         feature_df = feature_df, 
                         overwrite = TRUE)
# visualize pseudotime on the surface
plotSurface(spata_obj, color_by = "Pseudotime", smooth = TRUE, smooth_span = 0.2, pt_size = 1.8) 

spata_plot <- 
  plotSurface(spata_obj, color_by = "seurat_clusters", pt_clrp = "npg", pt_size = 1.4) + 
  legendBottom()
# SPATAS transform functions pass the feature data to the compiled object which makes 
# all features computed by yourself available, e.g. for monocle3::plot_cells()
monocle_plot <- 
  plot_cells(sp_cds,
             reduction_method = "UMAP",
             color_cells_by = "seurat_clusters" # visualize a transferred variable
  ) +
  scale_color_add_on(aes = "color", variable = "discrete", clrp = "npg") 

# output 1 & 2 
spata_plot|monocle_plot
# example 1: extract all variables that follow the linear descending trend while moving towards the hypoxic area (See Figure 2.1)
# with an auc-evaluation equal to or lower than 2
descending_genes <-
  filterTrajectoryTrends(atdf = atdf_genes,
                         limit = 2,
                         trends = "Gradient descending", 
                         variables_only = FALSE) # return a data.frame

descending_genes

# example 2: extract variables that featured a peak along the trajectory 
# with auc-evaluation equal to or lower than 4
peaking_genes <-
  filterTrajectoryTrends(atdf = atdf_genes,
                         limit = 4,
                         trends = c("Early peak", "Late peak", "One peak"),
                         variables_only = TRUE) # return a vector of variables

head(peaking_genes)

trajectory_length <- getTrajectoryLength(spata_obj, trajectory_name = "tumor-areas", binwidth = 5)

trajectory_length

trajectory_direction <- 1:trajectory_length

linear_ascending <- scales::rescale(1:trajectory_length, to = c(0,1))

plotTrajectoryGenes(object = spata_obj, trajectory_name = "tumor-areas", genes = "GFAP") + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = trajectory_direction, y = linear_ascending),
                     color = "blue",
                     size = 1)

traj_df <- getTrajectoryDf(object = spata_obj,
                           trajectory_name = "tumor-areas",
                           variables = c("HM_HYPOXIA", # hypoxia gene set
                                         "METRN"), # METRN gene 
                           binwidth = 5, 
                           shift_wider = TRUE)

dplyr::select(traj_df, trajectory_part, trajectory_order, HM_HYPOXIA) 

atdf_cust <- assessTrajectoryTrendsCustomized(object = spata_obj,
                                              trajectory_name = "tumor-areas",
                                              customized_trends = dplyr::select(traj_df, HM_HYPOXIA, METRN), 
                                              variables = all_genes
)

atdf_cust

similar_to_HM_HYPOXIA <- filterTrajectoryTrends(atdf_cust, limit = 0.8, trends = "HM_HYPOXIA")

# print names of similar genes
similar_to_HM_HYPOXIA

# plot both in comparison 
plotTrajectoryGeneSets(object = spata_obj, 
                       trajectory_name = "tumor-areas", 
                       gene_sets = c("HM_HYPOXIA"), 
                       display_facets = TRUE, 
                       clrp = "default")

plotTrajectoryGenes(object = spata_obj, 
                    trajectory_name = "tumor-areas", 
                    genes = similar_to_HM_HYPOXIA[1:4], 
                    display_facets = TRUE)

similar_to_METRN <- filterTrajectoryTrends(atdf_cust, limit = 2, trends = "METRN")

# print names of similar gene sets
similar_to_METRN

# plot both in comparison 
plotTrajectoryGenes(object = spata_obj, 
                    trajectory_name = "tumor-areas", 
                    genes = c("METRN"), 
                    display_facets = TRUE, 
                    clrp = "default")

plotTrajectoryGenes(object = spata_obj, 
                    trajectory_name = "tumor-areas", 
                    genes = similar_to_METRN, 
                    display_facets = TRUE)

descending_genes_vec <- descending_genes$variables

hm_colors <- viridis::inferno(n = 100)

plotTrajectoryHeatmap(object = spata_obj, 
                      trajectory_name = "tumor-areas",
                      variables = descending_genes_vec,
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_rownames = TRUE,
                      split_columns = FALSE, 
                      smooth_span = 0.5)

plotTrajectoryHeatmap(object = spata_obj, 
                      trajectory_name = "tumor-areas",
                      variables = peaking_genes, 
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_row_names = FALSE, 
                      split_columns = TRUE, # splits the heatmap to highlight the trajectory parts
                      smooth_span = 0.5)



hm_colors <- viridis::inferno(n = 100)

plotTrajectoryHeatmap(object = spata_obj, 
                      trajectory_name = "tumor-areas",
                      variables = descending_genes_vec,
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_rownames = TRUE,
                      split_columns = FALSE, 
                      smooth_span = 0.5)

plotTrajectoryHeatmap(object = spata_obj, 
                      trajectory_name = "tumor-areas",
                      variables = peaking_genes, 
                      arrange_rows = "maxima",
                      colors = hm_colors,
                      show_row_names = FALSE, 
                      split_columns = TRUE, # splits the heatmap to highlight the trajectory parts
                      smooth_span = 0.5)

# visualize pseudotime on the surface
sp_seurat <- transformSpataToSeurat(object = spata_obj)
Seurat::SpatialFeaturePlot(object = sp_seurat, features = "Pseudotime")

#Autoencoder Denoising
# all expression matrices before denoising
getExpressionMatrixNames(object = spata_obj)

# active expression matrix before denoising
getActiveMatrixName(object = spata_obj)

# denoising your data 
spata_obj <-
  runAutoencoderDenoising(
    object = spata_obj, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )

# all expression matrices after denoising
getExpressionMatrixNames(object = spata_obj)

# active expression matrix after denoising
getActiveMatrixName(object = spata_obj)















# store example genes of interest as character vectors
genes_a <- c("TUBA1B", "HOPX", "PLP1", "ACTB")
genes_b <- c("CARTPT", "OLIG1", "GFAP", "SYNPR", "HOPX", "CCK")
# plot a cluster feature
p1 <- 
  plotSurface(object = spata_obj,
              color_by = "seurat_clusters",
              pt_clrp = "npg",
              pt_size = 1.2) +
  ggplot2::labs(color = "Clusters")

# plot gene expression 
p2 <- 
  plotSurface(object = spata_obj, 
              color_by = "TUBA1B",
              pt_size = 1.8,
              pt_clrsp = "magma")

# combine with patchwork 
p1 + legendTop() +
  p2 + legendTop() 
# plot gene expression 
p1 <- 
  plotSurface(object = spata_obj, 
              color_by = "TUBA1B",
              pt_size = 1.8,
              pt_clrsp = "magma"
  )

# plot gene expression (spatially smoothed)
p2 <- 
  plotSurface(object = spata_obj, 
              color_by = "TUBA1B",
              pt_size = 1.8,
              pt_clrsp = "magma",
              smooth = TRUE, 
              smooth_span = 0.5)

# combine with patchwork 
p1 + legendNone() +
  p2 + legendTop()

plots <- plotSurfaceInteractive(object = spata_obj)

input_list <- list(Bcell = c('CD19', 'CD79A', 'MS4A1', 'CD22'),
                   Tcell = c('CD3D', 'CD3E', 'CD3G','IL7R')       ,   
                   NKcell = c('NKG7', 'FCGR3A', 'IFNG', 'GNLY', 'GZMB','KLRD1','KLRF1') ,          
                   HSC_MFB =  c('RGS5', 'COL1A1','ACTA2','PDGFRB')   ,
                   Endothelialcell =   c('PECAM1','CDH5','STC1','TM4SF1'),          
                   Monocyte = c('CST3', 'LYZ','CD1C', 'FCER1A', 'FCER1G'),
                   Macro_Kupffercells  = c('CST3', 'LYZ', 'CD68', 'CD163', 'FCER1G'),
                   Hepatocytes = c('APOA2','ALB','APOC1','APOC3'),                                
                   DC = c('CST3','LYZ', 'CD1C','IDO1'),
                   Cycling = c('STMN1','TOP2A','NUSAP1','HMGB2'))

plotSurfaceAverage(object = spata_obj, 
                   color_by = input_list, 
                   smooth = TRUE,
                   pt_size = 1.5,
                   pt_clrsp = "viridis")
plotSurface(object = spata_obj, color_by = "seurat_clusters")


# pca 
spata_obj <- runPca(object = spata_obj)

# tsne 
spata_obj <- runTsne(object = spata_obj, n_pcs = 20) 

# umap 
spata_obj <- runUmap(object = spata_obj, n_pcs = 20)
# pca
getPcaDf(object = spata_obj, n_pcs = 20)
# tsne
getTsneDf(object = spata_obj)
# umap
getUmapDf(object = spata_obj)
spata_df <- 
  getSpataDf(object = spata_obj)

with_gfap <- joinWith(object = spata_obj,
                      spata_df = spata_df,
                      genes = "GFAP")

# output 
with_gfap
with_pca <- joinWithPca(object = spata_obj,
                        spata_df = with_gfap, 
                        n_pcs = 20)

# output 
with_pca

plotPcaVariation(object = spata_obj, n_pcs = 20)

plotPca(object = spata_obj, 
        color_by = "seurat_clusters", 
        pt_alpha = 0.75,
        pt_clrp = "npg", 
        pt_size = 0.75,
        n_pcs = 12, 
        nrow = 3
)

plotTsne(object = spata_obj,
         color_by = "seurat_clusters", 
         pt_alpha = 0.75,
         pt_clrp = "npg", 
         pt_size = 0.75,
)

plotUmap(object = spata_obj,
         color_by = "seurat_clusters", 
         pt_alpha = 0.75,
         pt_clrp = "npg", 
         pt_size = 0.75,
)

monocle_clusters <- findMonocleClusters(object = spata_obj, 
                                        preprocess_method = "PCA", 
                                        reduction_method = c("UMAP", "PCA", "tSNE"), 
                                        cluster_method = c("leiden", "louvain"), 
                                        k = 5, 
                                        num_iter = 5)

# output
monocle_clusters
# output
examineClusterResults(data = monocle_clusters)
# feature names before adding
getFeatureNames(spata_obj)
# add the cluster results
spata_obj <- 
  addFeatures(object = spata_obj, 
              feature_names = c("cluster_leiden_UMAP_k5", "cluster_leiden_tSNE_k5","cluster_louvain_PCA_k5"), 
              feature_df = monocle_clusters,
              overwrite = TRUE,
              key = "barcodes")

# feature names afterwards
getFeatureNames(spata_obj)
plotSurface(object = spata_obj,
            color_by = "cluster_leiden_UMAP_k5",
            pt_size = 1.9,
            pt_clrp = "jama") +
  ggplot2::labs(color = "Leiden UMAP")

# rename the clustering variable 
spata_obj <- renameFeatures(object = spata_obj, "Leiden_UMAP" = "cluster_leiden_UMAP_k5")

# print features names
getFeatureNames(spata_obj)


spata_obj <- 
  runDeAnalysis(object = spata_obj,
                across = "Leiden_UMAP", # across which identity groups
                method_de = c("wilcox", "bimod") # with which methods
  )

# get an overview about the de-analysis results stored in your object
printDeaOverview(spata_obj)



cluster_of_interest <- c("2", "3", "5")

genes_of_interest1 <- getDeaGenes(object = spata_obj, 
                                  across = "seurat_clusters", 
                                  across_subset = cluster_of_interest, 
                                  n_highest_lfc = 1 # return the top 1 gene for every group
)

genes_of_interest3 <- getDeaGenes(object = spata_obj, 
                                  across = "seurat_clusters", 
                                  across_subset = cluster_of_interest, 
                                  n_highest_lfc = 3 # return the top 3 genes for every group
)

# output
genes_of_interest1

plotHistogram(object = spata_obj, variables = genes_of_interest3, clrp = "default")
