work_dir <- c("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/monocle/")
project_dir <- file.path(work_dir,"monocle3")

library(monocle3)
################################# monocle3_order ###############################
monocle3_order <- function(object,
                           root_cell = NULL,
                           ...){
  
  data <- Seurat::GetAssayData(object, assay = 'RNA', slot = 'counts')
  cell_metadata <- object@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- monocle3::new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
  cds <- monocle3::preprocess_cds(cds, num_dim = 30)
  #pca
  cds <- monocle3::reduce_dimension(cds, preprocess_method = "PCA")
  # load seurat UMAP
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Seurat::Embeddings(object, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
  #cluster cell
  cds <- monocle3::cluster_cells(cds = cds, reduction_method = "UMAP")
  #learn the trajectory graph
  cds <- monocle3::learn_graph(cds)
  if(!is.null(root_cell)){
    root_cells = row.names(cds@colData[cds@colData$cell_type == root_cell,])
    cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells) 
  }else{
    cds <- monocle3::order_cells(cds, reduction_method = "UMAP")
  }
  return(cds)
}
################################# monocle3 plot ################################
monocle3_plot(){
  
}
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p
##细胞按拟时排序
# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
p + geom_vline(xintercept = seq(-7,-6,0.25)) + geom_hline(yintercept = seq(0,1,0.25))
embed <- data.frame(Embeddings(scRNAsub, reduction = "umap"))
embed <- subset(embed, UMAP_1 > -6.75 & UMAP_1 < -6.5 & UMAP_2 > 0.24 & UMAP_2 < 0.25)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)
   
##寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#基因表达趋势图
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="predicted.id", 
                         min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
##寻找共表达模块
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$predicted.id)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")

save(cds,file = "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/monocle/monocle3/monocle3_cds.rds")

