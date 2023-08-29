##load
suppressMessages(library(ggplot2))
suppressMessages(library(monocle))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(glue))
#suppressMessages(detach("package:monocle3",unload = TRUE))
#'@import monocle
#'@export monocle2_order
################################ monocle2 order ################################
monocle2_order <- function(object,
                           cellname = "all",
                           type = "Variablegene",
                           markergene = NULL,
                           definegene = NULL){
  DefaultAssay(object)<-'RNA'
  data <- as(as.matrix(object@assays$RNA@data), "sparseMatrix")
  pd <- new("AnnotatedDataFrame", data = object@meta.data)
  fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit =1,
                                expressionFamily = negbinomial.size())
  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)
  monocle_cds <- detectGenes(monocle_cds, min_expr = 1)
  expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= 10))
  monocle_cds <-  monocle_cds[expressed_genes, ]
  disp_table <- dispersionTable(monocle_cds)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 1)
  monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
  if(type == "Variablegene"){
    Variablegenes <- VariableFeatures(object,assay = "RNA")
    if(length(Variablegenes) == 0){
      selectgene <- VariableFeatures(object,assay = "SCT")
    }else{
      selectgene <- VariableFeatures(object,assay = "RNA")
    }
  }else if(type == "markergene"){
    selectgene <- markergene
  }else if(type == "definegene"){
    selectgene <- definegene
  }
  selectgene <- intersect(selectgene, rownames(monocle_cds))
  print(length(selectgene))
  monocle_cds <- reduceDimension(monocle_cds[selectgene, ], max_components = 2, 
                                 num_dim = 6, reduction_method = "tSNE")
  monocle_cds <- reduceDimension(monocle_cds[selectgene, ], max_components = 2,
                                 reduction_method = "DDRTree")
  monocle_cds <- orderCells(monocle_cds)
  return(monocle_cds)
}
#################################### monocle2 plot #############################
#'@import 
#'@export monocle2_plot
monocle2_plot <- function(monocle_cds,
                          output_dir,
                          genes = NULL,
                          width = 10,
                          height = 10,
                          ...){
  
  #gene
  if(!is.null(genes)){
    genes <- row.names(subset(fData(monocle_cds),
                                     gene_short_name %in% genes))
    gene_dir <- glue("{output_dir}/plot_gene") %>% checkdir()
    cds_subset <- monocle_cds[genes[1:10],]
    pdf(glue("{gene_dir}/gene_in_pseudotime.pdf"),width = 10,height = 10)
    p <- plot_genes_in_pseudotime(cds_subset, color_by = "cell_types")
    print(p)
    dev.off()
    pdf(glue("{gene_dir}/gene_branched_pseudotime.pdf"),width = 10,height = 10)
    p <- plot_genes_branched_pseudotime(cds_subset,
                                        branch_point = 1,
                                        color_by = "cell_types",
                                        ncol = 1)
    print(p)
    dev.off()
  }
 
  #cell
  cell_dir <- glue("{output_dir}/plot_cell") %>% checkdir()
  
  plot_cell_trajectory(monocle_cds, color_by = "State") %>% ggsave(filename = glue("{cell_dir}/State.png"),width = 8,height = 6)
  plot_cell_trajectory(monocle_cds, color_by = "Pseudotime") %>% ggsave(filename = glue("{cell_dir}/Pseudotime.png"),width = 8,height = 6)
  plot_cell_trajectory(monocle_cds, color_by = "cell_types") %>% ggsave(filename = glue("{cell_dir}/cell_types.png"),width = 8,height = 6)
  plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters") %>% ggsave(filename = glue("{cell_dir}/seurat_clusters.png"),width = 8,height = 6)
  
  #heatmap
  heatmap_dir <- glue("{output_dir}/plot_heatmamp") %>% checkdir()
  
  select_genes <- row.names(subset(fData(monocle_cds),
                                   gene_short_name %in% genes))
  diff_test_res <- differentialGeneTest(monocle_cds[select_genes,],
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
  
  pdf(glue("{heatmap_dir}/pseudotime_heatmap.pdf"),width = 10,height = 15)
  plot_pseudotime_heatmap(monocle_cds[sig_gene_names,],
                          num_clusters = 5,
                          cores = 1,
                          show_rownames = T)
  dev.off()
  
  BEAM_res <- BEAM(monocle_cds, branch_point = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  pdf(glue("{heatmap_dir}/branched_heatmap.pdf"),width = 10,height = 15)
  
  branch_gene <- row.names(subset(BEAM_res,qval < 1e-4))
  monocle::plot_genes_branched_heatmap(monocle_cds[branch_gene,],
                              branch_point = 1,
                              num_clusters = 4,
                              cores = 1,
                              use_gene_short_name = T,
                              show_rownames = T
                               )
  
  dev.off()
}
################################### main #######################################

  #sce_list<- SplitObject(sce,split.by = "orig.ident")
  #monocle_cds <- monocle2_order(object = ,project_dir = project_dir ,cellname = "HSC",
  #                              type = "markergene",
  #                              markergene=unique(top50$gene))
  #monocle2_plot(monocle_cds = monocle_cds, proejct_dir = project_dir)
  

