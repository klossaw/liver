################################## package #####################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(pheatmap)
library(glue)
library(SingleR)
library(jhtools)
library(scCATCH)
################################### QC #########################################
#'@export run_seurat_qc
#'@param filter_RBgene 是否去除核糖体基因
#'@param filter_MTgene 是否去除线粒体基因
run_seurat_qc <- function(object,
                          output_dir,
                          filter_RBgene = FALSE,
                          filter_MTgene = FALSE,
                          percent_mt = 20,
                          max_nFeature_RNA = 2500,
                          min_nFeature_RNA = 200
                       ){
  output_dir %>% checkdir()
  #QC
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object[["percent.rb"]] <- PercentageFeatureSet(object, pattern = "^RP[SL]")
  rb_genes <- rownames(object)[grep("^RP[SL]",rownames(object))]
  mt_genes <- rownames(object)[grep("^MT-",rownames(object))]
  VlnPlot(object, 
          features = c("percent.mt","percent.rb"),
          group.by = "orig.ident", 
          ncol = 3,
          pt.size = 0) %>% 
    ggsave(filename = glue("{output_dir}/percent.png"),width = 10,height = 10)
  #cell-cycle
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  object <- RunPCA(object, features = c(s.genes, g2m.genes))
  DimPlot(object,group.by = "Phase") %>% 
    ggsave(filename = glue("{output_dir}/Phase.png"),width = 10,height = 10)
  Idents(object) <- "seurat_clusters"
  p1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(p1|p2,filename = glue("{output_dir}/Scatter.png"),width = 10,height = 10)
  #select 
  object <- subset(object, 
                   subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < percent_mt)
  if(filter_RBgene == TRUE){
    object <- object[!rownames(object)%in%rb_genes,]
  }
  if(filter_MTgene == TRUE){
    object <- object[!rownames(object)%in%mt_genes,]
  }
  write_rds(object,file = glue("{output_dir}/sce.rds"))
  return(object)
}
################################ singleR  ######################################
#'@export run_singleR
run_singleR <- function(object,
                        output_dir,
                        type = "seurat_clusters",
                        ...){
  output_dir %>% checkdir()
  singleRdata <- object
  Idents(object) <- type
  ref <- read_rds("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/ref/ref.HumanPrimaryCellAtlasData.Rds")
  pred.exp <- GetAssayData(object,slot = "data")
  pred <- SingleR(pred.exp, ref = ref, clusters = object@active.ident,
                        labels = ref$label.main)
  singleR_result <- data.frame(cluster = rownames(pred),
                               cell.pruned.labels = pred$pruned.labels,
                               cell.labels = pred$labels,
                               cell.first.labels = pred$first.labels)
  df2excel(singleR_result,fn.excel = glue("{output_dir}/singleR.xlsx"))
  
  singleRdata[["singleR.pruned.labels"]] <- "NA"
  singleRdata[["singleR.labels"]] <- "NA"
  singleRdata[["singleR.first.labels"]] <- "NA"
  
  for (i in 1:nrow(singleR_result)) {
    singleRdata$singleR.pruned.labels[which(Idents(singleRdata) == singleR_result[i,"cluster"])] <- 
      as.character(singleR_result[i,"cell.pruned.labels"])
    
    singleRdata$singleR.labels[which(Idents(singleRdata) == singleR_result[i,"cluster"])] <- 
      as.character(singleR_result[i,"cell.labels"])
    
    singleRdata$singleR.first.labels[which(Idents(singleRdata) == singleR_result[i,"cluster"])] <- 
      as.character(singleR_result[i,"cell.first.labels"])
    
  }
  DimPlot(singleRdata,reduction = "umap",group.by = "singleR.pruned.labels",label = TRUE)%>% 
    ggsave(filename = glue("{output_dir}/singleR_umap.png"),
           width = 10,height = 10)
  DimPlot(singleRdata,reduction = "tsne",group.by = "singleR.pruned.labels",label = TRUE)%>% 
    ggsave(filename = glue("{output_dir}/singleR_tsne.png"),
           width = 10,height = 10)
  write_rds(singleRdata,file = glue("{output_dir}/singleRdata.rds"))
  return(singleRdata)
}
################################## scCATCH ####################################
#'@export run_scCATCH
run_scCATCH <- function(object,
                        species = "Human",
                        tissue = "Liver",
                        cluster = "seurat_clusters",
                        output_dir){
  output_dir %>% checkdir()
  Idents(object) <- cluster
  assay_default <- DefaultAssay(object)
  message(glue("DefaultAssay:{DefaultAssay(object)}"))
  
  obj_scCATCH <- createscCATCH(data = object@assays[[assay_default]]@data, 
                               cluster = as.character(object@meta.data[,cluster]))
  obj_scCATCH <- findmarkergene(object = obj_scCATCH, 
                                species = species, 
                                marker = cellmatch, 
                                tissue = tissue)
  obj_scCATCH <- findcelltype(object = obj_scCATCH)
  
  new.cluster.ids <- obj_scCATCH@celltype$cell_type
  names(new.cluster.ids) <- obj_scCATCH@celltype$cluster
  object <- RenameIdents(object, new.cluster.ids)
  object@meta.data$scCATCH<-Idents(object)
  pdf(glue("{output_dir}/scCATCH.pdf"),width = 10,height = 10)
  DimPlot(object = object, 
          group.by = "scCATCH",
          label = T)
  dev.off()
  df2excel(dat.excel = obj_scCATCH@celltype,
           fn.excel = glue("{output_dir}/scCATCH.xlsx"),
           showRow = T)
  return(object)
  
}
############################## marker annotation ###############################
run_seurat_marker <- function(object,
                              cluster = "seurat_clusters",
                              output_dir){
  output_dir %>% checkdir()
  #find marker
  if(!file.exists(glue("{output_dir}/allmarker.csv"))){
      Idents(object) <- cluster
      options(future.globals.maxSize = 14000 * 1024^2)
      plan("multiprocess", workers = 64)
      all.marker <- FindAllMarkers(object)
      top5 <- all.marker %>% group_by(cluster) %>% top_n(5,wt = avg_log2FC) 
      top10 <- all.marker %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC) 
      top20 <- all.marker %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
      top30 <- all.marker %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
      top50 <- all.marker %>% group_by(cluster) %>% top_n(50,wt = avg_log2FC)
      save_marker <- list(all.marker,top50,top30,top20,top10)
      names(save_marker) <- c('AllMarker','top50','top30','top20','top10')
      write_csv(all.marker,file = glue("{output_dir}/allmarker.csv"))
      list2excel(save_marker,fn.excel = glue("{output_dir}/marker.xlsx"))
      exp <- AverageExpression(object,assays = "RNA",slot = "data")$RNA %>% 
        as.data.frame() %>% 
        drop_na()%>% 
        filter_all(all_vars(. != Inf))  
      df2excel(dat.excel = exp,fn.excel = glue("{output_dir}/exp.xlsx"))
  }else{
    all.marker <- read_csv(glue("{output_dir}/allmarker.csv"))
    top5 <- all.marker %>% group_by(cluster) %>% top_n(5,wt = avg_log2FC) 
    top10 <- all.marker %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC) 
    top20 <- all.marker %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
    top30 <- all.marker %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
    top50 <- all.marker %>% group_by(cluster) %>% top_n(50,wt = avg_log2FC)
    }
  #featue plot
  top.marker.feature <- split(top10,top10$cluster)
  feature.dir <- file.path(output_dir,"FeaturePlot") %>% checkdir()
  lapply(top.marker.feature, function(x){
    FeaturePlot(object, reduction = "tsne",features = x$gene) %>% 
      ggsave(filename = glue("{feature.dir}/{cluster}{unique(x$cluster)}tsne.png"),width = 20,height = 15)
    FeaturePlot(object, reduction = "umap",features = x$gene) %>% 
      ggsave(filename = glue("{feature.dir}/{cluster}{unique(x$cluster)}umap.png"),width = 20,height = 15)
  })
  
  #vln plot
  VlnPlot.dir <- file.path(output_dir,"VlnPlot") %>% checkdir()
  lapply(top.marker.feature, function(x){
    VlnPlot(object, group.by = cluster,features = x$gene) %>% 
      ggsave(filename = glue("{VlnPlot.dir}/{cluster}{unique(x$cluster)}vln.png"),width = 10,height = 10)
  })
  
  #dot plot
  Dot.dir <- file.path(output_dir,"Dotplot") %>% checkdir()
  lapply(top.marker.feature, function(x){
    DotPlot(object,
            features = x$gene,
            group.by = cluster,
            cols = c("lightgrey", "blue")) %>% 
      ggsave(filename = glue("{Dot.dir}/{cluster}{unique(x$cluster)}dot.png"),width = 10,height = 10)
  })
  #Dimplot
  Dimplot.dir <- file.path(output_dir,"Dimplot") %>% checkdir()
  DimPlot(object,reduction = "umap",group.by = cluster,label = TRUE)%>% 
    ggsave(filename = glue("{Dimplot.dir}/{cluster}_umap.png"),
           width = 10,height = 10)
  DimPlot(object,reduction = "tsne",group.by = cluster,label = TRUE)%>% 
    ggsave(filename = glue("{Dimplot.dir}/{cluster}_tsne.png"),
           width = 10,height = 10)
  DimPlot(object,reduction = "umap",group.by = cluster,split.by = "orig.ident",label = TRUE) %>% 
    ggsave(filename = glue("{Dimplot.dir}/{cluster}_sample_umap.png"),
           width = 4*length(unique(object$orig.ident)),height = 8)
  DimPlot(object,reduction = "tsne",group.by = cluster,split.by = "orig.ident",label = TRUE) %>% 
    ggsave(filename = glue("{Dimplot.dir}/{cluster}_sample_tsne.png"),
           width = 4*length(unique(object$orig.ident)),height = 8)
  if(!is.null(object@meta.data$orig.ident)){
  DimPlot(object,reduction = "umap",group.by = cluster,split.by = "group",label = TRUE) %>% 
    ggsave(filename = glue("{Dimplot.dir}/{cluster}_group_umap.png"),
           width = 4*length(unique(object$group)),height = 8)
  DimPlot(object,reduction = "tsne",group.by = cluster,split.by = "group",label = TRUE) %>% 
    ggsave(filename = glue("{Dimplot.dir}/{cluster}_group_tsne.png"),
           width = 4*length(unique(object$group)),height = 8)
  }
  return(object)
}
############################## marker plot #####################################
run_marker_plot <- function(object,
                            input_marker,
                            output_dir,
                            type = "seurat_clusters"){
  output_dir %>% checkdir()
  gene <- intersect(row.names(object),input_marker)
  if(length(gene)>10){
    gene_matrix <- matrix(gene,ncol = 4) 
    lapply(1:nrow(gene_matrix),function(x){
      gene <- gene_matrix[x,]
      FeaturePlot(object, reduction = "tsne",features = gene) %>% 
        ggsave(filename = glue("{output_dir}/{type}{x}tsne.png"),width = 10,height = 10)
      FeaturePlot(object, reduction = "umap",features = gene) %>% 
        ggsave(filename = glue("{output_dir}/{type}{x}umap.png"),width = 10,height = 10)
      #vln plot
      VlnPlot(object, group.by = type,features = gene) %>% 
        ggsave(filename = glue("{output_dir}/{type}{x}vln.png"),width = 10,height = 10)
      #dot plot
      DotPlot(object,features = gene,group.by = type,cols = c("lightgrey", "blue")) %>% 
        ggsave(filename = glue("{output_dir}/{type}{x}dot.png"),width = 10,height = 10)  
    })
  }
  FeaturePlot(object, reduction = "tsne",features = gene) %>% 
      ggsave(filename = glue("{output_dir}/{type}tsne.png"),width = 10,height = 10)
  FeaturePlot(object, reduction = "umap",features = gene) %>% 
      ggsave(filename = glue("{output_dir}/{type}umap.png"),width = 10,height = 10)
  #vln plot
  VlnPlot(object, group.by = type,features = gene) %>% 
      ggsave(filename = glue("{output_dir}/{type}vln.png"),width = 10,height = 10)
  #dot plot
  DotPlot(object,features = gene,group.by = type,cols = c("lightgrey", "blue")) %>% 
      ggsave(filename = glue("{output_dir}/{type}dot.png"),width = 10,height = 10)
  #heatmap plot
  exp <- AverageExpression(object,assays = "RNA",slot = "data")$RNA %>% 
    as.data.frame() %>% 
    drop_na()%>% 
    filter_all(all_vars(. != Inf))  
  df2excel(dat.excel = exp,fn.excel = glue("{output_dir}/{type}exp.xlsx"))
  
  exp[gene,] %>% 
      as.matrix() %>% 
      pheatmap(scale = "row",
               show_rownames = FALSE,
               filename = glue("{output_dir}/{type}heatmap.png"))
  return(object)
}
############################## add label #######################################
#'@export
#'@import 
sce_add_label <- function(object,
                          label,
                          idents_cluster = 'seurat_clusters'){
  Idents(object) <- idents_cluster
  object@meta.data$cell_types <- 'NA'
  for (i in 1:nrow(label)) {
    object@meta.data$cell_types[which(Idents(object) == as.character(label[i,1]))] <- 
      as.character(label[i,2])
  }
  Idents(object) <- 'cell_types'
  object@meta.data$cell_types <- factor(Idents(object),levels = unique(label$cell_types))
  return(object)
}


label <- read_csv('05AddLabel/label.csv')
sce_addlabel <- sce_add_label(object = sce,label = label)
DimPlot(sce_addlabel,group.by = 'cell_types',label = T)
################################### plot_label #################################
run_plot_label <- function(object,output_dir){
  DimPlot(object,reduction = "umap",group.by = "cell_type",label = TRUE)%>% 
    ggsave(filename = glue("{output_dir}/cell_type_umap.png"),width = 10,height = 10)
  DimPlot(object,reduction = "tsne",group.by = "cell_type",label = TRUE)%>% 
    ggsave(filename = glue("{output_dir}/cell_type_tsne.png"),width = 10,height = 10)
  if(!is.null(object@meta.data$group)){
    DimPlot(object,reduction = "umap",group.by = "cell_type",split.by = "group",label = TRUE) %>% 
      ggsave(filename = glue("{output_dir}/cell_type_group_umap.png"),width = 8*length(unique(object$group)),height = 8)
    DimPlot(object,reduction = "tsne",group.by = "cell_type",split.by = "group",label = TRUE) %>% 
      ggsave(filename = glue("{output_dir}/cell_type_group_tsne.png"),width = 8*length(unique(object$group)),height = 8)
  }
  DimPlot(object,reduction = "umap",group.by = "cell_type",split.by = "orig.ident",label = TRUE) %>% 
    ggsave(filename = glue("{output_dir}/cell_type_sample_umap.png"),width = 6*length(unique(object$orig.ident)),height = 12)
  DimPlot(object,reduction = "tsne",group.by = "cell_type",split.by = "orig.ident",label = TRUE) %>% 
    ggsave(filename = glue("{output_dir}/cell_type_sample_tsne.png"),width = 6*length(unique(object$orig.ident)),height = 12)
  #plot proportion
  colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", 
              "#8491B4FF", "#91D1C2FF",  "#7E6148FF", "#B09C85FF", 
              colorRampPalette(brewer.pal(8,'Dark2'), alpha = 1)(8))
  meta.data <- object@meta.data
  cluster.color <- colors[1:length(unique(meta.data$cell_type))]
  
  p1 <- ggplot(data = meta.data,mapping = aes(x = orig.ident, fill=cell_type))+
    geom_bar(stat = "count",width = 0.7, position = 'fill') +
    scale_fill_manual(values = cluster.color) +
    labs(y = "proportion of cells", x = "") +
    theme(text = element_text(family = "serif"),
          axis.text = element_text(family = "serif", size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
          panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line.x = element_line (colour = "black"), axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = NA, color = NA), legend.position = "left") 
  ggsave(p1,filename = glue("{output_dir}/sample_cell_type_stat_bar.png"), width = 8, height = 6)
  return(object) 
}
############################### run sub celltype ###############################
run_sub_type <- function(object,
                         sub_cluster = 'T cell',
                         type = 'cell_type',
                         dims = 1:15,
                         resolution = c(0.5,1.5,0.8)){
    Idents(object) <- type
    sub_object <- subset(object,idents = sub_cluster)
    sub_object <- SCTransform(sub_object) %>% 
      RunPCA() %>% 
      RunUMAP(dims = dims) %>% 
      RunTSNE(dims = dims) %>% 
      FindNeighbors(reduction = "pca", dims = dims) %>% 
      FindClusters(resolution = resolution)
    return(sub_object)
}
################################ celltypeist ###################################
run_celltypist <- function(object){
  adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(object[['RNA']]@counts)))),
                         obs = pandas$DataFrame(object@meta.data),
                         var = pandas$DataFrame(data.frame(gene = rownames(object[['RNA']]@counts),
                                                           row.names = rownames(object[['RNA']]@counts)))
  )
  scanpy$pp$normalize_total(adata, target_sum=1e4)
  scanpy$pp$log1p(adata)
  for(x in model_name){
    predictions = celltypist$annotate(adata, model = model_list[[x]], majority_voting = T)
    object  = AddMetaData(object, predictions$predicted_labels$majority_voting, col.name =x) 
  }
  return(object)
}
################################ run ###########################################
if(sys.nframe()==0){
  work_dir <- c("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/")
  project_dir <- file.path(work_dir,"SingleCellAnnotation") %>% checkdir()
  setwd(project_dir)
  #'@example
  sce <- readRDS("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/seurat/merge/harmony/merged_harmony.rds")
  # add metadata
  Normal <- c("hcc1C","hcc3C","hcc4C")
  Tumor <- c("hcc1T","hcc2T","hcc3T","hcc4T","hcc5T")
  sce@meta.data$group[sce$orig.ident %in% Normal] <- "Normal"
  sce@meta.data$group[sce$orig.ident %in% Tumor] <- "Tumor"
  #filter
  qc_dir <- file.path(project_dir,"01QC") %>% checkdir()
  sce <- run_seurat_qc(object = sce, 
                       output_dir = file.path(qc_dir,"filter0.2"),
                       filter_MTgene = FALSE,
                       filter_RBgene = FALSE,
                       percent_mt = 20)
  #'@example scCATCH
  #'https://github.com/ZJUFanLab/scCATCH/wiki/human_tissues
  scCATCH_dir = file.path(project_dir,"02scCATCH")
  run_scCATCH(object = sce,cluster = "seurat_clusters",output_dir = scCATCH_dir)
  
  #'@example
  singleR_dir = file.path(project_dir,"03singleR")
  run_singleR(object = sce,output_dir = singleR_dir)
  
  #'@example
  markerannotation_dir <- file.path(project_dir,"04markerannotation")
  output_dir <- file.path(markerannotation_dir,"seurat_clusters")
  run_seurat_marker(object = sce,
                    cluster = "seurat_clusters",
                    output_dir = output_dir)
  #'@example 
  immuneMarker <- read.csv('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/markerannotation/immune_marker.csv',header = FALSE)
  immune.Marker <- strsplit(immuneMarker$V2,split = ",")
  immune.Marker <- lapply(immune.Marker, function(x){
    str_trim(x)
  })
  names(immune.Marker) <- immuneMarker$V1
  markerannotation_dir <- file.path(project_dir,"markerannotation")
  
  lapply(names(immune.Marker), function(x){
    run_marker_plot(object = filter_0.2,input_marker = immune.Marker[[x]],output_dir = file.path(markerannotation_dir,x))
  })
  
  #'@example add label
  output_dir <- file.path(project_dir,"05AddLabel")
  label <- read_csv('05AddLabel/label.csv')
  sce_addlabel <- sce_add_label(object = sce,label = label)
  DimPlot(sce_addlabel,group.by = 'cell_type',label = T)
  #'@example select subcluster annotation
  library(reticulate)
  reticulate::use_condaenv('/cluster/home/yzy_jh/.conda/envs/scell/')
  reticulate::py_module_available('celltypeist')
  scanpy = import("scanpy")
  celltypist = import("celltypist")
  pandas <- import("pandas")
  numpy = import("numpy")
  
  #### 下载参考数据集
  setwd('ref/')
  celltypist$models$download_models(force_update = T) 
  celltypist$models$celltypist_path
  
  # 加载所有的参考数据集
  model_dir = list.files("/cluster/home/yzy_jh/.celltypist/data/models",full.names = T)
  model_name = list.files("/cluster/home/yzy_jh/.celltypist/data/models",full.names = F) %>% 
    str_replace_all(c('.pkl'=''))
  model_list = lapply(model_dir, function(x){
    celltypist$models$Model$load(model = x)})
  names(model_list ) = model_name
  ####. 2.seurat to celltypist
  sub_T <- subset(sce_annotation,ident=c('Tcell','NKcell','NKcell/Tcell'))
  seurat.data <- run_celltypist(sub_T)
  
  p0 = DimPlot(seurat.data,group.by = "cell_type", reduction = "umap", label = TRUE)
  p1 = DimPlot(seurat.data,group.by = "Immune_All_High", reduction = "umap", label = TRUE) 
  p2 = DimPlot(seurat.data,group.by = "Immune_All_Low", reduction = "umap", label = TRUE) 
  p0 + p1 + p2 
  
  ################################################################################
  plan("multiprocess", workers = 64)
  Idents(sce_annotation) <- "cell_type"
  marker <- FindAllMarkers(sce_annotation)
  top20 <- marker %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
  
  Bcell <- c('CD19', 'CD79A', 'MS4A1', 'CD22')
  Tcell <- c('CD3D', 'CD3E', 'CD3G')          
  NKcell <-c('NKG7', 'FCGR3A', 'IFNG', 'GNLY', 'GZMB')           
  HSC_MFB <-  c('RGS5', 'COL1A1','ACTA2','PDGFRB')   
  Endothelialcell <-   c('PECAM1','CDH5','STC1','TM4SF1')          
  Monocyte <- c('CST3', 'LYZ','CD1C', 'FCER1A', 'FCER1G')
  Macro_Kupffercells  <- c('CST3', 'LYZ', 'CD68', 'CD163', 'FCER1G')
  Hepatocytes <- c('APOA2','ALB','APOC1','APOC3')                                
  DC <- c('CST3','LYZ', 'CD1C','IDO1')
  library(RColorBrewer)
  display.brewer.pal(n = 8, name = 'Dark2')
  display.brewer.all()
  col=brewer.pal(n = 3, name = "YlGn")
  
  FeaturePlot(sce_annotation,features = Tcell) %>% ggsave(filename = "Tfeature.png",width = 10,height = 10)
  FeaturePlot(sce_annotation,features = Bcell) %>% ggsave(filename = "Bfeature.png",width = 10,height = 10)
  FeaturePlot(sce_annotation,features = NKcell)%>% ggsave(filename = "NKcell.png",width = 10,height = 10)
  FeaturePlot(sce_annotation,features = HSC_MFB)%>% ggsave(filename = "HSC_MFB.png",width = 10,height = 10)
  FeaturePlot(sce_annotation,features = Endothelialcell) %>% ggsave(filename = "Endothelialcell.png",width = 10,height = 10)
  FeaturePlot(sce_annotation,features = Monocyte) %>% ggsave(filename = "Monocyte.png",width = 10,height = 10)
  FeaturePlot(sce_annotation,features = Macro_Kupffercells) %>% ggsave(filename = "Macro_Kupffercells.png",width = 10,height = 10)
  FeaturePlot(sce_annotation,features = Hepatocytes) %>% ggsave(filename = "Hepatocytes.png",width = 10,height = 10)
  FeaturePlot(sce_annotation,features = DC)%>% ggsave(filename = "HSC_MFB.png",width = 10,height = 10)
  
  genes <- c(Tcell,NKcell,HSC_MFB,Endothelialcell,Monocyte,Hepatocytes,Bcell,DC)
  DotPlot(sce_annotation,features = unique(genes),cols = c('red','blue')) %>% ggsave(filename = "dot.png",width = 25,height = 10)
}
