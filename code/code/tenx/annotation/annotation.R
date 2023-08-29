packages <- c('Seurat','ggplot2','patchwork','dplyr','RColorBrewer','pheatmap','glue',
              'SingleR','jhtools','scCATCH','reticulate')
suppressMessages(lapply(packages,library,character.only = T))

sce_scCATCH <- function(object,species = "Human",tissue = "Liver",
                        cluster = "seurat_clusters",output_dir){
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

sce_add_label <- function(object,label,idents_cluster = 'seurat_clusters'){
  Idents(object) <- idents_cluster
  object@meta.data$cell_types <- 'Unknow'
  for (i in 1:nrow(label)) {
    object@meta.data$cell_types[which(Idents(object) == as.character(label[i,1]))] <- 
      as.character(label[i,2])
  }
  Idents(object) <- 'cell_types'
  object@meta.data$cell_types <- factor(Idents(object),levels = unique(label$cell_types))
  if('delete'%in%label[,2]){
    object <- object[,!object@meta.data$cell_types%in%'delete']
  }
  return(object)
}

sce_celltypist <- function(object,
                           model_dir = '/cluster/home/yzy_jh/.celltypist/data/models/Immune_All_Low.pkl'){
  reticulate::use_condaenv('/cluster/home/yzy_jh/.conda/envs/scell/')
  scanpy = reticulate::import("scanpy")
  celltypist = reticulate::import("celltypist")
  pandas <- reticulate::import("pandas")
  numpy = reticulate::import("numpy")
  adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(object[['RNA']]@counts)))),
                         obs = pandas$DataFrame(object@meta.data),
                         var = pandas$DataFrame(data.frame(gene = rownames(object[['RNA']]@counts),
                                                           row.names = rownames(object[['RNA']]@counts)))
  )
  scanpy$pp$normalize_total(adata, target_sum=1e4)
  scanpy$pp$log1p(adata)
  model <- celltypist$models$Model$load(model = model_dir)
  predictions = celltypist$annotate(adata, model = model, majority_voting = TRUE)
  object  = AddMetaData(object, predictions$predicted_labels$majority_voting, col.name ='prediction') 
  return(object)
}
#test
scCATCH_dir = file.path(project_dir,"02scCATCH")

sce <- sce_scCATCH(object = sce,cluster = "seurat_clusters",output_dir = scCATCH_dir)
DimPlot(sce,group.by = 'seurat_clusters',label = T)

sce_celltypist <- sce_celltypist(sce)

table(sce_celltypist$prediction,sce$seurat_clusters)

label <- read_csv('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/05AddLabel/label.csv')
sce <- sce_add_label(object = sce,label = label)

sub_T <- subset(sce,ident=c('Tcell','NKcell','NKcell/Tcell'))
seurat.data <- sce_celltypist(sub_T)
DimPlot(seurat.data,group.by = '')
