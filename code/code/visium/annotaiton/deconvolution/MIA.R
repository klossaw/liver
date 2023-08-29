library(magrittr)
setwd('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/deconvolution/MIA/')
st_data <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/seurat/merged.rds')
sc_data <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/05AddLabel/sce_addlabel.rds')
load('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/seurat/Sample/patient.Rdata')

sp_marker <- parallel::mclapply(paitent_list,function(x){
  
  marker <- NormalizeData(x) %>% 
            FindAllMarkers(assay = 'Spatial', only.pos = T) %>% as_tibble
  return(marker)
  
  },mc.cores = 16)

### find all markers of sc and st.
Idents(sc_data) <- "cell_subtypes"
Idents(st_data) <- "SCT_snn_res.0.5"
DefaultAssay(st_data) <- 'Spatial'
Idents(sce) <- 'cell_types'
  
sc.markers <- FindAllMarkers(object = sc_data, only.pos = T) %>% as_tibble
st.markers <- st_data %>% NormalizeData() %>% FindAllMarkers(assay = 'Spatial', only.pos = T) %>% as_tibble
### 按照FDR过滤
MIA <- function(st.markers,sc.markers){
  sc.markers %<>% filter(p_val < 0.05)
  st.markers %<>% filter(p_val < 0.05)
  
  sc.gene <- GetAssayData(sc_data) %>% rownames
  st.gene <- GetAssayData(st_data) %>% rownames
  sc.num <- length(sc.gene)
  st.num <- length(st.gene)
  
  background.gene <- intersect(sc.gene,st.gene)
  background.num <- length(background.gene)
  
  ### 使用背景基因过滤markers
  sc.markers %<>% filter(gene %in% background.gene)
  st.markers %<>% filter(gene %in% background.gene)
  
  scst.markers <- inner_join(sc.markers,st.markers,by="gene",suffix=c("_SC","_ST"))
  scst.markers %<>% mutate(cluster_SCST=paste0(cluster_SC,",",cluster_ST))
  sc.count.cluster <- sc.markers %>% dplyr::count(cluster)
  st.count.cluster <- st.markers %>% dplyr::count(cluster)
  scst.count.cluster <- scst.markers %>% dplyr::count(cluster_SCST) %>% separate(cluster_SCST, sep=",", into=c("cell_types","cluster"))
  scst.cross <- scst.count.cluster %>% spread(cluster,n,fill=0)
  
  scst.count.pvalue <- scst.count.cluster %>% inner_join(sc.count.cluster,by=c(cell_types="cluster"),suffix=c("","_sc"))  %>% inner_join(st.count.cluster,by=c(cluster="cluster"),suffix=c("","_st")) %>% mutate(n.bk = background.num, p_val = phyper(n-1, n_sc, background.num-n_sc, n_st, lower.tail = FALSE))
  scst.count.pvalue$p_val %<>% sapply(function(x)max(x,1e-200))
  
  scst.count.pvalue.table <- scst.count.pvalue %>% dplyr::select(cell_types,cluster,p_val) %>% data.frame %>% spread(cluster,p_val,fill=1)
  rownames(scst.count.pvalue.table) <- scst.count.pvalue.table$cell_types
  scst.count.pvalue.table <- scst.count.pvalue.table %>% dplyr::select(-cell_types)
  
  scst.count.pvalue.matrix <- scst.count.pvalue.table %>% as.matrix
  scst.count.score.matrix <- scst.count.pvalue.matrix %>% -log10(.)
  header <- c("cell_types", colnames(scst.count.pvalue.matrix))
  return(list(scst.count.score.matrix,
              scst.count.pvalue.matrix))
}

patient <- c('shcc1','shcc2','shcc3','shcc4','shcc5')
paitent_list <- list(SP_shcc1,SP_shcc2,SP_shcc3,SP_shcc4,SP_shcc5)

names(sp_marker) <- patient
paitent_list <- lapply(paitent_list, function(x){
  return(FindClusters(x,resolution = 0.3))
})
names(paitent_list) <- patient

SpatialDimPlot(paitent_list$shcc1)
SpatialDimPlot(paitent_list$shcc2)
SpatialDimPlot(paitent_list$shcc3)
SpatialDimPlot(paitent_list$shcc4)
SpatialDimPlot(paitent_list$shcc5,alpha = c(0.1,1.6))

mia_data <- lapply(patient,function(x){
  mia_data <- MIA(sc.markers = sc.markers,st.markers = sp_marker[[x]])
  p.out.score <- pheatmap::pheatmap(mia_data[[1]],scale = 'row',
                                    cluster_rows = FALSE,
                                    cluster_cols = FALSE,
                                    display_numbers = TRUE)
  ggsave(p.out.score,
         filename = glue("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/deconvolution/MIA/fig/{x}.png"),
         width = 10,height = 10)
  })

RColorBrewer::brewer.pal.info
cols <- c('0'="#E41A1C",
          '1'="#377EB8",
          '2'="#4DAF4A",
          '3'="#984EA3",
          '4'="#FF7F00",
          '5'="#FFFF33",
          '6'="#A65628",
          '7'="#F781BF",
          '8'="#999999")

samples <- c("shcc1C","shcc1T", "shcc2T", "shcc3C", "shcc3T", "shcc4C", "shcc4T", "shcc5T")
names(st_data@images) <- samples
lapply(samples,function(x){
  SpatialDimPlot(st_data,images = x,label = T,cols  =  cols) %>% 
    ggsave(filename = glue("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/deconvolution/MIA/cluster/{x}_cluster.png"),width=10,height=10)
  
})

DimPlot(st_data)
SpatialDimPlot(st_data)


st.markers  <- marker_sp5
sc.markers <-  

for(i in unique(sub_Myeloid_label$cluster)){
  if(sub_Myeloid_marker$cluster==i){
    sub_Myeloid_marker$cell=sub_Myeloid_label$cluster[i,]
  }
}

sub_Myeloid_marker$cell <- 'NA'
for (i in 1:nrow(sub_Myeloid_label)) {
  sub_Myeloid_marker[sub_Myeloid_marker$cluster == as.character(sub_Myeloid_label[i,1]),]$cluster <- 
    as.character(sub_Myeloid_label[i,2])
}
sub_Myeloid_marker$cluster <- sub_Myeloid_marker$cell

sub_Myeloid_marker$cell <- 'NA'
for (i in 1:nrow(sub_Myeloid_label)) {
  sub_Myeloid_marker[sub_Myeloid_marker$cluster == as.character(sub_Myeloid_label[i,1]),]$cluster <- 
    as.character(sub_Myeloid_label[i,2])
}
sub_Myeloid_marker$cluster <- sub_Myeloid_marker$cell

NK_T_marker$cell <- 'NA'
for (i in 1:nrow(NKT_label)) {
  NK_T_marker[NK_T_marker$cluster == as.character(NKT_label[i,1]),]$cluster <- 
    as.character(NKT_label[i,2])
}
NK_T_marker$cluster <- NK_T_marker$cell

names(st_data@images) <- samples

Idents(sce) <- 'cell_subtypes'
all_markers <- FindAllMarkers(sce)

mia_data_all <- MIA(sc.markers = all_marker,st.markers = st.markers)
mia_data_Myeloid <- MIA(sc.markers = sub_Myeloid_marker,st.markers = st.markers)
mia_data_NK_T <- MIA(sc.markers = NK_T_marker,st.markers = st.markers)
mia_data_Hepatocytes <- MIA(sc.markers = sub_Hepatocytes_marker,st.markers = st.markers)

p_alltype <- pheatmap::pheatmap(mia_data_all[[1]],
                                  cluster_rows = FALSE,
                                  cluster_cols = FALSE,
                                  display_numbers = TRUE,scale = 'row')

p_Myeloid <- pheatmap::pheatmap(mia_data_Myeloid[[1]],
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                display_numbers = TRUE,scale = 'row')
p_NK_T <- pheatmap::pheatmap(mia_data_NK_T[[1]],
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                display_numbers = TRUE,scale = 'row')
p_Hepatocytes <- pheatmap::pheatmap(mia_data_Hepatocytes[[1]],
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                display_numbers = TRUE,scale = 'row')


table(st_data$orig.ident,st_data$seurat_clusters)

save(sp_marker,sub_T,sub_Hepatocytes,file = '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/deconvolution/MIA/MIA.Rdata')
p.out.score <- pheatmap::pheatmap(scst.count.score.matrix,
                                  color = ,
                                  cluster_rows = FALSE,
                                  cluster_cols = FALSE,
                                  display_numbers = TRUE)

p.out.pvalue <- pheatmap::pheatmap(scst.count.pvalue.matrix, 
                                   colorRampPalette(c("GhostWhite", "#B0C4DE", "#4682B4"))(100),
                                   cluster_rows = FALSE,
                                   cluster_cols = FALSE,
                                   display_numbers = TRUE,
                                   cellwidth = 45, 
                                   cellheight = 45)
p.out.pvalue 

SpatialDimPlot(st_data)
FeaturePlot(st_data,features = )

table(st_data$seurat_clusters,st_data$orig.ident)
### 在空间中展示单细胞中各种细胞类型的top3 marker基因表达量。
dir.create(file.path(cfg$od,"BMK2_MIA","top3Marker"),recursive=T)
sc.markers %>% split.data.frame(sc.markers$cluster) %>% 
  lapply(function(tb){
    top.gene <- slice_max(tb,avg_log2FC,n=3)$gene
    p.out <- SpatialFeaturePlot(st_data, features = top.gene)
    celltype <- tb$cluster[1] %>% str_replace_all(" ","_") 
    ggsave(filename = file.path(cfg$od,"BMK2_MIA","top3Marker",paste(celltype, "top3Marker.png", sep = "_")), plot = p.out, width = 6, height = 4, scale = 1.3)
    ggsave(filename = file.path(cfg$od,"BMK2_MIA","top3Marker",paste(celltype, "top3Marker.pdf", sep = "_")), plot = p.out, width = 6, height = 4, scale = 1.3)
  })


## AddModulescore
p1 <- FeaturePlot(st_data,features = names(features))
p2 <- DimPlot(st_data,label = T)
p <- p2|p1
p  
ggsave(p,filename = glue("{AddModulescore_dir}/spatialann.png"),width = 10,height = 8)

Seurat::RidgePlot(st_data,features = names(features))


Idents(sce) <- 'cell_types'
all_marker <- FindAllMarkers(sce)

features <- all_marker %>% group_by(cluster) %>% top_n(20,avg_log2FC) %>% split.data.frame(.$cluster) %>% lapply(function(tb){tb$gene})
st_data <- AddModuleScore(st_data, features = features, name = "ModuleScore")


a <- length(st_data@meta.data)-length(features)+1
b <- length(st_data@meta.data)
colnames(st_data@meta.data)[a:b] <- names(features)
p.out <- st_data %>% VlnPlot(names(features),pt.size=0) & NoLegend() + theme(axis.title = element_blank(), axis.text = element_text(size = 8), title = element_text(size = 10))
p.out

AddModulescore_dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/deconvolution/AddModulescore/shcc5T"
checkdir(AddModulescore_dir)
for (i in names(features)) {
  p1 <- SpatialFeaturePlot(st_data,features = i,alpha = c(0,1.5)) 
  ggsave(p1,filename = glue("{AddModulescore_dir}/{i}.png"),width = 12,height = 6)
 
  p2 <- FeaturePlot(st_data,features = i)
  ggsave(p2,filename = glue("{AddModulescore_dir}/{i}_feature.png"),width = 6,height = 6)
}
for (i in names(features)) {
p0 <- SpatialFeaturePlot(st_data,features = i,alpha = c(0,1.5),images = 'shcc5T') 
ggsave(p0,filename = glue("{AddModulescore_dir}/{i}shcc5T.png"),width = 6,height = 6)
}
ggsave(filename = file.path(cfg$od,"BMK2_MIA", paste("marker_score.png", sep = "_")), plot = p.out, width = 6, height = 4, scale = 1.3)
ggsave(filename = file.path(cfg$od,"BMK2_MIA", paste("marker_score.pdf", sep = "_")), plot = p.out, width = 6, height = 4, scale = 1.3)



save(sc.markers,st.markers,sp_marker,st_data,sc_data,file = 'marker.Rdata')
