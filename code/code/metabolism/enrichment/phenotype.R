
HallMarkers <- GSEABase::getGmt('/cluster/home/yzy_jh/ref/MSigDB/h.all.v7.5.1.symbols.gmt')
KEGG <- GSEABase::getGmt('/cluster/home/yzy_jh/ref/MSigDB/c2.cp.kegg.v7.5.1.symbols.gmt')

names(HallMarkers)

exp1 <- AverageExpression(SP_shcc1)
GSEAres <- GSVA::gsva(expr = exp1$Spatial,types,method='ssgsea')
pheatmap::pheatmap(GSEAres,scale = 'row')
table(SP_shcc1$seurat_clusters,SP_shcc1$orig.ident)

p1 <- SpatialDimPlot(SP_shcc1,label = T)
p2 <- SpatialFeaturePlot(SP_shcc1,features = 'nCount_Spatial')

p1|p2

p1 <- SpatialDimPlot(SP_shcc2,label = T)
p2 <- SpatialFeaturePlot(SP_shcc2,features = 'nCount_Spatial')

p1|p2

p1 <- SpatialDimPlot(SP_shcc4,label = T)
p2 <- SpatialFeaturePlot(SP_shcc4,features = 'nCount_Spatial')
p3 <- FeaturePlot(SP_shcc4,features = 'nCount_Spatial')
p1|p2|p3

SpatialFeaturePlot(SP_shcc3,features = c('HIF1A','C3','GSTA1'))
SpatialFeaturePlot(SP_shcc3,features = c('IGHG4','IGLC2','IGLC1','IGHG1'))
SpatialFeaturePlot(SP_shcc3,features = c('FTL','APOE','SOD2','FTH1'))


Idents(SP_shcc1) <- 'orig.ident'
DefaultAssay(SP_shcc1) <- ''
diff_in_n <- FindAllMarkers(SP_shcc1,assay = "Spatial")
top10 <- diff_in_n %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
VlnPlot(SP_shcc1,features = unique(top10$gene),group.by = 'orig.ident',pt.size = 0)

Idents(SP_shcc4) <- 'orig.ident'
DefaultAssay(SP_shcc4) <- 'Spatial'
diff_patient4 <- FindAllMarkers(SP_shcc4,assay = "Spatial")
top10 <- diff_patient4 %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)

spata_obj <- read_rds('obj_list.rds')
names(spata_obj) <- samples
spata_obj <- parallel::mclapply(samples,function(x){
  GSEAres <- GSVA::gsva(expr = as.matrix(spata_obj[[x]]@data[[x]]$counts),list(HALLMARK_GLYCOLYSIS),method='ssgsea')
  spata_obj[[x]]@fdata[[x]]$GLYCOLYSIS <- as.numeric(GSEAres[1,])
  print(spata_obj[[x]]@fdata[[x]]$GLYCOLYSIS)
  return(spata_obj[[x]])
},mc.cores = 8)

fig_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/Integrative/SPATA2/GSEA'
checkdir(fig_dir)

ann_obj <- seurat_obj_shcc5T

HALLMARK_HYPOXIA <- HallMarke[['HALLMARK_HYPOXIA']]@geneIds
names(HALLMARK_HYPOXIA) <- 'HALLMARK_HYPOXIA'

HALLMARK_GLYCOLYSIS

HALLMARK_GLYCOLYSIS <- HallMarke[['HALLMARK_GLYCOLYSIS']]@geneIds
names(HALLMARK_GLYCOLYSIS) <- 'HALLMARK_GLYCOLYSIS'


HALLMARK_G2M_CHECKPOINT <- HallMarke[['HALLMARK_G2M_CHECKPOINT']]@geneIds
names(HALLMARK_G2M_CHECKPOINT) <- 'HALLMARK_G2M_CHECKPOINT'

HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY <- HallMarke[['HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY']]@geneIds
names(HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY) <- 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'

HALLMARK_OXIDATIVE_PHOSPHORYLATION <- HallMarke[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']]@geneIds
names(HALLMARK_OXIDATIVE_PHOSPHORYLATION) <- 'HALLMARK_OXIDATIVE_PHOSPHORYLATION'

types <- lapply(names(HallMarkers),function(x){
  type <- HallMarke[[x]]@geneIds 
  return(type)
})
names(types) <- names(HallMarkers)

GSEAres <- GSVA::gsva(expr = exp$Spatial,types,method='ssgsea')

all_cluster <- AverageExpression(st_data)

GSEAres <- GSVA::gsva(expr = all_cluster$Spatial,types,method='ssgsea')
pheatmap::pheatmap(GSEAres,scale = 'row')

table(st_data$seurat_clusters,st_data$orig.ident)

GSEAres <- GSVA::gsva(expr = as.matrix(ann_obj@assays$Spatial@counts),list(HALLMARK_GLYCOLYSIS),method='ssgsea')
ann_obj$GLYCOLYSIS <- as.numeric(GSEAres[1,])

GSEAres <- GSVA::gsva(expr = as.matrix(ann_obj@assays$Spatial@counts),list(HALLMARK_HYPOXIA),method='ssgsea')
ann_obj$HYPOXIA <- as.numeric(GSEAres[1,])

GSEAres <- GSVA::gsva(expr = as.matrix(ann_obj@assays$Spatial@counts),list(HALLMARK_G2M_CHECKPOINT),method='ssgsea')
ann_obj$G2M_CHECKPOINT <- as.numeric(GSEAres[1,])

GSEAres <- GSVA::gsva(expr = as.matrix(ann_obj@assays$Spatial@counts),list(HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY),method='ssgsea')
ann_obj$REACTIVE_OXYGEN_SPECIES_PATHWAY <- as.numeric(GSEAres[1,])

GSEAres <- GSVA::gsva(expr = as.matrix(ann_obj@assays$Spatial@counts),list(HALLMARK_OXIDATIVE_PHOSPHORYLATION),method='ssgsea')
ann_obj$OXIDATIVE_PHOSPHORYLATION <- as.numeric(GSEAres[1,])


SpatialFeaturePlot(ann_obj,features = 'GLYCOLYSIS') 
ggsave(filename = 'phenotype/GLYCOLYSIS.png',width = 5,height = 5)
SpatialFeaturePlot(ann_obj,features = 'HYPOXIA')
SpatialFeaturePlot(ann_obj,features = 'G2M_CHECKPOINT')
SpatialFeaturePlot(ann_obj,features = 'HYPOXIA')
SpatialFeaturePlot(ann_obj,features = 'REACTIVE_OXYGEN_SPECIES_PATHWAY')
SpatialFeaturePlot(ann_obj,features = 'OXIDATIVE_PHOSPHORYLATION')

exp <- AverageExpression(ann_obj)
exp$Spatial

st_data@meta.data %<>% select(-names(features))

marker_select <- de %>% group_by(seurat_clusters) %>% filter(avg_logFC>0.3)
marker_select_list <- split(marker_select,marker_select$seurat_clusters)

ids = bitr(de$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
markers_bitr = merge(de,ids,by.x='gene',by.y='SYMBOL')
markers_bitr <- split(markers_bitr$ENTREZID,marker_select$seurat_clusters)

KEGG_cluster <- compareCluster(markers_bitr, fun="enrichKEGG",
                               organism="hsa", pvalueCutoff=0.2)
p_KEGG <- dotplot(KEGG_cluster) + theme(axis.text.x = element_text(angle = 45, 
                                                                   vjust = 0.5, hjust=0.5,size=15),
                                        axis.text.y = element_text(size=10,face = 'bold')) + 
  ggtitle("KEGG")

ggsave(p_KEGG,filename = 'kegg.png',width = 8,height = 12)
GO_cluster_BP <- compareCluster(markers_bitr, fun="enrichGO",
                                OrgDb =org.Hs.eg.db, ont = "BP",pvalueCutoff=0.) 
GO_cluster_CC <- compareCluster(markers_bitr, fun="enrichGO",
                                OrgDb =org.Hs.eg.db, ont = "CC",pvalueCutoff=0.1)
GO_cluster_MF <- compareCluster(markers_bitr, fun="enrichGO",
                                OrgDb =org.Hs.eg.db, ont = "MF",pvalueCutoff=0.1)

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

lq_gene <- c('IFNG','IL1B','CXCL9','IL15','TNFSF123','CD40','PTPRC')#'TNF','HMGB1',
receiver_gene <- c('CCL3','CCL4','CCL5','ALB','AREG')
receiver_gene1 <- c('KLRFK1','JUND','LTB','CD69','FTL','DDIT4','NKFBIA','XCL1')
lq_gene <- c('C3','SAA1','APOA1')
receiver_gene <- c('CCL20','ALB','APOE')
SpatialFeaturePlot(ann_obj,features = lq_gene)
SpatialFeaturePlot(ann_obj,features = receiver_gene)
SpatialFeaturePlot(ann_obj,features = receiver_gene1)

SpatialFeaturePlot(ann_obj,features =c('GSTA2','LDHA'))

marker_select_list$`4`$gene
ids = bitr(marker_select_list$`4`$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
markers_bitr = merge(marker_select_list$`4`,ids,by.x='gene',by.y='SYMBOL')

y <- clusterProfiler::gseGO(,org.Hs.eg.db,
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)

head(y)

library(ReactomePA)
data(geneList, package="DOSE"
)
de <- names(geneList)[abs(geneList) > 1.5]
head(de)


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
                    pt_clrp = "npg",
                    pt_size = 1.5) +
  ggplot2::labs(color = "Clusters")
ggsave(p|p1|p2,filename = glue("{fig_dir}/{x}_HYPOXIA.png"),width = 15,height = 8)
ggsave(p3,filename = glue("{fig_dir}/{x}_cluster.png"),width = 8,height = 8)
},mc.cores = 8)

plots <- plotSurfaceInteractive(object = spata_obj2)

write_rds(spata_obj,file = 'spata_obj.rds')
p <- plotSurface(spata_obj2,
            color_by = "ALB",
            pt_size = 1.9) 
p

SpatialDimPlot(SP_shcc1)


exculde

clod_marker <- c('CXCR4','CXCL12','','','','')
hold_marker <- c('')
checkpoint <- c('PDCD1','CTLA4','LAG3','TIGIT')
SpatialFeaturePlot(SP_shcc1,features = 'PDCD1',alpha = c(0,1.5))
SpatialFeaturePlot(SP_shcc1,checkpoint[4],alpha = c(0,1.5))







