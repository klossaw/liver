dir <- paste(list.files("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/counts",full.names = TRUE),"/outs",sep = "")
names(dir) = c("shcc1C", "shcc1T", "shcc2T", "shcc3C", "shcc3T", "shcc4C", "shcc4T", "shcc5T")
liver <- list()
for(i in 1:length(dir)){
  liver[[i]] <-Seurat::Load10X_Spatial(data.dir = dir[i])
  liver[[i]]@meta.data$orig.ident <-names(dir)[i]
}
for (i in 1:length(liver)) {
  liver[[i]] <- SCTransform(liver[[i]], assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
  } ##SCT标准化

liver.merge <- merge(liver[[1]], y=c(liver[[2]], liver[[3]], liver[[4]],liver[[5]],liver[[6]],liver[[7]],liver[[8]]))
dim(liver.merge)
table(liver.merge@meta.data$orig.ident)
DefaultAssay(liver.merge) <- "SCT"
VariableFeatures(liver.merge) <- c(VariableFeatures(liver[[1]]), 
                                   VariableFeatures(liver[[2]]), 
                                   VariableFeatures(liver[[3]]), 
                                   VariableFeatures(liver[[4]]),
                                   VariableFeatures(liver[[5]]),
                                   VariableFeatures(liver[[6]]),
                                   VariableFeatures(liver[[7]]),
                                   VariableFeatures(liver[[8]])) 

liver.merge <- RunPCA(liver.merge, verbose = FALSE)
liver.merge <- FindNeighbors(liver.merge, dims = 1:30)
liver.merge <- FindClusters(liver.merge,resolution = 0.8, verbose = FALSE)
liver.merge <- RunUMAP(liver.merge, dims = 1:30)
pdf("merge_umap.pdf",width=10,height=5)
DimPlot(liver.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()
pdf("merge_umap_slice.pdf",width=10,height=5)
SpatialDimPlot(liver.merge) 
dev.off()
pdf("merge_umap_gene.pdf",width=10,height=5)
SpatialFeaturePlot(liver.merge, features = c("", ""))
dev.off()

saveRDS(liver.merge,"merge.rds") #保存数据
liver.markers <- FindAllMarkers(liver.merge, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
write.table(liver.markers,"marker.txt",row.names=TRUE,col.names=TRUE,sep="\t")
topgene<-liver.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
pdf("merge_heatmap.pdf",width=10,height=8)
DoHeatmap(liver.merge, features = topgene$gene,size = 2) + NoLegend()
dev.off()

for (i in 1:length(liver)) {
  liver[[i]] <- NormalizeData(liver[[i]]) #未使用SCT标准化
  liver[[i]] <- FindVariableFeatures(liver[[i]], selection.method = "vst", nfeatures = 2000)}
##
features <- SelectIntegrationFeatures(object.list = liver)
liver.anchors <- FindIntegrationAnchors(object.list = liver, anchor.features = features)
liver_in <- IntegrateData(anchorset = liver.anchors)
dim(liver_in)
DefaultAssay(liver_in) <- "integrated"
liver_in <- ScaleData(liver_in, verbose = FALSE)
liver_in <- RunPCA(liver_in, npcs = 30, verbose = FALSE)
liver_in <- RunUMAP(liver_in, reduction = "pca", dims = 1:30)
liver_in <- FindNeighbors(liver_in, reduction = "pca", dims = 1:30)
liver_in <- FindClusters(liver_in, resolution = 0.8,verbose = FALSE)
pdf("merge_umap1.pdf",width=10,height=5)
DimPlot(liver_in, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()
pdf("merge_umap_slice1.pdf",width=10,height=5)
SpatialDimPlot(liver_in)
dev.off()

for (i in 1:length(liver)) {
  liver[[i]] <- SCTransform(liver[[i]], assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)} #使用SCT标准化
features <- SelectIntegrationFeatures(object.list = liver, nfeatures = 3000)
liver <- PrepSCTIntegration(object.list = liver, anchor.features = features)
liver.anchors <- FindIntegrationAnchors(object.list = liver, normalization.method = "SCT", anchor.features = features)
liver.sct <- IntegrateData(anchorset = liver.anchors, normalization.method = "SCT")
dir.create("Integrate_SCT_data")
setwd("Integrate_SCT_data")
dim(liver.sct)
liver.sct <- RunPCA(liver.sct, npcs = 30, verbose = FALSE)
liver.sct <- RunUMAP(liver.sct, reduction = "pca", dims = 1:30)
liver.sct <- FindNeighbors(liver.sct, reduction = "pca", dims = 1:30)
liver.sct <- FindClusters(liver.sct, resolution = 0.8,verbose = FALSE)
pdf("merge_umap2.pdf",width=10,height=5)
DimPlot(liver.sct, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()
pdf("merge_umap_slice2.pdf",width=10,height=5)
SpatialDimPlot(liver.sct)
dev.off()
save(list = ls,file = "merge.Rdata")

setwd("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/")
liver.merge <- read_rds("visium/seurat/merge/merge.rds")
SpatialDimPlot(liver.merge,group.by = "seurat_clusters")
SpatialPlot(liver.merge,images = c("slice1.3","slice1.4"))
DimPlot(liver.merge,reduction = "umap",split.by = "orig.ident")

de_marker <- FindAllMarkers(liver.merge)
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
sce_reference <- singleRdata
DefaultAssay(sce_reference) <- "SCT"
DefaultAssay(liver.merge) <- "SCT"

anchors <- FindTransferAnchors(reference = sce_reference, 
                               query = liver.merge, 
                               normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = sce_reference$singleR.pruned.labels, 
                                  prediction.assay = TRUE,
                                  weight.reduction = liver.merge[["pca"]], 
                                  dims = 1:30)

liver.merge[["predictions"]] <- predictions.assay

DefaultAssay(liver.merge) <- "predictions"

SpatialFeaturePlot(liver.merge, 
                   features = c("Smooth-muscle-cells", "Hepatocytes"), 
                   pt.size.factor = 1.6,
                   images = c("slice1.2","slice1.3"), 
                   ncol = 2, 
                   crop = TRUE
                   )

liver.merge <- FindSpatiallyVariableFeatures(liver.merge, 
                                             assay = "predictions", 
                                             selection.method = "markvariogram",
                                             features = rownames(liver.merge), 
                                             r.metric = 5, 
                                             slot = "data"
                                             )

top.clusters <- head(SpatiallyVariableFeatures(liver.merge), 4)

SpatialPlot(object = liver.merge, 
            features = top.clusters, 
            ncol = 2)

SpatialFeaturePlot(liver.merge, 
                   features = c("5","6"), 
                   pt.size.factor = 1, 
                   ncol = 2, 
                   crop = FALSE, 
                   alpha = c(0.1, 1)
                   )
