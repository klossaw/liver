library(scales)
rna_velocity_dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/velocity/rna_velocity/"

setwd(file.path(project_work_dir, "rna_velocity"))

seurat_names <- substr(colnames(sce_addlabel),7,28) 

seurat_obj <- sce_addlabel
row.names(seurat_obj@meta.data) <- seurat_names
row.names(seurat_obj@reductions$umap@cell.embeddings) <- seurat_names

# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
write.csv(Embeddings(seurat_obj, reduction = "umap"), file = "cell_embeddings.csv")
# 获取每个细胞的barcode
write.csv(seurat_names, file = "cellID_obs.csv", row.names = FALSE)
# 提取每个细胞的cluster信息
write.csv(seurat_obj@meta.data[, 'seurat_clusters', drop = FALSE], file = "cell_clusters.csv")
# 提取每个细胞的celltype信息
write.csv(seurat_obj@meta.data[, 'cell_types', drop = FALSE], file = "cell_celltype.csv")
# 获取celltype的颜色信息
hue_pal()(length(levels(seurat_obj$cellType)))
# 获取cluster的颜色信息
hue_pal()(length(levels(seurat_obj$seurat_clusters)))

# 绘制umap图，与RNA速率图对比看
Idents(seurat_obj) <- "cell_types"
alpha.use <- 0.8
p1 <- DimPlot(object = seurat_obj, reduction = "umap", label = TRUE, label.size = 5, pt.size=0.4, raster=FALSE)
p1$layers[[1]]$mapping$alpha <- alpha.use
