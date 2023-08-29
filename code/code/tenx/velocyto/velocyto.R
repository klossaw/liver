library(Seurat)
library(velocyto.R)
library(tidyverse)
library(SeuratWrappers)
##数据基础分析
# 读取loom文件

load("projects/liver/analysis/qzhang/human/tenx/tmp.Rdata")
loom.dir <- list.files("projects/liver/analysis/qzhang/human/tenx/velocity/loom",full.names = TRUE)


spliced.list <- lapply(loom.dir, function(x){
  loom <- read.loom.matrices(x,engine = "hdf5r")
  return(loom)
  })



spliced <- cbind(spliced.list[[1]])

sce[["spliced"]] <- CreateAssayObject(counts = spliced.list[["spliced"]])#object为Seurat object

object[["unspliced"]] <- CreateAssayObject(counts = x[["unspliced"]])

object[["ambiguous"]] <- CreateAssayObject(counts = x[["ambiguous"]]) 


velo <- as.Seurat(x = velo)
# 降维聚类
velo <- velo %>% SCTransform(assay="spliced") %>% RunPCA(verbose=F) 
ElbowPlot(velo, ndims = 50)
nPC=1:20
velo <- FindNeighbors(velo, dims = nPC) %>% FindClusters() %>% 
  RunUMAP(dims = nPC) %>% RunTSNE(dims = nPC)

##给细胞分配颜色
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = velo)))
names(x = ident.colors) <- levels(x = velo)
cell.colors <- ident.colors[Idents(object = velo)]
names(x = cell.colors) <- colnames(x = velo)

##速率分析
velo <- RunVelocity(velo, deltaT = 1, kCells = 25, fit.quantile = 0.02, 
                    spliced.average = 0.2, unspliced.average = 0.05, ncores = 18)
#kCells：用于斜率平滑度计算最近邻细胞数量，越大越平滑，越小越能反映局部特征
#fit.quantile：0.02代表对基因表达量最高2%与最低2%的值执行gamma拟合
#spliced.average：过滤低表达丰度基因的标准，计算的是基因在cluster内的平均counts值
#unspliced.average：同上

##全局速率可视化
emb = Embeddings(velo, reduction = "umap")
vel = Tool(velo, slot = "RunVelocity")
show.velocity.on.embedding.cor(emb = emb, vel = vel, n = 200, scale = "sqrt", 
                               cell.colors = ac(cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, 
                               show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
                               arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)

##特定基因速率可视化
gene = "Camp"
RunVelocity(velo, deltaT=1, kCells=25, fit.quantile=0.02, old.fit=vel, 
            cell.emb=emb, cell.colors=cell.colors, show.gene=gene, do.par=T)

