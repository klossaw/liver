work.dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/pyscenic/"
mkdir <- function(dir.list){
  lapply(dir.list,function(x){
    if (!dir.exists(x)) {
      dir.create(x)
    }
  })
}
setwd(work.dir)
file.path(work.dir,c("data","code","plot")) %>% mkdir() 

suppressMessages(library(ggplot2))
suppressMessages(library(SCENIC))
suppressMessages(library(SCopeLoomR))
suppressMessages(library(AUCell))

od <- file.path(work.dir,"plot/")
projectname <- "liver"

# function
P.plot <- function(data,od,filename,width = 5, height = 15 ){
  file.pdf <- paste(filename,'.pdf',sep = "")
  file.png <- paste(filename,'.png',sep = "")
  ggsave(data,filename = paste(od, file.pdf,sep = ""), width = width, height = height)
  ggsave(data,filename = paste(od, file.png,sep = ""), width = width, height = height)
}

cellInfo <- sce@meta.data
rownames(cellInfo)<-gsub('-','.',rownames(cellInfo))#将sce和auc两个文件基因名格式统一

scenicLoomPath="data/sample_SCENIC.loom"
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 


if(!is.null(sce@meta.data$seurat_clusters)){
  regulonActivity_bycluster <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                      function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  colnames(regulonActivity_bycluster)<-factor(colnames(regulonActivity_bycluster))
  
  # plot heatmap 注意转录因子数量一般情况是200到300，如果少于100则修改参数
  p.all <- pheatmap::pheatmap(regulonActivity_bycluster,
                              scale ='row',
                              cluster_cols = F,
                              show_rownames = T, #fontsize_row=3, 
                              color=colorRampPalette(c("blue","white","red"))(100),
                              breaks=seq(-3, 3, length.out = 100),
                              treeheight_row=10, treeheight_col=10,
                              border_color=NA)
  p.50 <- pheatmap::pheatmap(regulonActivity_bycluster[sample(1:length(row.names(regulonActivity_bycluster)),50,replace = F),],
                             scale ='row',
                             cluster_cols = F,
                             show_rownames = T, #fontsize_row=3, 
                             color=colorRampPalette(c("blue","white","red"))(100),
                             breaks=seq(-3, 3, length.out = 100),
                             treeheight_row=10, treeheight_col=10,
                             border_color=NA)
  p.100 <- pheatmap::pheatmap(regulonActivity_bycluster[sample(1:length(row.names(regulonActivity_bycluster)),100,replace = F),],
                              scale ='row',
                              cluster_cols = F,
                              show_rownames = T, #fontsize_row=3, 
                              color=colorRampPalette(c("blue","white","red"))(100),
                              breaks=seq(-3, 3, length.out = 100),
                              treeheight_row=10, treeheight_col=10,
                              border_color=NA)
  # save image 
  filename = paste(projectname,'regulonAUC_50',sep = '_')
  P.plot(data = p.50,od=od,filename = filename,width =5 ,height =15)
  
  filename = paste(projectname,'regulonAUC_100',sep = '_')
  P.plot(data = p.100,od=od,filename = filename,width =5 ,height =20 )
  
  filename = paste(projectname,'regulonAUC_all',sep = '_')#转录因子总数量不确定，200左右 width：5,height：30合适。
  P.plot(data = p.all,od=od,filename = filename,width =5 ,height =20 )
  
}

if(!is.null(sce@meta.data$celltype)){
  regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                       function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  colnames(regulonActivity_byCellType)<-factor(colnames(regulonActivity_byCellType))
  # plot heatmap 注意转录因子数量一般情况是200到300，如果少于100则修改参数
  p.all <- pheatmap::pheatmap(regulonActivity_byCellType,
                              scale ='row',
                              cluster_cols = F,
                              show_rownames = T, #fontsize_row=3, 
                              color=colorRampPalette(c("blue","white","red"))(100),
                              breaks=seq(-3, 3, length.out = 100),
                              treeheight_row=10, treeheight_col=10,
                              border_color=NA)
  p.50 <- pheatmap::pheatmap(regulonActivity_byCellType[sample(1:length(row.names(regulonActivity_byCellType)),50,replace = F),],
                             scale ='row',
                             cluster_cols = F,
                             show_rownames = T, #fontsize_row=3, 
                             color=colorRampPalette(c("blue","white","red"))(100),
                             breaks=seq(-3, 3, length.out = 100),
                             treeheight_row=10, treeheight_col=10,
                             border_color=NA)
  p.100 <- pheatmap::pheatmap(regulonActivity_byCellType[sample(1:length(row.names(regulonActivity_byCellType)),100,replace = F),],
                              scale ='row',
                              cluster_cols = F,
                              show_rownames = T, #fontsize_row=3, 
                              color=colorRampPalette(c("blue","white","red"))(100),
                              breaks=seq(-3, 3, length.out = 100),
                              treeheight_row=10, treeheight_col=10,
                              border_color=NA)
  # save image
  filename = paste(projectname,'regulonAUC_50',sep = '_')
  P.plot(data = p.50,od=od,filename = filename)
  
  filename = paste(projectname,'regulonAUC_100',sep = '_')
  P.plot(data = p.100,od=od,filename = filename,width =5 ,height =20 )
  
  filename = paste(projectname,'regulonAUC_all',sep = '_')#转录因子总数量不确定，200左右 width：5,height：30合适。
  P.plot(data = p.all,od=od,filename = filename,width =5 ,height =20 )
}



