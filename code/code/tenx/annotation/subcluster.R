#packages
library(Seurat)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(glue)
library(jhtools)
#dir
project_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/06subcluster/'
setwd(project_dir)
subT_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/06subcluster/subT'
checkdir(subT_dir)
subMyeloid_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/06subcluster/subMyeloid'
checkdir(subMyeloid_dir)
subHepatocytes_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/06subcluster/subHepatocytes'
checkdir(subHepatocytes_dir)
#load data
load('ALL.Rdata')
sce <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/05AddLabel/sce_addlabel.rds')
#analysis
Idents(sce) <- 'cell_types'
sub_NK_T <- subset(sce,ident='NK&T cell')

sub_T <- SCTransform(sub_NK_T)
#filter mt rp
sub_NK_T <- FindVariableFeatures(sub_NK_T, selection.method = "vst", nfeatures = 2000)
Vargene <- VariableFeatures(sub_NK_T)
Vargene_filter <- Vargene[-c(grep('^RP[SL]',Vargene),grep('^MT-',Vargene))] %>% as.character()
#cluster

sub_NK_T <- RunPCA(sub_NK_T,features =Vargene_filter,npcs = 50) %>% 
  RunHarmony(assay.use = "SCT",group.by.vars = 'orig.ident') %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = c(0.5,0.3), verbose = FALSE)

DimPlot(sub_NK_T,pt.size = 0.5,label = T)

filtergene <- c(grep('^RP[SL]',row.names(sub_NK_T)),grep('^MT-',row.names(sub_NK_T)))
NK_T_marker_filter <- FindAllMarkers(sub_NK_T[-filtergene,],assay = "RNA")
DefaultAssay(sub_NK_T) <- 'RNA'
NK_T_marker <- NormalizeData(sub_NK_T) %>% FindAllMarkers(assay = "RNA")

top10 <- NK_T_marker %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
top20 <- NK_T_marker %>% dplyr::filter(pct.1 > 0.4 & pct.2 < 0.4) %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
top30 <- NK_T_marker %>% dplyr::filter(pct.1 > 0.4 & pct.2 < 0.4) %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
top50 <- NK_T_marker %>% dplyr::filter(pct.1 > 0.4 & pct.2 < 0.4) %>% group_by(cluster) %>% top_n(50,wt = avg_log2FC)


sub_NK_T_exp <- AverageExpression(sub_NK_T,assays = 'RNA',slot = "scale.data")$RNA
annotation_marker <- c(Tmarker,naiveTcells,exhaustedTcells,regulatoryTcells,Cytokines,resident,TF,Thelper)
pheatmap::pheatmap(t_exp[annotation_marker,],scale = 'row',cluster_rows = F,cluster_cols = F)

NKT_label <- data.frame('cluster'=c('0','1','3','4','5','6','7','9','10','11','12','2','8'),
                        'cell_types'=c('CD4-S100A4','CD4-CCR7','NK-GNLY','CD8-CCL4','CD8-HSPA6','NK-XCL1',
                                       'CD4-FOXP3','CD8-IFNG','CD8-CCL20','NK-GZMB','CD4-LAG3','other','other'))
NKT_label <- NKT_label[order(NKT_label$cell_types),]
sub_NK_T <- sce_add_label(sub_NK_T,label = NKT_label)
DimPlot(sub_NK_T,group.by = 'cell_types',label = T) %>% ggsave(filename = '06subcluster/subT_ann.png',width = 10,height = 10)

for (i in unique(sub_NK_T$cell_types)) {
  cells <- row.names(sub_NK_T@meta.data)[sub_NK_T$cell_types==i]
  sce@meta.data[cells,]$cell_subtypes <- i
}
FeaturePlot(sub_NK_T,features = Tmarker)|DimPlot(sub_NK_T,label = T)
DotPlot(sub_NK_T,features = Tmarker)|DimPlot(sub_NK_T,label = T)

library(ComplexHeatmap)
Heatmap(sub_NK_T_exp[annotation_marker,])

Tmarker <- c('CD4','CD8A','CD8B','CD3D','CD3E','CD3G')
naiveTcells <- c('CCR7', 'LEF1', 'SELL', 'TCF7', 'S1PR1')
exhaustedTcells <- c('HAVCR2', 'PDCD1', 'GZMB', 'ITGAE','LAG3')
regulatoryTcells <- c('FOXP3', 'CTLA4', 'IL2RA')#5
Cytokines <- c('IL2','IL17A','GNLY','GZMK','GZMB','GZMA','NKG7','IFNG','LAMTOR3','PRF1')
resident <- c('CD69','RUNX3','NR4A1')
co_stimulatory <- c('TNFRSF14','CD28','ICOS','TNFRSF9')
TF <- c('TBX21','ZNF683','ZEB2','ID2','EOMES','HIF1A','TOX')
Thelper <- c('IL17A','NFKBIA','CD40LG')
nk_marker <- c('NKG7','KLRF1','KLRD1')


allT_marker <- c(Tmarker,naiveTcells,exhaustedTcells,regulatoryTcells,Cytokines,co_stimulatory,TF,Thelper,nk_marker) %>% unique()

DotPlot(sub_NK_T,features = unique(allT_marker)) + coord_flip() + 
  scale_color_viridis_c() + 
  RotatedAxis() + 
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0,hjust = 1,vjust = 0.5))

pheatmap::pheatmap(sub_NK_T_exp[allT_marker,],scale = 'row',cluster_rows = F,cluster_cols = T)
DimPlot(sce,group.by = 'cell_subtypes',label = T)
DimPlot(sub_Hepatocytes,split.by = 'group',label = TRUE)
DimPlot(sub_Hepatocytes,split.by = 'orig.ident',label = TRUE)
DimPlot(sub_Hepatocytes,group.by = 'seurat_clusters',label = TRUE)
DimPlot(sub_Hepatocytes,group.by = 'orig.ident',label = TRUE)
################################## sub Hepatocytes #############################
sce <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/05AddLabel/sce_addlabel.rds')
Idents(sce) <- 'cell_types'
sub_Hepatocytes <- subset(sce,ident='Hepatocytes')

DefaultAssay(sub_Hepatocytes) <- 'RNA'
sub_Hepatocytes@assays$SCT <- NULL
sub_Hepatocytes <- SCTransform(sub_Hepatocytes)
#filter mt rp
Vargene <- VariableFeatures(sub_Hepatocytes)
Vargene_filter <- Vargene[-c(grep('^RP[SL]',Vargene),grep('^MT-',Vargene))] %>% as.character()
sub_Hepatocytes <- FindVariableFeatures(sub_Hepatocytes, selection.method = "vst", nfeatures = 1500)
#cluster
sub_Hepatocytes <- RunPCA(sub_Hepatocytes,features =Vargene_filter,npcs = 50) %>% 
  RunHarmony(assay.use = "SCT",group.by.vars = 'orig.ident') %>% 
  RunUMAP(reduction = "harmony", dims = dims) %>% 
  RunTSNE(reduction = "harmony", dims = dims) %>% 
  FindNeighbors(reduction = "harmony", dims = dims) %>% 
  FindClusters(resolution = 0.2, verbose = FALSE)


sub_Hepatocytes_marker <- NormalizeData(sub_Hepatocytes) %>% FindAllMarkers(assay = "RNA")

top10 <- sub_Hepatocytes_marker %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)

DefaultAssay(sub_Hepatocytes) <- 'SCT'
DoHeatmap(sub_Hepatocytes,features = top10$gene)

sub_Hepatocytes@meta.data$cell_types <- paste('Hepatocytes-C',sub_Hepatocytes@meta.data$seurat_clusters,sep = '')

for (i in unique(sub_Hepatocytes$cell_types)) {
  cells <- row.names(sub_Hepatocytes@meta.data)[sub_Hepatocytes$cell_types==i]
  sce@meta.data[cells,]$cell_subtypes <- i
}

DimPlot(sub_Hepatocytes,group.by = 'cell_types',label = T) %>% 
  ggsave(filename = '06subcluster/sub_Hepatocytes_umap.png',width = 10,height = 10)

DimPlot(sce,group.by = 'cell_subtypes')

Normal <- c("hcc1C","hcc3C","hcc4C")
Tumor <- c("hcc1T","hcc2T","hcc3T","hcc4T","hcc5T") 
output_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/infercnv/Hepatocytes' 
run_infercnv(object = sub_Hepatocytes,ref = 'malignant',output_dir = output_dir,Normal = Normal,Tumor = Tumor)

sub_Hepatocytes_cds <- monocle2_order(sub_Hepatocytes_harmony)
output_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/06subcluster/sub_Hepatocytes'
monocle2_plot(monocle_cds = sub_Hepatocytes_cds,output_dir = output_dir,genes = unique(top10$gene))

#/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/infercnv/Hepatocytes
###################################### sub Myeloid cell ########################
sub_Myeloid <- subset(sce,ident=c('Myeloid cell'))

DefaultAssay(sub_Myeloid) <- 'SCT'
sub_Myeloid <- SCTransform(sub_Myeloid,variable.features.n = 1500) %>% 
  RunPCA() %>% 
  RunHarmony(assay.use = "SCT",group.by.vars = 'orig.ident') %>% 
  RunUMAP(reduction = "harmony", dims = 1:15) %>% 
  RunTSNE(reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.1, verbose = FALSE)

DimPlot(sub_Myeloid,label = T)

sub_Myeloid_marker <- NormalizeData(sub_Myeloid,assay = 'RNA') %>% FindAllMarkers(assay = "RNA")

top10 <- sub_Myeloid_marker %>% dplyr::filter(pct.1 > 0.4 & pct.2 < 0.4) %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
top20 <- sub_Myeloid_marker %>% dplyr::filter(pct.1 > 0.4 & pct.2 < 0.4) %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
top30 <- sub_Myeloid_marker %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
top50 <- sub_Myeloid_marker %>% group_by(cluster) %>% top_n(50,wt = avg_log2FC)

Macro_marker <- c('CST3', 'LYZ', 'CD68', 'CD163', 'CD169')
Mono_marker <- c('CST3', 'LYZ', 'CD14', 'VCAN', 'S100A9', 'RNASE2', 'S100A12')
mast_marker <- c('CST3', 'S100A9', 'AIF1', 'CD79A', 'IL7R')
write_csv(sub_Myeloid_marker,file = glue("{subMyeloid_dir}/marker.csv"))
#plot
Myeloid_marker <- c('CST3', 'LYZ','FCGR3A','CD16','RNF213','S100A8','BAG3','C1QC','C1QB','IL32','CCL5','LAMP3','IDO1')

sub_Myeloid_label <- data.frame(cluster=c('0','1','2','3','4'),
                                cell_types=c('Mono-RNF213','Macro-CXCL8','Macro-C1QC','Mono-IL32','DC-LAMP3'))
sub_Myeloid_label <- sub_Myeloid_label[order(sub_Myeloid_label$cell_types),]
sub_Myeloid <- sce_add_label(sub_Myeloid,label = sub_Myeloid_label)

for (i in unique(sub_Myeloid$cell_types)) {
  cells <- row.names(sub_Myeloid@meta.data)[sub_Myeloid$cell_types==i]
  sce@meta.data[cells,]$cell_subtypes <- i
}

p <- DimPlot(sub_Myeloid,group.by = 'cell_types',label = T)
p
DimPlot(sce,group.by = 'cell_subtypes',label = T)

ggsave(p,filename = glue("06subcluster/subMyeloid_label.png"),width = 10,height = 8)

p <- DotPlot(sub_Myeloid,features = Myeloid_marker,group.by = 'cell_types')+ coord_flip() + 
     scale_color_viridis_c() + 
     RotatedAxis() + 
     theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 0,hjust = 1,vjust = 0.5))
p
ggsave(p,filename = glue("{subMyeloid_dir}/marker.png"),width = 10,height = 10)

sub_Myeloid <- ScaleData(sub_Myeloid)
p <- DoHeatmap(object = sub_Myeloid,assay = 'RNA',slot = 'data',features = top10$gene,group.by = 'seurat_clusters') 
p
ggsave(p,filename = glue("{subMyeloid_dir}/label.png"),width = 10,height = 10)
write_rds(sub_Myeloid,file = glue("{subMyeloid_dir}/sub_Myeloid.rds"))

save(sce,sub_NK_T,sub_Myeloid,sub_Hepatocytes,NK_T_marker,file = '06subcluster/ALL.Rdata')

library(RColorBrewer)
colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", 
            "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", 
            colorRampPalette(brewer.pal(8,'Dark2'), alpha = 1)(8))
meta.data <- sub_NK_T@meta.data
cluster.color <- colors[1:length(unique(meta.data$cell_types))]
p1 <- ggplot(data = meta.data,mapping = aes(x = group, fill=cell_types))+
  geom_bar(stat = "count",width = 0.7, position = 'fill') +
  scale_fill_manual(values = cluster.color) +
  labs(y = "proportion of cells", x = "") +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(family = "serif", size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line.x = element_line (colour = "black"), axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = NA, color = NA), legend.position = "left") 
p1
p2 <- DimPlot(sub_NK_T,group.by = 'cell_types',cols = cluster.color,label = T,label.size = 3)
p2
ggsave(p1,filename = glue("06subcluster/subTprop.png"), width = 4, height = 4)
ggsave(p2,filename = glue("06subcluster/subTumap.png"), width = 6, height = 5)

meta.data <- sub_Myeloid@meta.data
cluster.color <- colors[1:length(unique(meta.data$cell_types))]
p1 <- ggplot(data = meta.data,mapping = aes(x = group, fill=cell_types))+
  geom_bar(stat = "count",width = 0.7, position = 'fill') +
  scale_fill_manual(values = cluster.color) +
  labs(y = "proportion of cells", x = "") +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(family = "serif", size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line.x = element_line (colour = "black"), axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = NA, color = NA), legend.position = "left") 
p1
p2 <- DimPlot(sub_Myeloid,group.by = 'cell_types',cols = cluster.color,label = T)
p1/p2
ggsave(p1,filename = glue("06subcluster/subMeyloidprop.png"), width = 4, height = 4)
ggsave(p2,filename = glue("06subcluster/subMeyloidumap.png"), width = 6, height = 5)
##
meta.data <- sub_Hepatocytes@meta.data
cluster.color <- colors[1:length(unique(meta.data$cell_types))]
p1 <- ggplot(data = meta.data,mapping = aes(x = group, fill=cell_types))+
  geom_bar(stat = "count",width = 0.7, position = 'fill') +
  scale_fill_manual(values = cluster.color) +
  labs(y = "proportion of cells", x = "") +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(family = "serif", size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line.x = element_line (colour = "black"), axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = NA, color = NA), legend.position = "left") 
p1
p2 <- DimPlot(sub_Hepatocytes,group.by = 'cell_types',cols = cluster.color,label = F)
p1/p2
ggsave(p1,filename = glue("06subcluster/Hepatocytesprop.png"), width = 4, height = 4)
ggsave(p2,filename = glue("06subcluster/Hepatocytesumap.png"), width = 6, height = 5)

library(monocle)
detach("package:monocle3",unload = TRUE)
sub_Myeloid_cds <- monocle2_order(sub_Myeloid)
sub_NK_T_cds <- monocle2_order(sub_NK_T[,sample(1:31509,10000)])
sub_Hepatocytes_cds <- monocle2_order(sub_Hepatocytes)

monocle::plot_cell_trajectory(sub_Myeloid_cds,color_by = 'cell_types')
monocle2_plot(monocle_cds = sub_Myeloid_cds,output_dir = '06subcluster/subMyeloid/',genes = unique(sub_Myeloid_marker$gene))


monocle::plot_cell_trajectory(sub_Myeloid_cds,color_by = 'cell_types')
monocle2_plot(monocle_cds = sub_Hepatocytes_cds,output_dir = '06subcluster/sub_Hepatocytes',genes = unique(sub_Hepatocytes_marker$gene))



monocle::plot_cell_trajectory(sub_Myeloid_cds,color_by = 'cell_types')
monocle2_plot(monocle_cds = sub_NK_T_cds,output_dir = '06subcluster/subT/',genes = unique(NK_T_marker$gene))
