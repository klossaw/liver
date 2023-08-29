library(Seurat)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
sce <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/05AddLabel/sce_addlabel.rds')
Idents(sce) <- 'cell_types'
sub_T <- subset(sce,ident='Tcell')
DefaultAssay(sub_T) <- 'RNA'
sub_T@assays$SCT <- NULL
sub_T <- SCTransform(sub_T)
#filter mt rp
Vargene <- VariableFeatures(sub_T)
Vargene_filter <- Vargene[-c(grep('^RP[SL]',Vargene),grep('^MT-',Vargene))] %>% as.character()
sub_T <- FindVariableFeatures(sub_T, selection.method = "vst", nfeatures = 2000)
#cluster
dims <- 1:15
resolution <- c(0.4,0.3,0.2) 

sub_T <- RunPCA(sub_T,features =Vargene_filter,npcs = 30) %>% 
  RunHarmony(assay.use = "SCT",group.by.vars = 'orig.ident') %>% 
  RunUMAP(reduction = "harmony", dims = dims) %>% 
  RunTSNE(reduction = "harmony", dims = dims) %>% 
  FindNeighbors(reduction = "harmony", dims = dims) %>% 
  FindClusters(resolution = resolution, verbose = FALSE)

filtergene <- c(grep('^RP[SL]',row.names(sub_T)),grep('^MT-',row.names(sub_T)))
T_marker <- FindAllMarkers(sub_T[-filtergene,],assay = "RNA")

top10 <- T_marker %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
top20 <- T_marker %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
top30 <- T_marker %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)

top50 <- T_marker %>% group_by(cluster) %>% top_n(50,wt = avg_log2FC)

sub_T@meta.data$subtypes <- paste('Tcell-C',sub_T@meta.data$seurat_clusters,sep = '')
sce@meta.data$cell_subtypes <- as.character(sce@meta.data$cell_types)
for (i in unique(sub_T$subtypes)) {
  cells <- row.names(sub_T@meta.data)[sub_T$subtypes==i]
  sce@meta.data[cells,]$cell_subtypes <- i
}

pal <- paletteer::paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5,2,6,8,10)]
pal <- paletteer::paletteer_d("ggsci::nrc_npg")
DimPlot(sub_T,group.by = 'seurat_clusters',cols = pal,label = T)
table(sub_T$orig.ident,sub_T$seurat_clusters)

t_exp <- AverageExpression(sub_T)$RNA
Idents(sub_T) <- 'seurat_clusters'

annotation_marker <- c(Tmarker,naiveTcells,exhaustedTcells,regulatoryTcells,Cytokines,resident,TF,Thelper)
pheatmap::pheatmap(t_exp[annotation_marker,],scale = 'row',cluster_rows = F)

Tmarker <- c('CD4','CD8A','CD8B','CD3D','CD3E','CD3G')
naiveTcells <- c('CCR7', 'LEF1', 'SELL', 'TCF7', 'S1PR1')
exhaustedTcells <- c('HAVCR2', 'PDCD1', 'GZMB', 'ITGAE','LAG3')
regulatoryTcells <- c('FOXP3', 'CTLA4', 'IL2RA')#5
Cytokines <- c('IL2','IL17A','GNLY','GZMK','GZMB','GZMA','NKG7','IFNG','LAMTOR3','PRF1')
resident <- c('CD69','RUNX3','NR4A1')
co_stimulatory <- c('TNFRSF14','CD28','ICOS','TNFRSF9')
TF <- c('TBX21','ZNF683','ZEB2','ID2','EOMES','HIF1A','TOX')
Thelper <- c('IL17A','NFKBIA','CD40LG')

select_sub_T <- subset(sub_T_,downsample=500)
DoHeatmap(object = select_sub_T,features = top10$gene)
VlnPlot(sub_T,assay = 'RNA',features = naiveTcells,pt.size = 0)

proliferative_Tcells <- c('MKI67', 'STMN1', 'CMC1')
central_memoryTcells <- c('CD8', 'IL7R', 'CD74', 'TYROBP')
effector_Tcells <- c('CX3CR1', 'KLRG1', 'FCGR3A', 'FGFBP2', 'S1PR1', 'GZMH')#3
central_memory_Tcells <- c('CCR7', 'SELL', 'ANXA1', 'ANXA2', 'S1PR1')
effector_memory_Tcells <- c('GZMK', 'CCL4', 'CCL5', 'NKG7')
tissue_resident_memory_T_cells <- c('KLRD1', 'KLRF1', 'GNLY', 'IL32')
MAIT <- c('CD8A','SLC4A10', 'KLRB1', 'ZBTB16', 'NCR3')
NK
CD56_high_NKcells <- c('NKG7', 'SELL', 'NCAM1', 'KLRD1', 'GNLY', 'CCL3')
circulating_NK_cells_n <- c('NKG7', 'FCGR3A', 'IFNG', 'GNLY', 'GZMB', 'CD69')
circulating_NK_cells_t <- c('NKG7', 'FCGR3A', 'IFNG', 'GNLY', 'GZMB', 'CD69', 'HSPA1A')
CD56_int_NKcells_n <- c('NKG7', 'CD69', 'XCL1')
CD56_int_NKcells_t <- c('NKG7', 'IL7R', 'AREG', 'CCL3', 'CCL5')
tissue_resident_NKcells_tissue <- c('NKG7', 'CD160', 'TIGIT', 'CCL4', 'CCL3','IL32', 'CXCR4', 'ZNF331')
tissue_resident_NKcells_tumor <- c('NKG7', 'CD160', 'HSPA1A', 'DNAJB1', 'HSPA1B', 'HSPA6')
proliferative_NKcells <- c('NKG7', 'MKI67', 'TUBB', 'STMN')

subT_marker <- c(naiveTcells,exhaustedTcells,
                 regulatoryTcells,CD8,effector_Tcells,
                 central_memory_Tcells,tissue_resident_memory_T_cells,
                 MAIT)

T_marker_list <- split(T_marker,f = T_marker$cluster)


deg <- T_marker_list$`0`
deg <- deg[is.finite(deg$avg_log2FC),]
deg <- deg[order(deg$avg_log2FC, decreasing = T),]
genelist <- structure(deg$avg_log2FC, names = deg$gene)
res <- clusterProfiler::GSEA(genelist,TERM2GENE = gsea.human$C5,pvalueCutoff = 1,maxGSSize = 50000)
enrichplot::gseaplot2(res, geneSetID = 'HP_DELAYED_GROSS_MOTOR_DEVELOPMENT', pvalue_table = TRUE)

T_label <- c('0'='resident memory T cells',
             '1'='naive T cells',
             '2'='effector memory T cells',
             '3'='CD',
             '4'='',
             '5'='CD4 T regulatory',
             '6'='',
             '7'='CD8 exhusted')


VlnPlot(sub_T,assay = 'RNA',features = c('HAVCR2', 'PDCD1', 'GZMB', 'ITGAE', 'CXCL13'))#exhausted T cells
VlnPlot(sub_T,features = naiveTcells)#naive T cells
VlnPlot(sub_T,features = regulatoryTcells)
VlnPlot(sub_T,features = c('NKG7','NKTR','IFNG', 'GNLY', 'GZMB','CD69'))
VlnPlot(sub_T,features = c('CD160'))#
VlnPlot(sub_T,features = c('CD3D', 'CD3E', 'CD3G'))
VlnPlot(sub_T,features = c('NKTR'))#
VlnPlot(sub_T,features = c('IL7R','CD40LG'),assay = 'RNA')

#################################sub NK&T###############################################
sub_TNK <- subset(sce,ident=c('Tcell','NKcell'))
sub_TNK <- NormalizeData(sub_TNK) %>% ScaleData() %>% 
  RunPCA(features =Vargene_filter,npcs = 30) %>% 
  RunHarmony(assay.use = "RNA",group.by.vars = 'orig.ident') %>% 
  RunUMAP(reduction = "harmony", dims = dims) %>% 
  RunTSNE(reduction = "harmony", dims = dims) %>% 
  FindNeighbors(reduction = "harmony", dims = dims) %>% 
  FindClusters(resolution = resolution, verbose = FALSE)

filtergene <- c(grep('^RP[SL]',row.names(sub_TNK_)),grep('^MT-',row.names(sub_TNK_)))
NK_T_marker <- FindAllMarkers(sub_TNK_[-filtergene,],assay = "RNA")

top10 <- NK_T_marker %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
top20 <- NK_T_marker %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
top30 <- NK_T_marker %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
top50 <- NK_T_marker %>% group_by(cluster) %>% top_n(50,wt = avg_log2FC)


NKT_label <- data.frame('0'='NKT cell',
                        '1'='memory T cell',
                        '2'='NKT cell',
                        '3'='CD8A exhusted',
                        '4'='Cytotoxic cell',
                        '5'='other',
                        '6'='Treg',
                        '7'='NK cell',
                        '8'='exhusted T cell')
NKT_label <- data.frame('cluster'=c('0','1','2','3','4','5','6','7','8'),
                        'cell_types'=c('NKT cell','memory T cell','NKT cell','CD8A exhusted','Cytotoxic T cell','other','Treg','NK cell','exhusted T cell'))
NK_T <- sce_add_label(object = sub_TNK_,label = NKT_label)
DimPlot(NK_T,label = T)
table(NK_T$cell_types)

sce@meta.data$cell_subtypes <- as.character(sce$cell_types)
for (i in unique(NK_T$cell_types)) {
  cells <- row.names(NK_T@meta.data)[NK_T$cell_types==i]
  sce@meta.data[cells,]$cell_subtypes <- as.character(i)
}

table(sce$cell_subtypes)

DimPlot(sce,group.by = 'cell_subtypes',label = T)


sub_TNK_exp <- AverageExpression(sub_TNK_)$RNA
Idents(sub_TNK_) <- 'seurat_clusters'

annotation_marker <- c(Tmarker,naiveTcells,exhaustedTcells,regulatoryTcells,Cytokines,resident,TF,Thelper)
pheatmap::pheatmap(sub_TNK_exp[annotation_marker,],scale = 'row',cluster_rows = F)

Tmarker <- c('CD4','CD8A','CD8B','CD3D','CD3E','CD3G')
naiveTcells <- c('CCR7', 'LEF1', 'SELL', 'TCF7', 'S1PR1')
exhaustedTcells <- c('HAVCR2', 'PDCD1', 'GZMB', 'ITGAE','LAG3')
regulatoryTcells <- c('FOXP3', 'CTLA4', 'IL2RA')#5
Cytokines <- c('IL2','IL17A','GNLY','GZMK','GZMB','GZMA','NKG7','IFNG','LAMTOR3','PRF1')
resident <- c('CD69','RUNX3','NR4A1')
co_stimulatory <- c('TNFRSF14','CD28','ICOS','TNFRSF9')
TF <- c('TBX21','ZNF683','ZEB2','ID2','EOMES','HIF1A','TOX')
Thelper <- c('IL17A','NFKBIA','CD40LG')

sub_TNK_filter <- subset(sub_TNK_filter,downsample=1000)
DoHeatmap(object = sub_TNK_filter,features = top10$gene)

VlnPlot(sub_TNK_filter,assay = 'RNA',features = naiveTcells,pt.size = 0)

sub_TNK_filter <- subset(sub_TNK_,ident=c('0','1','2','3','4','6','7'))
sub_TNK_filter <- NormalizeData(sub_TNK_filter) %>% ScaleData() %>% 
  RunPCA(features =Vargene_filter,npcs = 30) %>% 
  RunHarmony(assay.use = "RNA",group.by.vars = 'orig.ident') %>% 
  RunUMAP(reduction = "harmony", dims = dims) %>% 
  RunTSNE(reduction = "harmony", dims = dims) %>% 
  FindNeighbors(reduction = "harmony", dims = dims) %>% 
  FindClusters(resolution = c(0.6,0.4), verbose = FALSE)

filtergene <- c(grep('^RP[SL]',row.names(sub_TNK_filter)),grep('^MT-',row.names(sub_TNK_filter)))
NK_T_marker_filter <- FindAllMarkers(sub_TNK_filter[-filtergene,],assay = "RNA")

top10 <- NK_T_marker_filter %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
top20 <- NK_T_marker_filter %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)
top30 <- NK_T_marker_filter %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
top50 <- NK_T_marker_filter %>% group_by(cluster) %>% top_n(50,wt = avg_log2FC)
DimPlot(sub_TNK_filter,label = T)
NKT_label <- data.frame('0'='NKT cell',
               '1'='memory T cell',
               '2'='CD8A exhusted',
               '3'='NKT cell',
               '4'='Cytotoxic cell',
               '5'='T helper cell',
               '6'='Treg',
               '7'='NK cell')
NKT_label <- data.frame('cluster'=c('0','1','2','3','4','5','6','7'),
                        'cell_types'=c('NKT cell','memory T cell','CD8A exhusted','NKT cell','Cytotoxic T cell','T helper cell','Treg','NK cell'))
NK_T_filter <- sce_add_label(object = sub_TNK_filter,label = NKT_label)
DimPlot(NK_T_filter,label = T)
table(NK_T_filter$cell_types)

for (i in unique(NK_T_filter$cell_types)) {
  cells <- row.names(NK_T_filter@meta.data)[NK_T_filter$cell_types==i]
  sce@meta.data[cells,]$cell_subtypes <- i
}

table(sce$cell_subtypes)

DimPlot(sce,group.by = 'cell_subtypes',label = T)

DimPlot(sub_Hepatocytes,split.by = 'group',label = TRUE)
DimPlot(sub_Hepatocytes,split.by = 'orig.ident',label = TRUE)
DimPlot(sub_Hepatocytes,group.by = 'seurat_clusters',label = TRUE)
DimPlot(sub_Hepatocytes,group.by = 'orig.ident',label = TRUE)
################################## sub Hepatocytes #############################
library(Seurat)
library(harmony)
sce <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/05AddLabel/sce_addlabel.rds')
Idents(sce) <- 'cell_types'
sub_Hepatocytes <- subset(sce,ident='Hepatocytes')

DefaultAssay(sub_Hepatocytes) <- 'RNA'
sub_Hepatocytes@assays$SCT <- NULL
sub_Hepatocytes <- SCTransform(sub_Hepatocytes)
top10 <- head(VariableFeatures(sub_Hepatocytes), 10)
#filter mt rp
Vargene <- VariableFeatures(sub_Hepatocytes)
Vargene_filter <- Vargene[-c(grep('^RP[SL]',Vargene),grep('^MT-',Vargene))] %>% as.character()
sub_Hepatocytes <- FindVariableFeatures(sub_Hepatocytes, selection.method = "vst", nfeatures = 1500)
#cluster
dims <- 1:20
resolution <- c(0.1) 
sub_Hepatocytes <- RunPCA(sub_Hepatocytes,features =Vargene_filter,npcs = 50) %>% 
  RunHarmony(assay.use = "SCT",group.by.vars = 'orig.ident') %>% 
  RunUMAP(reduction = "harmony", dims = dims) %>% 
  RunTSNE(reduction = "harmony", dims = dims) %>% 
  FindNeighbors(reduction = "harmony", dims = dims) %>% 
  FindClusters(resolution = resolution, verbose = FALSE)

marker_He <- FindAllMarkers(sub_Hepatocytes)
top10 <- marker_He %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
DefaultAssay(sub_Hepatocytes) <- 'SCT'
DoHeatmap(sub_Hepatocytes,features = top10$gene)



sub_Hepatocytes@meta.data$subtypes <- paste('Hepatocytes-C',sub_Hepatocytes@meta.data$seurat_clusters,sep = '')

for (i in unique(sub_Hepatocytes$subtypes)) {
  cells <- row.names(sub_Hepatocytes@meta.data)[sub_Hepatocytes$subtypes==i]
  sce@meta.data[cells,]$cell_subtypes <- i
}

table(sce$cell_subtypes)

DimPlot(sce,group.by = 'cell_subtypes',label = T)


DimPlot(sub_Hepatocytes,split.by = 'group',label = TRUE)
DimPlot(sub_Hepatocytes,split.by = 'orig.ident',label = TRUE)
DimPlot(sub_Hepatocytes,group.by = 'seurat_clusters',label = TRUE)
DimPlot(sub_Hepatocytes,group.by = 'orig.ident',label = TRUE)





table(sub_Hepatocytes_harmony$orig.ident,sub_Hepatocytes_harmony$seurat_clusters)

source('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/infercnv/singlecell_infercnv.R')
source('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/monocle/monocle2/monocle2.R')

Normal <- c("hcc1C","hcc3C","hcc4C")
Tumor <- c("hcc1T","hcc2T","hcc3T","hcc4T","hcc5T") 
output_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/infercnv/Hepatocytes' 
run_infercnv(object = sub_Hepatocytes,ref = 'malignant',output_dir = output_dir,Normal = Normal,Tumor = Tumor)

sub_Hepatocytes_cds <- monocle2_order(sub_Hepatocytes_harmony)
output_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/SingleCellAnnotation/06subcluster/sub_Hepatocytes'
monocle2_plot(monocle_cds = sub_Hepatocytes_cds,output_dir = output_dir,genes = unique(top10$gene))

#/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/infercnv/Hepatocytes