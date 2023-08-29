suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(enrichR))
suppressMessages(library(future))
suppressMessages(library(jhtools))
suppressMessages(library(glue))
library(org.Hs.eg.db)
############################ FUNCTION ##########################################
enrich_analysis <- function(genes,
                            info,
                            ...
                            ){
  TERM2GENE <- data.frame(info[,1], info[,3])
  TERM2NAME <- data.frame(info[,1], info[,2])
  TERM2NAME <- unique(TERM2NAME[order(TERM2NAME[,1]),])
  enrich <- clusterProfiler::enricher(genes,
                                      TERM2GENE = TERM2GENE,
                                      TERM2NAME = TERM2NAME,
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1,
                                      pAdjustMethod = "fdr")
  return(enrich)
}

###################################### GO/KEGG ######################################
run_gokegg <- function(genes,GO = TRUE,KEGG = TRUE,
                       output_dir, width = 10, 
                       height = 10,name = 'sce'){
  output_dir %>% checkdir()
  if(GO == TRUE){
    enrich_GO <- enrich_analysis(genes = genes,
                   info = GO.info
                  )
    if(!is.null(enrich_GO)){
      df2excel(enrich_GO@result, glue("{output_dir}/{name}GO_result.xlsx"))
      if(length(enrich_GO@result$ID) > 20){
        p.GO <- barplot(enrich_GO,showCategory = 20,title = "Enrichment GO")
        ggsave(p.GO, filename = glue("{output_dir}/{name}GO.png"), width = width, height = height)
      }else{
        p.GO <- barplot(enrich_GO,title = "Enrichment GO")
        ggsave(p.GO, filename = glue("{output_dir}/{name}GO.png"), width = width, height = height)
      }
    }
  }
  
  if(KEGG == TRUE){
    enrich_KEGG <- enrich_analysis(genes = genes,
                    info = KEGG.info
                    )
    if(!is.null(enrich_KEGG)){
      df2excel(enrich_KEGG@result, glue("{output_dir}/{name}KEGG_result.xlsx"))
      if(length(enrich_KEGG@result$ID) > 20){
        p.KEGG <- barplot(enrich_KEGG,showCategory = 20,title = "Enrichment KEGG")
        ggsave(p.KEGG, filename = glue("{output_dir}/{name}KEGG.png"), width = width, height = height)
      }else{
        p.KEGG <- barplot(enrich_GO,title = "Enrichment KEGG")
        ggsave(p.KEGG, filename = glue("{output_dir}/{name}KEGG.png"), width = width, height = height)
      }
    } 
  }
}

##################################### GSEA #####################################
run_gsea <- function(geneset,
                     deg,
                     output_dir,
                     width = 10, 
                     height = 10,
                    ...){
  output_dir %>% checkdir()
  deg <- deg[is.finite(deg$avg_log2FC),]
  deg <- deg[order(deg$avg_log2FC, decreasing = T),]
  genelist <- structure(deg$avg_log2FC, names = deg$gene)
  res <- clusterProfiler::GSEA(genelist,TERM2GENE = geneset,pvalueCutoff = 1,maxGSSize = 50000)
  df2excel(res@result, fn.excel = glue("{output_dir}/gsea_result.xlsx"),showRow = TRUE)
  if(length(row.names(res@result)) > 1000 ){
    resultID <- res@result %>% 
                dplyr::filter(pvalue < 0.01,abs(NES) > 1) %>% 
                dplyr::select(ID) 
  }else if(length(row.names(res@result)) > 100){
    resultID <- res@result %>% 
                dplyr::filter(pvalue < 0.05,abs(NES) > 1) %>% 
                dplyr::select(ID) 
  }else{resultID <- res@result %>% select(ID)}
  for(j in resultID$ID){
      print(j)
      p.gsea <- enrichplot::gseaplot2(res, geneSetID = j, title = j,pvalue_table = TRUE)
      ggsave(p.gsea, filename = glue("{output_dir}/{j}.png"), width = width, height = height)
  }
  return(res)
 }
##################################### GSVA #####################################
run_gsva <- function(geneset,
                    exp,
                    method = "ssgsea",
                    output_dir,
                    width = 10, 
                    height = 10,
                    name = 'C7'){
  output_dir %>% checkdir()
  #methods=c("gsva", "ssgsea", "zscore", "plage")
  message(paste("run",method,sep = "---"))
  gsva.res <- GSVA::gsva(expr = exp, geneset, method = method)
  gsva.res %>% data.frame() %>% 
  df2excel(glue("{output_dir}/result.xlsx"),showRow = T)
  gsva.res <- gsva.res[order(apply(gsva.res,1,var),decreasing = TRUE),]
  p.gsva<-pheatmap::pheatmap(gsva.res[1:100,],scale = 'row',cluster_rows = TRUE,
                             fontsize_row = 7,cluster_cols = FALSE)
  ggsave(p.gsva, filename = glue("{output_dir}/{name}gsva.png"), width = width, height = height)
  return(gsva.res)
} 

################################## run_enrichment ##############################
run_singlecell_enrichment <- function(object,
                                      project_dir,
                                      split_object = "seurat_clusters",
                                      ...){
  project_dir %>% checkdir()
  if(file.exists(glue("{project_dir}/marker.xlsx"))){
    df <- readxl::read_xlsx(glue("{project_dir}/marker.xlsx"))
  }else{
    options(future.globals.maxSize = 20 * 1024^3)
    plan(multisession, workers = 20) 
    Idents(object) <- split_object
    df <- FindAllMarkers(object,assay = "RNA", logfc.threshold = 0, min.pct = 0)
    plan("sequential")
    df2excel(df, fn.excel =  glue("{project_dir}/marker.xlsx"))
  }
  if(file.exists(glue("{project_dir}/exp.xlsx"))){
    exp <- readxl::read_xlsx(glue("{project_dir}/exp.xlsx")) %>% as.matrix()
  }else{
    exp <- AverageExpression(object,assays = "RNA")$RNA
    exp %>% data.frame() %>% 
    df2excel(fn.excel =  glue("{project_dir}/exp.xlsx"))
  }

  ##GO/KEGG
  df_list <- split(df, f = df$cluster)
  df_dir <- file.path(project_dir,"Go_KEGG") %>% checkdir()
  select_list <- df %>% group_by(cluster) %>% filter(avg_log2FC > 0.5, p_val < 0.05) %>% split(f = .$cluster)
  
  for (i in names(select_list)) {
    
    output_dir <- file.path(df_dir,paste("cluster",i,sep = "")) %>% checkdir()
    print(output_dir)
    run_gokegg(genes = select_list[[i]]$gene,
               output_dir = output_dir
    )
  } 
  ##GSEA
  gsea_dir <- file.path(project_dir,"GSEA") %>% checkdir()
  if(species=="human"){genesets <- gsea.human}else{genesets <- gsea.mouse}
  
  for (i in names(df_list)){
    output1_dir <- file.path(gsea_dir,paste("cluster",i,sep = "")) %>% checkdir()
    print(output1_dir)
    for (j in category) {
      output2_dir <- file.path(output1_dir ,paste("Geneset",j,sep = "")) %>% checkdir()
      print(output2_dir)
      run_gsea(geneset = genesets[[j]],
               deg = df_list[[i]],
               output_dir = output2_dir)
    }
  }
  ##GSVA
  gsva_dir <- file.path(project_dir,"GSVA") %>% checkdir()
  if(species=="human"){
    genesets <- gsva.human
  }else{
    genesets <- gsva.mouse
  }
  
  for (i in category) {
    output_dir <- file.path(gsva_dir ,paste("Geneset",i,sep = "")) %>% checkdir()
    run_gsva(geneset = genesets[[i]],
             exp = exp,
             method = "gsva",
             output_dir = output_dir)
  }
}
run_spatial_enrichment <- function(object,
                                   project_dir,
                                   split_object = "seurat_clusters"){
  project_dir %>% checkdir()
  if(file.exists(glue("{project_dir}/marker.xlsx"))){
    df <- readxl::read_xlsx(glue("{project_dir}/marker.xlsx"))
  }else{
    options(future.globals.maxSize = 20 * 1024^3)
    plan(multisession, workers = 20) 
    Idents(object) <- split_object
    df <- FindAllMarkers(object, logfc.threshold = 0, min.pct = 0)
    plan("sequential")
    df2excel(df, fn.excel =  glue("{project_dir}/marker.xlsx"))
  }
  if(file.exists(glue("{project_dir}/exp.xlsx"))){
    exp <- readxl::read_xlsx(glue("{project_dir}/exp.xlsx")) %>% as.matrix()
  }else{
    exp <- AverageExpression(object,assays = "RNA")$RNA
    exp %>% data.frame() %>% 
      df2excel(fn.excel =  glue("{project_dir}/exp.xlsx"))
  }
  
  ##GO/KEGG
  df_list <- split(df, f = df$cluster)
  df_dir <- file.path(project_dir,"Go_KEGG") %>% checkdir()
  select_list <- df %>% group_by(cluster) %>% filter(avg_log2FC > 0.5, p_val < 0.05) %>% split(f = .$cluster)
  
  for (i in names(select_list)) {
    output_dir <- file.path(df_dir,paste("cluster",i,sep = "")) %>% checkdir()
    print(output_dir)
    run_gokegg(genes = select_list[[i]]$gene,
               output_dir = output_dir
    )
  }
  ##GSEA
  gsea_dir <- file.path(project_dir,"GSEA") %>% checkdir()
  if(species=="human"){genesets <- gsea.human}else{genesets <- gsea.mouse}
  
  for (i in names(df_list)){
    output1_dir <- file.path(gsea_dir,paste("cluster",i,sep = "")) %>% checkdir()
    print(output1_dir)
    for (j in category) {
      output2_dir <- file.path(output1_dir ,paste("Geneset",j,sep = "")) %>% checkdir()
      print(output2_dir)
      run_gsea(geneset = genesets[[j]],
               deg = df_list[[i]],
               output_dir = output2_dir)
    }
  }
  ##GSVA
  gsva_dir <- file.path(project_dir,"GSVA") %>% checkdir()
  if(species=="human"){
    genesets <- gsva.human
  }else{
    genesets <- gsva.mouse
  }
  
  for (i in category) {
    output_dir <- file.path(gsva_dir ,paste("Geneset",i,sep = "")) %>% checkdir()
    run_gsva(geneset = genesets[[i]],
             exp = exp,
             method = "gsva",
             output_dir = output_dir)
  }
}
################################################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr) #提供MSigdb数据库基因集
library(fgsea)
library(dplyr)
library(tibble)
library(Seurat)

data(gcSample)
DefaultAssay(SP_shcc5) <- 'Spatial'
sp5 <- NormalizeData(SP_shcc5,assay = 'Spatial') %>% FindAllMarkers(assay = 'Spatial')

top30 <- sp5 %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
gcSample <- split(sp5,sp5$cluster)
cluster_genes <- gcSample[['0']]
gs <-bitr(row.names(cluster_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
markers1<-cbind(cluster_genes[gs[,1],],gs)
geneList = markers1$avg_log2FC
names(geneList) = markers1$ENTREZID
geneList = sort(geneList,decreasing = T)

ego2 <- enrichGO(gene         = de,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "all",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

kk <- enrichKEGG(gene         = de,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

head(ego2, 3)                

dotplot(ego2)

dotplot(kk)

de <- names(geneList)[abs(geneList) > 0.6]
de
edo <- DOSE::enrichDGN(de)
## convert gene ID to Symbol
edox <- clusterProfiler::setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

edox2 <- enrichplot::pairwise_termsim(edox)
p1 <- enrichplot::treeplot(edox2)
ggsave(p1,filename = 'tmp.png',width = 10,height = 10)

p2 <- enrichplot::treeplot(edox2, hclust_method = "average")
p1


aplot::plot_list(p1, p2, tag_levels='A')
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)                     
p1 <- emapplot(xx)
p2 <- emapplot(xx, legend_n=2) 
p3 <- emapplot(xx, pie="count")
p4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])



cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') 
     tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=",") %>% 
     dplyr::select(cellMarker, geneID) %>%
     dplyr::mutate(geneID = strsplit(geneID, ','))
cell_markers
sce.markers <- T_marker

ids=bitr(sp5$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(sp5,ids,by.x='gene',by.y='SYMBOL')
marker_list <- split(sp5$gene,sp5$cluster)
xx <- compareCluster(marker_list, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
dotplot(xx) + theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))

y <- compareCluster(marker_list, fun='enrichr',TERM2GENE = TERM2NAME,minGSSize=1,pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1)

dotplot(y, showCategory=10,includeAll=TRUE)

xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)                     
p1 <- emapplot(xx)
p2 <- emapplot(xx, legend_n=2) 
p3 <- emapplot(xx, pie="count")
p4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


####################################### main ###################################
if(sys.nframe()==0){
  
  BP <- GSEABase::getGmt('/cluster/home/yzy_jh/ref/MSigDB/C5GO/c5.go.bp.v7.5.1.symbols.gmt')
  GO.info <- read_rds("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/enrichment/extdata/GO.info.rds")
  KEGG.info <- read_rds("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/enrichment/extdata/KEGG.info.rds")
  gsea.human <- read_rds("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/enrichment/extdata/gsea.human.rds")
  gsva.human <- read_rds("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/enrichment/extdata/gsva.human.rds")
  category <- c('C2','C5','C7','C8','H')
  species <- "human"
  project_dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium//enrichment_test/"
  run_enrichment(object = sce,project_dir = project_dir)
}

