st_sc_enrichment_score <- function(object, 
                                 method = "gsva", 
                                 ref = '/cluster/home/yzy_jh/ref/MSigDB/c2.cp.kegg.v7.5.1.symbols.gmt',
                                 output_dir = './',
                                 ncores = 16) {
  
  exp <- object@assays$Spatial@counts
  exp_matrix <- data.frame(as.matrix(exp))
  geneSets <- GSEABase::getGmt(ref)
  #AUCell
  if (method == "AUCell") {
    geneSets <- subsetGeneSets(geneSets, rownames(exp_matrix)) 
    cells_rankings <- AUCell::AUCell_buildRankings(as.matrix(exp_matrix), nCores=ncores, plotStats=F) #rank
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets, cells_rankings) #calc
    signature_exp <- AUCell::getAUC(cells_AUC) %>% data.frame()
  }
  
  #ssGSEA
  if (method == "ssgsea") {
    gsva_es <- GSVA::gsva(as.matrix(exp_matrix), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  
  #GSVA
  if (method == "gsva") {
    gsva_es <- GSVA::gsva(as.matrix(exp_matrix), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp <- data.frame(t(gsva_es))
  }
  #GSVA
  if (method == "addmodule") {
    object <- AddModuleScore(object,features = ,name = )
  }
  colnames(signature_exp) <- gsub('.','-',colnames(signature_exp))
  object@assays$enrichment <- CreateAssayObject(signature_exp)
  readr::write_rds(obejct,file = glue("{output_dir}/seurat.rds"))
  return(object)
}


umap_plot <- object@reductions$umap@cell.embeddings

row.names(umap_plot)<-colnames(object)


library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
signature_ggplot<-data.frame(umap_plot, t(signature_exp))
p <- ggplot(data=signature_ggplot, aes(x=UMAP_1, y=UMAP_2, color = signature_ggplot[,4])) +  #this plot is great
  geom_point(size = 0.3) +
  scale_fill_gradientn(colours = pal) +
  scale_color_gradientn(colours = pal) +
  #labs(color = input.pathway) +
  #xlim(0, 2)+ ylim(0, 2)+
  xlab("UMAP 1") +ylab("UMAP 2") +
  theme(aspect.ratio=1)+
  #theme_bw()
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


p
object <- AddMetaData(object = tmp, 
                      metadata = signature_exp$KEGG_N_GLYCAN_BIOSYNTHESIS,
                      col.name = 'KEGG_N_GLYCAN_BIOSYNTHESIS')
FeaturePlot(object,features = 'KEGG_N_GLYCAN_BIOSYNTHESIS')
DimPlot(object,group.by = 'orig.ident')
tmp <- st_sc_enrichment_score(st_data, 
                       method = "gsva", 
                       ref = '/cluster/home/yzy_jh/ref/MSigDB/c2.cp.kegg.v7.5.1.symbols.gmt',
                       ncores = 64)
