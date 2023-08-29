library(STdeconvolve)
library(SpatialExperiment)
library(jhtools)
library(pbmcapply)
################################# run stdeconvlove #############################
run_stdeconvolve <- function(cellranger_dir,
                             samples,
                             output_dir,
                             mc.core = 32){
  object <- list()
  object <- pbmclapply(samples,function(x){
  if(!dir.exists(glue("{cellranger_dir}/{x}/outs"))){
      stop("probable reason 'No such file or directory'")}
    se <- SpatialExperiment::read10xVisium(samples = glue("{cellranger_dir}/{x}/outs"),
                                           type = "sparse",
                                           data = "filtered")
    ## this is the genes x barcode sparse count matrix
    cd <- se@assays@data@listData$counts
    pos <- SpatialExperiment::spatialCoords(se)
    
    ## change column names to x and y
    ## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
    colnames(pos) <- c("y", "x")
    counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
    corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
    ldas <- fitLDA(t(as.matrix(corpus)), Ks = c(15))
    optLDA <- optimalModel(models = ldas, opt = 15)
    
    results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
  }, mc.cores = mc.core)
  return(object)
}

work_dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/"
project_dir <- file.path(work_dir,"deconvolution") %>% checkdir()
setwd(project_dir)

cellranger_dir <-('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/counts/counts/')
samples <- c('shcc1C','shcc1T','shcc2T','shcc3T','shcc4C','shcc4T','shcc5C','shcc5T')  
resluts <- run_stdeconvolve(cellranger_dir = cellranger_dir,
                            samples = samples,
                            output_dir = output_dir)

deconProp <- results$theta
  deconGexp <- results$beta
  plt <- vizAllTopics(theta = deconProp,
                      pos = pos,
                      r = 45,
                      lwd = 0,
                      showLegend = TRUE,
                      plotTitle = NA) +
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +
    
    ## outer border
    ggplot2::geom_rect(data = data.frame(pos),
                       ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                    ymin = min(y)-90, ymax = max(y)+90),
                       fill = NA, color = "black", linetype = "solid", size = 0.5) +
    
    ggplot2::theme(
      plot.background = ggplot2::element_blank()
    ) +
    
    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")
  plt
  

se <- SpatialExperiment::read10xVisium(samples = cellranger_dir,
                                       type = "sparse",
                                       data = "filtered")

## this is the genes x barcode sparse count matrix
cd <- se@assays@data@listData$counts
pos <- SpatialExperiment::spatialCoords(se)

## change column names to x and y
## for this dataset, we will visualize barcodes using "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
colnames(pos) <- c("y", "x")
counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = c(15))
optLDA <- optimalModel(models = ldas, opt = 15)

results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
plt <- vizAllTopics(theta = deconProp,
                    pos = pos,
                    r = 45,
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +
  ## outer border
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +
  
  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +
  
  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")
plt
ps <- lapply(colnames(deconProp), function(celltype) {
  
  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 1, alpha = 0.5,
           low = "white",
           high = "red") +
    
    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")
  
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)
ps <- lapply(colnames(deconProp), function(celltype) {
  
  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 1, alpha = 0.5,
           low = "white",
           high = "red") +
    
    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")
  
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)