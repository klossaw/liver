suppressPackageStartupMessages(c(
  library(tcltk),
  library(glue),
  library(tidyverse))
)
################################# CARD #########################################
devtools::install_local('/cluster/home/yzy_jh/sbin/R/library/4.1.1/CARD-master.zip',
                         lib = "/cluster/home/yzy_jh/sbin/R/library/4.1.1")

# load package
library(CARD,lib.loc = "/cluster/home/yzy_jh/sbin/R/library/4.1.1")
sp_count <- getassa
sc_count <- 
sc_meta <- annotation_sce@meta.data

CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
## QC on scRNASeq dataset! ...
## QC on spatially-resolved dataset! ..
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

print(CARD_obj@Proportion_CARD[1:2,])

## set the colors. Here, I just use the colors in the manuscript, if the color is not provided, the function will use default color in the package. 
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,spatial_location = CARD_obj@spatial_location, colors = colors)
print(p1)

## select the cell type that we are interested
ct.visualize = c("Acinar_cells","Cancer_clone_A","Cancer_clone_B","Ductal_terminal_ductal_like","Ductal_CRISP3_high-centroacinar_like","Ductal_MHC_Class_II","Ductal_APOL1_high-hypoxic","Fibroblasts")
## visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 4)                                 ### number of columns in the figure panel
print(p2)


p3 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
print(p3)

CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)
## The rownames of locations are matched ...
## Make grids on new spatial locations ...

## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)
library(ggplot2)
p4 <- ggplot(location_imputation, 
             aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p4)

p5 <- CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("#440154FF","#21908CFF","#FDE725FF"),    
  NumCols = 4)                                  
print(p5)

p6 <- CARD.visualize.gene(
  spatial_expression = CARD_obj@refined_expression,
  spatial_location = location_imputation,
  gene.visualize = c("Tm4sf1","S100a4","Tff3","Apol1","Crisp3","CD248"),
  colors = NULL,
  NumCols = 6)
print(p6)

p7 <- CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = c("Tm4sf1","S100a4","Tff3","Apol1","Crisp3","CD248"),
  colors = NULL,
  NumCols = 6)
print(p7)
################################ SPOTlight #####################################
#BiocManager::install("ggcorrplot", lib = "/cluster/home/yzy_jh/sbin/R/library/4.1.1/")
#install_github("https://github.com/MarcElosua/SPOTlight",
#               lib = "/cluster/home/yzy_jh/sbin/R/library/4.1.1/")
library(ggcorrplot,lib.loc = "/cluster/home/yzy_jh/sbin/R/library/4.1.1/")
library(SPOTlight,lib.loc = "/cluster/home/yzy_jh/sbin/R/library/4.1.1/")
library(scater)
library(scran)
library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(NMF)
library(Seurat)
library(pbmcapply)

run_spotlight <- function(sc_data,
                          st_data){
  sc_data <- logNormCounts(sc_data)
  dec <- modelGeneVar(sc_data)
  hvg <- getTopHVGs(dec, n = 3000)
  colLabels(sc_data) <- colData(sc_data)$cell_types
  genes <- !grepl("^Rp[l|s]|Mt", rownames(sc_data))
  # Compute marker genes
  mgs <- scoreMarkers(sc_data)
  mgs_ls <- lapply(names(mgs), function(i){
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.6, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_ls)
  
  res <- SPOTlight(
    x = sc_data,
    y = SP_shcc1@assays$Spatial@counts,
    groups = sc_data$cell_types,
    mgs = mgs_df,
    weight_id = "mean.AUC",
    gene_id = "gene",
    group_id = "cluster"
  )
  return(res)
}
sc_data <- as.SingleCellExperiment(sce_addlabel[,sample(1:40000,5000)])



res <- pbmclapply(patient, function(x){
  st_data <- paitent_list[[x]]
  res <- SPOTlight(
    x = sc_data,
    y = st_data@assays$Spatial@counts,
    groups = sc_data$cell_types,
    mgs = mgs_df,
    weight_id = "mean.AUC",
    gene_id = "gene",
    group_id = "cluster")
  },mc.cores=16)

#run
sc_data <- SplitObject(sce_annotation,split.by = "orig.ident")
sc_list <- sc_data
st_list <- liver
names(sc_list) <- samples
names(st_list) <- samples
res_list <- pbmclapply(samples,function(x){
  sc_data <- as.SingleCellExperiment(sc_list[[x]])
  res <- run_spotlight(sc_data = sc_data,st_data = st_list[[x]])
},mc.cores=32)

names(res_list) <- samples

pbmclapply(samples,function(x){
  mat <- res_list[[x]]$mat
  imagecord <- st_list[[x]]@images$slice1@coordinates[,2:3]
  por <- mat[row.names(imagecord),]
  plt <- plotSpatialScatterpie(x = imagecord, y = por)
  ggsave(plt,filename = glue("{x}.png"))
},mc.cores = 32)

plt

# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$cell_type)
# downsample to at most 20 cells per identity
n_cells <- 20
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]


head(mat <- res$mat)[, seq_len(3)]
# Extract NMF model fit
mod <- res$NMF
plotTopicProfiles(
  x = mod,
  y = sce$cell_type,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)

plotTopicProfiles(
  x = mod,
  y = sce$cell_type,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 6)

sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)
plotCorrelationMatrix(mat)
plotInteractions(mat, which = "heatmap", metric = "prop")
plotInteractions(mat, which = "heatmap", metric = "jaccard")
plotInteractions(mat, which = "network")
# Scatterpie



ct <- colnames(mat)
mat[mat < 0.1] <- 0
# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct


 



################################## RTCD ########################################
work_dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/"
project_dir <- file.path(work_dir,"spacexr")
setwd(project_dir)
#devtools::install_local("spacexr-master.zip", 
#                         build_vignettes = FALSE,
#                         lib = "/cluster/home/yzy_jh/sbin/R/library/4.1.1")
library(spacexr,lib.loc = "/cluster/home/yzy_jh/sbin/R/library/4.1.1")
library(Matrix)
library(doParallel)
library(ggplot2)

### Load in/preprocess your data, this might vary based on your file type
if(!file.exists(file.path(savedir,'myRCTD.rds'))) {
  ### Load in/preprocess your data, this might vary based on your file type
  refdir <- system.file("extdata",'Reference/Vignette',package = 'spacexr') # directory for the reference
  counts <- read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
  rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
  meta_data <- read.csv(file.path(refdir,"meta_data.csv")) # load in meta_data (barcodes, clusters, and nUMI)
  cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list
  
  ### Create the Reference object
  reference <- Reference(counts, cell_types, nUMI)
  #> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match
  #> colSums of counts. If this is unintended, please correct this discrepancy. If
  #> this is intended, there is no problem.
  
  ## Examine reference object (optional)
  print(dim(reference@counts)) #observe Digital Gene Expression matrix
  #> [1] 384 475
  table(reference@cell_types) #number of occurences for each cell type
  datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') # directory for sample Slide-seq dataset
  counts <- read.csv(file.path(datadir,"MappedDGEForR.csv")) # load in counts matrix
  coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
  rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
  rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
  nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
  
  ### Create SpatialRNA object
  puck <- SpatialRNA(coords, counts, nUMI)
  
  ## Examine SpatialRNA object (optional)
  print(dim(puck@counts)) # observe Digital Gene Expression matrix
  hist(log(puck@nUMI,2)) # histogram of log_2 nUMI 
  
  print(head(puck@coords)) # start of coordinate data.frame
  barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
  
  # This list can be restricted if you want to crop the puck e.g. 
  # puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
  # on the plot:
  plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                       title ='plot of nUMI') 
  
  
  myRCTD <- create.RCTD(puck, reference, max_cores = 1)
  
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  
  results <- myRCTD@results
  # normalize the cell type proportions to sum to 1.
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
  
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  # Plots all weights for each cell type as in full_mode. (saved as 
  # 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  # Plots the weights for each cell type as in doublet_mode. (saved as 
  # 'results/cell_type_weights_doublets.pdf')
  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                       results$results_df) 
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
  # 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
  
  # 'results/all_doublets_type.pdf')
  plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names) 
  # a table of frequency of doublet pairs 
  doub_occur <- table(doublets$second_type, doublets$first_type) 
  plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 
}
############################## STdeconvolve ####################################
library(STdeconvolve)
library(SpatialExperiment)
library(jhtools)
library(pbmcapply)

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
#run
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
################################## destvi  ###################################
library(Seurat)
library(ggplot2)
library(reticulate)
#devtools::install_github("cellgeni/sceasy",lib='/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
library(sceasy,lib.loc ='/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
#install.packages("anndata",lib='/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
library(anndata,lib.loc='/cluster/home/yzy_jh/sbin/R/library/4.1.1/')

use_python("/cluster/home/yzy_jh/.conda/envs/scvi/bin/python3",required = TRUE)
use_condaenv("/cluster/home/yzy_jh/.conda/envs/spatial/")
py_module_available("scanpy")
py_module_available("scvi")
Sys.which("conda")
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
  
sc_list <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/deconvolution/destvi/sc_list.rds')
st_list <- read_rds('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/visium/deconvolution/destvi/st_list.rds')
cortex_sc_adata <- convertFormat(sc_list$shcc1T, 
                                 from="seurat", 
                                 to="anndata", 
                                 main_layer="counts", 
                                 drop_single_values=FALSE)
cortex_st_adata <- convertFormat(st_list$shcc1T, 
                                 from="seurat", 
                                 to="anndata", 
                                 assay="Spatial", 
                                 main_layer="counts", 
                                 drop_single_values=FALSE)
scvi$model$CondSCVI$setup_anndata(cortex_sc_adata, labels_key="cell_type")
sclvm <- scvi$model$CondSCVI(cortex_sc_adata, weight_obs=TRUE)
sclvm$train(max_epochs=as.integer(250))
# Make plot smaller.
saved <- options(repr.plot.width=6, repr.plot.height=5)

sclvm_elbo <- py_to_r(sclvm$history["elbo_train"]$astype("float64"))
p <- ggplot(data = sclvm_elbo, mapping = aes(x=as.numeric(rownames(sclvm_elbo)), y=elbo_train)) + geom_line() + xlab("Epoch") + ylab("ELBO") + xlim(10, NA)
ggsave(p,filename = 'sclvm_elbo.pnf',width = 10,height = 10)
# Revert plot settings.
options(saved)

scvi$model$DestVI$setup_anndata(cortex_st_adata)
stlvm <- scvi$model$DestVI$from_rna_model(cortex_st_adata, sclvm)
stlvm$train(max_epochs=as.integer(2500))
################################## BayesSpace #############################

BiocManager::install("BayesSpace",lib='/cluster/home/yzy_jh/sbin/R/library/R5.19/')
c('DirichletReg','RcppDist')
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
melanoma <- qTune(melanoma, qs=seq(2, 10), platform="ST", d=7)
qPlot(melanoma)
