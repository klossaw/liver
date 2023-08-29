##### load all needed packges #####
#library
# library(glue)
# library(pbmcapply)
# library(patchwork)
# library(pdftools)
# library(MetaboAnalystR)
# library(jhtools)
# library(imager)
# library(magick)
# install.packages('magick',lib = '/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
for(i in c(
  'glue',
  'pbmcapply',
  'patchwork',
  'pdftools',
  'jhtools',
  'imager',
  'magick',
  'ggsci',
  'ropls',
  'tidyverse'
)) {
  
  if (!require(i, character.only = T)) {
    install.packages(i, lib = '/cluster/home/yzy_jh/sbin/R/library/', )
  }else{
    print(paste0('Package ', "'", i, "'", ' has been installed and loaded, enjoy coding.'))
  }
  
}
# devtools::install_github('ropensci/magick', lib = '/cluster/home/yzy_jh/sbin/R/library/4.1.1/')



##### sub images #####
subset_images <- function(sample = 'A',
                          np = 'neg',
                          mz = '100.00293',
                          fig = 'pdf',
                          merge = FALSE,
                          output_dir,
                          figname) {
  output_dir %>% jhtools::checkdir()
  all_images_dir <-
    list.files(
      glue(
        "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/images_smooth/{sample}/{np}/"
      ),
      full.names = TRUE
    )
  print(all_images_dir)
  exists_images_dir <-
    all_images_dir[grep(glue("{mz}.{fig}"), all_images_dir)]
  if (length(exists_images_names) > 0) {
    if (merge == TRUE) {
      pdf_combine(input = exists_images_dir,
                  output = glue("{output_dir}/{figname}merged.pdf"))
    } else{
      file.copy(exists_images_dir, output_dir, overwrite = TRUE)
    }
  } else{
    stop("no match")
  }
}


##### combine the images selected #####

select_all_sample_image <- function(mz, np = 'neg', output_dir) {
  output_dir %>% jhtools::checkdir()
  samples <- c('A', 'A1', 'B', 'B1', 'C', 'C1', 'D', 'E')
  all_images_dir <-
    list.files(
      glue(
        "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/images_smooth/{samples}/{np}/"
      ),
      full.names = TRUE
    )
  exists_images_dir <-
    all_images_dir[grep(glue("{mz}.jpg"), all_images_dir)]
  print(exists_images_dir)
  
  imgA <- image_read(exists_images_dir[1])
  imgA1 <- image_read(exists_images_dir[2])
  imgB <- image_read(exists_images_dir[3])
  imgB1 <- image_read(exists_images_dir[4])
  imgC <- image_read(exists_images_dir[5])
  imgC1 <- image_read(exists_images_dir[6])
  imgD <- image_read(exists_images_dir[7])
  imgE <- image_read(exists_images_dir[8])
  
  img1 <- image_append(c(imgC1, imgC), stack = FALSE)#patient1
  img2 <- image_append(c(imgA, imgA1), stack = FALSE)#patient2
  img3 <- image_append(c(imgB1, imgB), stack = FALSE)#patient3
  img4 <- image_append(c(imgD, imgE), stack = FALSE)#patient4,5
  image_append(c(img1, img2, img3, img4), stack = TRUE) %>% image_write(path = glue("{output_dir}/{mz}.jpg"))
}


samples <- c('A', 'A1', 'B', 'B1', 'C', 'C1', 'D', 'E')
all_sample_sscc <-
  list.files(
    glue(
      "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/1.成像图/SSCC/{samples}"
    ),
    full.names = TRUE
  )

imgA <-
  image_read(all_sample_sscc[grep(glue("Cluster-sscc-A-neg.jpg"), all_sample_sscc)])
imgA1 <-
  image_read(all_sample_sscc[grep(glue("Cluster-sscc-A1-neg.jpg"), all_sample_sscc)])
imgB <-
  image_read(all_sample_sscc[grep(glue("Cluster-sscc-B-neg.jpg"), all_sample_sscc)])
imgB1 <-
  image_read(all_sample_sscc[grep(glue("Cluster-sscc-B1-neg.jpg"), all_sample_sscc)])
imgC <-
  image_read(all_sample_sscc[grep(glue("Cluster-sscc-C-neg.jpg"), all_sample_sscc)])
imgC1 <-
  image_read(all_sample_sscc[grep(glue("Cluster-sscc-C1-neg.jpg"), all_sample_sscc)])
imgD <-
  image_read(all_sample_sscc[grep(glue("Cluster-sscc-D-neg.jpg"), all_sample_sscc)])
imgE <-
  image_read(all_sample_sscc[grep(glue("Cluster-sscc-E-neg.jpg"), all_sample_sscc)])
sscc_dir <-
  checkdir("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/SSCC")
img1 <- image_append(c(imgC1, imgC), stack = FALSE)#patient1
img2 <- image_append(c(imgA, imgA1), stack = FALSE)#patient2
img3 <- image_append(c(imgB1, imgB), stack = FALSE)#patient3
img4 <- image_append(c(imgD, imgE), stack = FALSE)#patient4,5
image_append(c(img1, img2, img3, img4), stack = TRUE) %>%
  image_write(path = glue("{sscc_dir}/allsample.jpg"))



##### Unfinished: kegg related enrichment #####

#load data
project_dir <-
  '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/'
#exp
sample_exp <-
  list.files('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/counts/')
sample_exp <- sample_exp[-grep('mean', sample_exp)]
sample_exp_dir <- glue("{project_dir}/counts/{sample_exp}")
#rawdata
sample_rawdata <-
  '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/rawdata/imzML/'
sample_rawdata_dir <-
  glue("{project_dir}/rawdata/imzML/{sample_rawdata}")
#成像图
sample_images <-
  '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/images_smooth/'


#subset
output_dir <- file.path(project_dir, 'test1')
samples <- c('A', 'A1', 'B', 'B1', 'C', 'C1', 'D', 'E')
names(samples) <-
  c('shcc4T',
    'shcc4C',
    'shcc3T',
    'shcc3C',
    'shcc1T',
    'shcc1C',
    'shcc2T',
    'shcc5T')
subset_images(sample = samples,
              merge = TRUE,
              output_dir = output_dir)
#results
Qualitative <-
  readxl::read_xlsx(
    '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/2.定性结果/Qualitative.xlsx'
  )
Qualitative_Metabolites_filter <-
  Qualitative[!is.na(Qualitative$Metabolites), ]
Qualitative_mz_filter <- Qualitative[!is.na(Qualitative$mz), ]
Qualitative_Formula_filter <-
  Qualitative[!is.na(Qualitative$Formula), ]
Qualitative_kegg_filter <- Qualitative[!is.na(Qualitative$KEGG), ]


output_dir <- file.path(project_dir, 'merge_mz_smooth')
samples <- c('A', 'A1', 'B', 'B1', 'C', 'C1', 'D', 'E')
pbmclapply(Qualitative_mz_filter$mz, function(x) {
  subset_images(
    sample = samples,
    mz = x,
    merge = TRUE,
    output_dir = output_dir,
    figname = x
  )
}, mc.cores = 64)

output_dir <- file.path(project_dir, 'mz_all_smooth')

pbmclapply(Qualitative_mz_filter$mz, function(x) {
  subset_images(
    sample = 'all',
    mz = x,
    merge = FALSE,
    output_dir = output_dir,
    figname = x
  )
}, mc.cores = 64)


output_dir <- file.path(project_dir, 'merge_kegg_smooth')
samples <- c('A', 'A1', 'B', 'B1', 'C', 'C1', 'D', 'E')
pbmclapply(Qualitative_kegg_filter$mz, function(x) {
  figname <- Qualitative_kegg_filter %>%
    filter(mz == x) %>%
    select(KEGG) %>%
    as.character() %>%
    str_replace_all(c(' ' = '', ';' = '_', '\n' = ''))
  print(figname)
  subset_images(
    sample = samples,
    mz = x,
    merge = TRUE,
    output_dir = output_dir,
    figname = figname
  )
}, mc.cores = 64)

output_dir <- file.path(project_dir, 'merge_Metabolites_smooth')
samples <- c('A', 'A1', 'B', 'B1', 'C', 'C1', 'D', 'E')
pbmclapply(Qualitative_Metabolites_filter$mz, function(x) {
  figname <- Qualitative_Metabolites_filter %>%
    filter(mz == x) %>%
    select(Metabolites) %>%
    as.character() %>%
    str_replace_all(
      c(
        ' ' = '',
        ';' = '',
        '\n' = '',
        '-' = ''
      ),
      '(' = '',
      ')' = '',
      ':' = '',
      '+' = '',
      '-' = '',
      ',' = ''
    )
  print(figname)
  subset_images(
    sample = samples,
    mz = x,
    merge = TRUE,
    output_dir = output_dir,
    figname = figname
  )
}, mc.cores = 64)

kegg_select <- c(
  'hsa00310',
  'hsa00410',
  'hsa00592',
  'hsa00052',
  'hsa00330',
  'hsa04974',
  'hsa00970',
  'hsa00770',
  'hsa00471',
  'hsa04977',
  'hsa02010',
  'hsa00480',
  'hsa04964'
)
work_dir <-
  c('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/')
setwd(work_dir)
diff <-
  readxl::read_xlsx('空间代谢组报告/4.比较分析/HCC_control-mean/neg/4.差异代谢物/差异代谢物.xlsx')
output_dir <- file.path(work_dir, 'diff_images')






##### Cancerous vs. Para-cancerous #####

# the Prerequisites

## confirm the locations of matrix and qualitative results
## attention: the anno_meta data.frame is focused on MSI- mode, if you want to get the MSI+ mode, please add the argument 'sheet = 3'
matrix_dir <-
  "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/counts"
anno_meta <-
  openxlsx::read.xlsx(
    '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/2.定性结果/Qualitative.xlsx'
  )
anno_meta <- anno_meta[!is.na(anno_meta$Metabolites), ]

## define a function to get each sample's m/z abundance
## Please confirm your file directory and check the file format.
trans_exp <- function(x) {
  matrix <- readxl::read_xlsx(x)
  matrix <- matrix[!is.na(matrix$Metabolites), ]
  anno_meta <- matrix[, 1:17]
  exp <- matrix[, -(1:17)] %>% apply(2, as.numeric)
  colnames(exp) <- colnames(matrix[, -(1:17)])
  row.names(exp) <- matrix$mz
  return(exp)
}

# the Formal process

## get the expression of metabolites in each sample of specific MSI mode
## note: get the paired-comparison first,so that sample 'D' and 'E' are excluded

matrix_A_neg <- trans_exp(glue("{matrix_dir}/A-ALL-neg-all.xlsx"))
matrix_A1_neg <-
  trans_exp(glue("{matrix_dir}/A1-ALL-neg-all.xlsx"))
matrix_B_neg <- trans_exp(glue("{matrix_dir}/B-ALL-neg-all.xlsx"))
matrix_B1_neg <-
  trans_exp(glue("{matrix_dir}/B1-ALL-neg-all.xlsx"))
matrix_C_neg <- trans_exp(glue("{matrix_dir}/C-ALL-neg-all.xlsx"))
matrix_C1_neg <-
  trans_exp(glue("{matrix_dir}/C1-ALL-neg-all.xlsx"))
# matrix_D_neg <- trans_exp(glue("{matrix_dir}/D-ALL-neg-all.xlsx"))
# matrix_E_neg <- trans_exp(glue("{matrix_dir}/E-ALL-neg-all.xlsx"))

merge_mz <- Reduce(
  intersect,
  list(
    row.names(matrix_A_neg),
    row.names(matrix_A1_neg),
    row.names(matrix_B_neg),
    row.names(matrix_B1_neg),
    row.names(matrix_C_neg),
    row.names(matrix_C1_neg)
  )
)
#                  row.names(matrix_D_neg),
#                  row.names(matrix_E_neg)))
#
#pca
# combine_N <- cbind(matrix_C1_neg[merge_mz,],
#                    matrix_A_neg[merge_mz,],
#                    matrix_B1_neg[merge_mz,])
# combine_T <- cbind(matrix_C_neg[merge_mz,],
#                    matrix_A1_neg[merge_mz,],
#                    matrix_B_neg[merge_mz,],
#                    matrix_D_neg[merge_mz,],
#                    matrix_E_neg[merge_mz,])
#  combine_N <- cbind(matrix_A1_neg, matrix_B1_neg, matrix_C1_neg)
#  combine_T <- cbind(matrix_A_neg, matrix_B_neg, matrix_C_neg)
#                     # ,matrix_D_neg[merge_mz,], matrix_E_neg[merge_mz,])
# comb_all <- cbind(combine_N,combine_T)

# avg <- aggregate(t(comb_all),by=list(group),mean)
# avg <- aggregate(t(comb_all),by=list(group),mean)

comb_all <- cbind(
  matrix_A1_neg[merge_mz, ],
  matrix_B1_neg[merge_mz, ],
  matrix_C1_neg[merge_mz, ],
  matrix_A_neg[merge_mz, ],
  matrix_B_neg[merge_mz, ],
  matrix_C_neg[merge_mz, ]
)

gp_spot <- c(
  rep("Paracancerous", ncol(matrix_A1_neg)),
  rep("Paracancerous", ncol(matrix_B1_neg)),
  rep("Paracancerous", ncol(matrix_C1_neg)),
  rep("Cancerous", ncol(matrix_A_neg)),
  rep("Cancerous", ncol(matrix_B_neg)),
  rep("Cancerous", ncol(matrix_C_neg))
)

sp_spot <- c(
  rep("A1", ncol(matrix_A1_neg)),
  rep("B1", ncol(matrix_B1_neg)),
  rep("C1", ncol(matrix_C1_neg)),
  rep("A", ncol(matrix_A_neg)),
  rep("B", ncol(matrix_B_neg)),
  rep("C", ncol(matrix_C_neg))
)

# dim(comb_all)

## normalizing and scaling of all spots
comb_all_log <- log_transform(comb_all)
comb_all_scale <- auto_scale(comb_all)

## get the average of both sample level and group level
avg_samp <- aggregate(t(comb_all_scale), by = list(sp_spot), mean)
avg_samp <- t(avg_samp)
colnames(avg_samp) <- avg_samp[1, ]
avg_samp <- avg_samp[-1, ]

avg_group <- aggregate(t(comb_all_scale), by = list(gp_spot), mean)
avg_group <- t(avg_group)
colnames(avg_group) <- avg_group[1, ]
avg_group <- avg_group[-1, ]
# note: if you have load the 'pheatmap' package, the package name
# 'pheatmap::' could be removed
pheatmap::pheatmap(avg_samp)



##### the PCA of samples and spots #####

## PCA on the sample level

### perform PCA calculation
### attention: the matrix 'avg_samp' is composed of characters,
### you should convert it into a numeric matrix, or it doesn't work
avg_samp <- apply(avg_samp, 2, as.numeric)
avg_group <- apply(avg_group, 2, as.numeric)

rownames(avg_samp) <- rownames(comb_all_log)
rownames(avg_group) <- rownames(comb_all_log)
## attention: variables for 'group' and 'sample' arguments should match your need
gp_samp <- rep(c('Cancerous', 'Paracancerous'), 3)

pca_samp <- prcomp(t(avg_samp))
sp_samp <- rownames(pca_samp$x)

pca_samp_x <-
  data.frame(pca_samp$x, group = gp_samp, sample = sp_samp)

### drawing PCA plot with ggplot2
### basic: ggplot(pca_x,aes(x = PC1,y = PC2,color = sample)) + geom_point()
### refined: for sample level
ggplot(pca_samp_x, aes(
  x = PC1,
  y = PC2,
  color = sample,
  shape = group
)) +
  geom_point(size = 4) + geom_hline(yintercept = 0, linetype = "dashed") + # add a horizonal line
  geom_vline(xintercept = 0, linetype = "dashed") + # add a vertical line
  scale_color_igv() + # set the color palette as IGV (Integrative Genomics Viewer), dependent on 'ggsci' package
  theme_bw() +
  stat_ellipse(level = .68) + # set the confidient ellipse
  labs(
    x = paste('PC1(', round(pca_samp$sdev[1], digits = 2), '%)', sep = ''),
    y = paste('PC2(', round(pca_samp$sdev[2], digits = 2), '%)', sep = '')
  ) +
  theme(
    legend.position = c(0.88, 0.48),
    legend.background = element_blank(),
    text = element_text(
      face = 'bold',
      size = 12,
      colour = 'black'
    )
  )

### another way for PCA: ropls::opls
### BUT the package found a ONE component model in PCA mode

### for PLS-DA:
### remember the argument 'crossvalI' should smaller than sample number
pls_samp <-
  ropls::opls(
    x = t(avg_samp),
    y = gp_samp,
    crossvalI = 5,
    orthoI = 0
  )

par(mfrow = c(1, 2))
ropls::plot(
  pls_samp,
  typeVc = 'x-score',
  parAsColFcVn = gp_samp,
  parPaletteVc = ggsci::pal_d3(alpha = .8)(2)
)
ropls::plot(pls_samp, typeVc = 'x-loading')
### for OPLS-DA:
opls_samp <-
  ropls::opls(
    x = t(avg_samp),
    y = gp_samp,
    crossvalI = 5,
    orthoI = NA,
    predI = 1
  )
### all the graphs produced above are all drawn with 'graphic' system,
### which is different from 'grid' system, the base of the 'ggplot2' package
### ploting with 'ggplot2' shows as the bottom

### perform PCA at spot level, with 'prcomp' function
pca_spot <- prcomp(t(comb_all_scale))

pca_spot_x <- data.frame(pca_spot$x,
                         group = gp_spot, sample = sp_spot)

### plot PCA at spot level, same as the above
ggplot(pca_spot_x, aes(
  x = PC1,
  y = PC2,
  color = sample,
  shape = group
)) +
  geom_point(size = 4, alpha = .6) + geom_hline(yintercept = 0, linetype =
                                                  "dashed") + # add a horizonal line
  geom_vline(xintercept = 0, linetype = "dashed") + # add a vertical line
  scale_color_aaas() + # set the color palette as IGV (Integrative Genomics Viewer), dependent on 'ggsci' package
  theme_bw() +
  #stat_ellipse(mapping = aes(group=sample))+
  labs(
    x = paste('PC1(', round(pca_spot$sdev[1], digits = 2), '%)', sep = ''),
    y = paste('PC2(', round(pca_spot$sdev[2], digits = 2), '%)', sep = '')
  ) +
  theme(
    legend.position = c(0.92, 0.2),
    legend.background = element_blank(),
    text = element_text(
      face = 'bold',
      size = 12,
      colour = 'black'
    )
  )

### remove the pca object at spot level, it takes too much memory
rm(list = c('pca_spot', 'pca_spot_x'))


## another way for PCA: ropls::opls
## this function requires the samples as columns and variables as rows
## and it takes a LONG LONG TIME~~~~~~~~~~~~~~~
## model constructing:
## pca_spot <- ropls::opls(t(apply(comb_all_scale),MARGIN = 2, as.numeric))

### choice 1: set the random seed, sub-sampling by ratio of 95:5
### please choose your favorite one, but DO NOT perform the both
set.seed(1)
A1_rdm <-
  sample(c(F, T),
         ncol(matrix_A1_neg),
         prob = c(95, 5),
         replace = T)
A_rdm <-
  sample(c(F, T),
         ncol(matrix_A_neg),
         prob = c(95, 5),
         replace = T)
B1_rdm <-
  sample(c(F, T),
         ncol(matrix_B1_neg),
         prob = c(95, 5),
         replace = T)
B_rdm <-
  sample(c(F, T),
         ncol(matrix_B_neg),
         prob = c(95, 5),
         replace = T)
C1_rdm <-
  sample(c(F, T),
         ncol(matrix_C1_neg),
         prob = c(95, 5),
         replace = T)
C_rdm <-
  sample(c(F, T),
         ncol(matrix_C_neg),
         prob = c(95, 5),
         replace = T)
spot_rdm <- c(A1_rdm, A_rdm, B1_rdm, B_rdm, C1_rdm, C_rdm)

### choice 2: randomly sub-sampling, similar to sub-sampling by ratio
### please choose your favorite one, but DO NOT perform the both
set.seed(1)
spot_rdm <- sample(ncol(comb_all_scale), size = 8000)

## PCA of the selected random variables
spot_sel <-
  apply(comb_all_scale[, spot_rdm], MARGIN = 2, as.numeric)
rownames(spot_sel) <- rownames(comb_all_scale)
gp_spot_sel <- gp_spot[spot_rdm]
sp_spot_sel <- sp_spot[spot_rdm]

pca_selSpot <- ropls::opls(t(spot_sel))

## as a result, an overview with 4 sub-images return, i, ii, iii and iv for the graphs
## left to right, then up to down
## i: (top left)the scree plot, suggests components capturing the inertia
## ii: (top right) outlier, showing distances within and orthogonal to the projection plane,
##    with names of the samples with a high value for at least one of the distances are indicated
## iii: (bottom left) x-score, variance along each axis equals variance captured by each component
##      , depending on total variance of the dataMatrix X and percentage of this variance captured
##      by this component
## iv: (bottom right) x-loading: variables with most extreme values for
##     each loading are colored and labeled

## plotting:
## ropls::plot(pca_spot, typeVc = "x-score",
##            parAsColFcVn = gp_spot,
##            parLabVc = sp_spot,
##            parPaletteVc = )
ropls::plot(
  pca_selSpot,
  typeVc = "x-score",
  parAsColFcVn = gp_spot,
  parLabVc = sp_spot_sel,
  parPaletteVc = ggsci::pal_aaas(alpha = .8)(8)[2:7]
)
## PLS(-DA): with ropls::opls()
## pls_spot <- ropls::opls(t(apply(comb_all_scale),MARGIN = 2, as.numeric),
##                        y = gp_spot)
pls_selSpot <-
  ropls::opls(t(spot_sel), y = gp_spot_sel, orthoI = 0)
## plot the pls results at spot level, annotated by samples
ropls::plot(
  pls_selSpot,
  typeVc = 'x-score',
  parAsColFcVn = sp_spot_sel,
  parLabVc = sp_spot_sel,
  parEllipsesL = F,
  parPaletteVc = ggsci::pal_aaas(alpha = .8)(8)[2:7]
)
## plot the PLS results at spot level, annotated by groups
ropls::plot(
  pls_selSpot,
  typeVc = 'x-score',
  parAsColFcVn = gp_spot_sel,
  parLabVc = sp_spot_sel,
  parEllipsesL = F,
  parPaletteVc = ggsci::pal_d3(alpha = .8)(8)[1:2]
)
## another choice for drawing pls plots: extract the results and drawing with ggplot2



##### tsne and umap trail #####


## get a pca of selected spots
comb_all_scal_sel <- comb_all_scale[, spot_sel]
pca_spot <- prcomp(t(comb_all_scal_sel))
tsne_data <- tsne::tsne(pca_spot$x[, 1:30])
tsne_ <- tsne_data$layout
colnames(tsne_) <- c('tsne1', 'tsne2')
ggplot(tsne_, aes(
  x = tsne1,
  y = tsne2,
  color = as.character(cluster_km$cluster)
)) + geom_point()

umap_data <- umap::umap(pca_spot$x[, 1:10])
head(umap_data$layout, 10)

umap_ <- umap_data$layout
colnames(umap_) <- c('umap1', 'umap2')
ggplot(umap_, aes(
  x = umap1,
  y = umap2,
  color = as.character(cluster_km$cluster)
)) + geom_point()



##### differential m/z with limma #####

# differential m/z with limma
diff_limma <- function(exp,group){
  design <- model.matrix(~group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(exp)
  fit <- lmFit(exp, design)
  fit <- eBayes(fit, trend=TRUE)
  results_limma <- topTable(fit, coef=2,n=Inf)
  return(results_limma)
}
# get the differential m/z, among all spots
results_limma <- diff_limma(comb_all_log,gp_spot)
# select the enriched m/zs, and get their KEGG compound number
diff_mz <- results_limma %>% filter(abs(logFC)>1&P.Value<0.05)
diff_mz_KEGG<- anno_meta[anno_meta$mz%in%row.names(diff_mz),]$KEGG
# Q1: ONE m/z usually has AT LEAST 2 COMPOUND IDs, 
# shall we get this result for enrichment?



##### e.g. differential m/z in a pair of samples #####
pair_dif <-
  function(sample,
           mode,
           matrix_dir,
           clr_range = NULL,
           method) {
    if (mode == 'neg') {
      anno_meta <- readxl::read_xlsx(
        '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/report/2.定性结果/Qualitative.xlsx',
        sheet = 1
      )
    } else if (mode == 'pos') {
      anno_meta <- readxl::read_xlsx(
        '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/report/2.定性结果/Qualitative.xlsx',
        sheet = 3
      )
    } else{
      print('Wrong input, please check your arguments.')
      return(NULL)
    }
    anno_meta <- anno_meta[!is.na(anno_meta$Metabolites), ]
    
    samp_1 <-
      trans_exp(glue({
        matrix_dir
      }, '/', glue({
        sample
      }, '1-ALL-', {
        mode
      }, '-all.txt')))
    samp <-
      trans_exp(glue({
        matrix_dir
      }, '/', glue({
        sample
      }, '-ALL-', {
        mode
      }, '-all.txt')))
    combine <- cbind(samp_1, samp)
    combine <- log_transform(combine)
    group <-
      c(rep("Paracancerous", ncol(samp_1)), rep("Cancerous", ncol(samp)))
    sample_name <-
      c(rep(glue(sample, '1'), ncol(samp_1)), rep(sample, ncol(samp)))
    
    ### PCA and plot
    pca_spot <- prcomp(t(combine))
    cluster_km <- kmeans(pca_spot$x[, 1:30], 10)
    pca_x <-
      data.frame(
        pca_spot$x,
        cluster = as.character(cluster_km$cluster),
        sample = sample_name,
        group = group
      )
    
    p <-
      ggplot(pca_x, aes(
        x = PC1,
        y = PC2,
        color = sample,
        shape = group
      )) +
      geom_point(size = 4, alpha = .6) +
      geom_hline(yintercept = 0, linetype = "dashed") + # add a horizonal line
      geom_vline(xintercept = 0, linetype = "dashed") + # add a vertical line
      scale_color_aaas() + # set the color palette as AAAS, dependent on 'ggsci' package
      theme_bw() +
      #stat_ellipse(mapping = aes(group=sample))+
      labs(
        x = paste('PC1(', round(pca_spot$sdev[1], digits = 2), '%)', sep = ''),
        y = paste('PC2(', round(pca_spot$sdev[2], digits = 2), '%)', sep = '')
      ) +
      theme(
        legend.position = c(0.9, 0.8),
        legend.background = element_blank(),
        text = element_text(
          face = 'bold',
          size = 12,
          colour = 'black'
        )
      )
    ggsave(
      p,
      filename = glue('./', paste0(
        sample, '1 vs ', sample,  ' ', mode, '_PCA.pdf'
      )),
      width = 8,
      height = 6
    )
    rm(list = 'p')
    
    if (method == 'limma') {
      results <- diff_limma(combine, group)
      sel_mz <-
        results %>% filter(abs(logFC) > 1 & P.Value < 0.05)
    } else if (method == 'wilcox') {
      ## The wilcoxon test function DID NOT work, the bug exists
      ## Please change to limma
      cells.1 <- grep(pattern = 'Cancerous',
                      x = group,
                      value = T)
      cells.2 <- grep(pattern = 'Paracancerous',
                      x = group,
                      value = T)
      results <- find_diff(
        object = combine,
        cells.1 = cells.1,
        cells.2 = cells.2,
        sample = sample
      )
      sel_mz <- res_wilcox %>% filter(abs(logFC) > 1 & P.value < 0.05)
    }
    sel_mz_mtb <-
      anno_meta[anno_meta$mz %in% row.names(sel_mz),]$Metabolites
    
    return(list(
      results = results,
      sel_mz = sel_mz,
      sel_mz_mtb = sel_mz_mtb
    ))
    
  }

diff_pairs_neg <-
  suppressMessages(parallel::mclapply(c('A', 'B', 'C'), function(x) {
    pair_dif(
      sample = {
        x
      },
      mode = 'neg',
      matrix_dir = matrix_dir,
      method = 'limma'
    )
  }, mc.cores = 20))

diff_pairs_pos <-
  suppressMessages(parallel::mclapply(c('A', 'B', 'C'), function(x) {
    pair_dif(
      sample = {
        x
      },
      mode = 'pos',
      matrix_dir = matrix_dir,
      method = 'limma'
    )
  }, mc.cores = 20))



###### unfinished: differential m/z in sub-regions of one sample #####

#划区分析
matrix <-
  readxl::read_xlsx(
    '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/3.选区数据/数据矩阵/E-ALL-neg-all.xlsx'
  )

matrix <- matrix[!is.na(matrix$Metabolites),]
anno_meta <- matrix[, 1:17]
exp <- matrix[,-(1:17)] %>% apply(2, as.numeric)
colnames(exp) <- colnames(matrix[,-(1:17)])
row.names(exp) <- matrix$mz

location <-
  read_table(
    '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/location/Esub1.txt',
    col_names = c('number', 'x', 'y', 'np')
  )
pre <- colnames(exp)[1] %>% str_split('-') %>% unlist()
pre <- pre[1:3] %>% str_c(collapse = '-')
coord <- glue("{pre}-{location$x}-{location$y}")
group <- rep('control', ncol(exp))
group[!is.na(match(colnames(exp), coord))] <- 'select'

table(group)


# limma
design <- model.matrix(~ group)
colnames(design) <- levels(group)#normal tumor
rownames(design) <- colnames(exp)
fit <- lmFit(exp, design)
fit <- eBayes(fit, trend = TRUE)
results_limma <- topTable(fit, coef = 2, n = Inf)
select_mz_up <-
  results_limma %>% filter(logFC > 1.5 & P.Value < 0.05)
select_mz_down <-
  results_limma %>% filter(logFC < (-1.5) & P.Value < 0.05)
select_mz_compund_up <-
  anno_meta[anno_meta$mz %in% row.names(select_mz_up),]$KEGG
select_mz_compund_down <-
  anno_meta[anno_meta$mz %in% row.names(select_mz_down),]$KEGG

compound_up <-
  select_mz_compund_up %>% str_replace_all(c(' ' = '', '\n' = '')) %>% str_split(';') %>% unlist() %>% unique()

compound_down <-
  select_mz_compund_down %>% str_replace_all(c(' ' = '', '\n' = '')) %>% str_split(';') %>% unlist() %>% unique()

enrich_compound_down <-
  clusterProfiler::enricher(gene = compound_down,
                            TERM2GENE = TERM2GENE,
                            TERM2NAME = TERM2NAME)

tmp <- enrich_compound_down@result$geneID %>% str_split('/')

enrich_compound_up <-
  clusterProfiler::enricher(gene = compound_up,
                            TERM2GENE = TERM2GENE,
                            TERM2NAME = TERM2NAME)
enrich_group_dir <-
  '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/enrichment/group/'
barplot(
  enrich_compound_down,
  showCategory = 20,
  x = "GeneRatio",
  font.size = 13
) %>%
  ggsave(
    filename = glue("{enrich_group_dir}/down_barplot.png"),
    width = 8,
    height = 10
  )
barplot(
  enrich_compound_up,
  showCategory = 20,
  x = "GeneRatio",
  font.size = 13
) %>%
  ggsave(
    filename = glue("{enrich_group_dir}/up_barplot.png"),
    width = 6,
    height = 4
  )
enrich_compound_down@result %>% write_csv(file = glue("{enrich_group_dir}/down.csv"))
enrich_compound_up@result %>% write_csv(file = glue("{enrich_group_dir}/up.csv"))

diff_mz_up <-
  '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/select_mz/diff_up'
diff_mz_down <-
  '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/select_mz/diff_down'
tumor_associated_dir <-
  '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/select_mz/associated'
all_mz_dir <-
  '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/select_mz/all_mz'
select_mz_up <- results_limma %>% filter(logFC > 7)
select_mz_down <- results_limma %>% filter(logFC < (-7))

tumor_associated <-
  c('89.02304', '182.00884', '87.00758', '89.02304')
lapply(row.names(select_mz_up), function(x) {
  select_all_sample_image(mz = x, output_dir = diff_mz_up)
})
lapply(row.names(select_mz_down), function(x) {
  select_all_sample_image(mz = x, output_dir = diff_mz_down)
})
lapply(tumor_associated, function(x) {
  select_all_sample_image(mz = x, output_dir = tumor_associated_dir)
})
lapply(row.names(comb_all), function(x) {
  select_all_sample_image(mz = x, output_dir = all_mz_dir)
})



##### ropls based PCA, PLS-DA & OPLS-DA #####
library(ropls, lib.loc = '/cluster/home/yzy_jh/sbin/R/library/R5.19/')

gp_spot <-
  c(rep("Paracancerous", ncol(combine_N)), rep("tumor", ncol(combine_T)))
sp_spot <- c(
  rep("C1", ncol(matrix_C1_neg)),
  rep("A", ncol(matrix_A_neg)),
  rep("B1", ncol(matrix_B1_neg)),
  rep("C", ncol(matrix_C_neg)),
  rep("A1", ncol(matrix_A1_neg)),
  rep("B", ncol(matrix_B_neg))
)

meta_matrix <-
  data.frame(loc = colnames(comb_all),
             group = group,
             sample = sample)
meta_matrix_select <- meta_matrix[sample(1:ncol(comb_all), 8000), ]
#PCA分析
select_matrix <- comb_all[, meta_matrix_select$loc]
dataMatrix <- t(select_matrix)
pca_spot <- opls(dataMatrix)

#根据group分类变量标记
group_NT <- meta_matrix_select$group %>% factor()
group_sample <- meta_matrix_select$sample %>% factor()
plot(pca_spot,
     typeVc = "x-score",
     parAsColFcVn = group_NT)

oplsda_NT <- opls(dataMatrix, group_NT)
oplsda_sample <- opls(dataMatrix, group_sample)
oplsda_sample
save(oplsda_NT, oplsda_sample, file = "opls.Rdata")
#检查训练子集的预测
trainVi <- getSubsetVi(oplsda)
table(group[trainVi], fitted(oplsda))
#计算测试子集的性能
table(group[-trainVi],
      predict(plsda_data, dataMatrix[-trainVi,]))
library(ggplot2)
library(ggsci)
library(tidyverse)
#提取样本在 OPLS-DA 轴上的位置
sample.score = oplsda@scoreMN %>%  #得分矩阵
  as.data.frame() %>%
  mutate(group = group,
         o1 = oplsda@orthoScoreMN[, 1]) #正交矩阵
head(sample.score)#查看
p <- ggplot(sample.score, aes(p1, o1, color = gender)) +
  geom_hline(yintercept = 0,
             linetype = 'dashed',
             size = 0.5) + #横向虚线
  geom_vline(xintercept = 0,
             linetype = 'dashed',
             size = 0.5) +
  geom_point() +
  #geom_point(aes(-10,-10), color = 'white') +
  labs(x = 'P1(5.0%)', y = 'to1') +
  stat_ellipse(
    level = 0.95,
    linetype = 'solid',
    size = 1,
    show.legend = FALSE
  ) + #添加置信区间
  scale_color_manual(values = c('#008000', '#FFA74F')) +
  theme_bw() +
  theme(
    legend.position = c(0.1, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(
      color = 'black',
      size = 12,
      family = 'Arial',
      face = 'plain'
    ),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(
      color = 'black',
      size = 15,
      family = 'Arial',
      face = 'plain'
    ),
    axis.title = element_text(
      color = 'black',
      size = 15,
      family = 'Arial',
      face = 'plain'
    ),
    axis.ticks = element_line(color = 'black')
  )
tmp <- data.frame(dataMatrix)
tmp1 <- t.test(tmp ~ group_NT, tmp)$p.value

tmp1 <- apply(array, margin, ...)
#6-1-展示VIP值分布情况
ageVn <- group_NT
pvaVn <- apply(dataMatrix, 2,
               function(feaVn)
                 cor.test(ageVn, feaVn)[["p.value"]])




vipVn <- getVipVn(opls(
  dataMatrix,
  ageVn,
  predI = 1,
  orthoI = NA,
  fig.pdfC = "none"
))
quantVn <- qnorm(1 - pvaVn / 2)
rmsQuantN <- sqrt(mean(quantVn ^ 2))
opar <- par(
  font = 2,
  font.axis = 2,
  font.lab = 2,
  las = 1,
  mar = c(5.1, 4.6, 4.1, 2.1),
  lwd = 2,
  pch = 16
)
plot(
  pvaVn,
  vipVn,
  col = "red",
  pch = 16,
  xlab = "p-value",
  ylab = "VIP",
  xaxs = "i",
  yaxs = "i"
)
box(lwd = 2)
curve(
  qnorm(1 - x / 2) / rmsQuantN,
  0,
  1,
  add = TRUE,
  col = "red",
  lwd = 3
)
abline(h = 1, col = "blue")
abline(v = 0.05, col = "blue")
par(opar)

#6-2-绘制棒棒糖图VIP值筛选差异代谢物
#VIP值帮助寻找重要的代谢物
vip <- getVipVn(oplsda)
vip_select <- vip[vip > 1]    #通常以VIP值>1作为筛选标准
head(vip_select)

vip_select <- cbind(oplsda[names(vip_select),], vip_select)
names(vip_select)[4] <- 'VIP'
vip_select <-
  vip_select[order(vip_select$VIP, decreasing = TRUE),]
head(vip_select)    #带注释的代谢物，VIP>1 筛选后，并按 VIP 降序排序
#对差异代谢物进行棒棒糖图可视化
#代谢物名字太长进行转换
vip_select$cat = paste('A', 1:nrow(vip_select), sep = '')
p2 <- ggplot(vip_select, aes(cat, VIP)) +
  geom_segment(aes(
    x = cat,
    xend = cat,
    y = 0,
    yend = VIP
  )) +
  geom_point(
    shape = 21,
    size = 5,
    color = '#008000' ,
    fill = '#008000'
  ) +
  geom_point(aes(1, 2.5), color = 'white') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = '', y = 'VIP value') +
  theme_bw() +
  theme(
    legend.position = 'none',
    legend.text = element_text(
      color = 'black',
      size = 12,
      family = 'Arial',
      face = 'plain'
    ),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(
      color = 'black',
      size = 15,
      family = 'Arial',
      face = 'plain'
    ),
    axis.text.x = element_text(angle = 90),
    axis.title = element_text(
      color = 'black',
      size = 15,
      family = 'Arial',
      face = 'plain'
    ),
    axis.ticks = element_line(color = 'black'),
    axis.ticks.x = element_blank()
  )
p2

