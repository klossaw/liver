#library
library(glue)
library(pbmcapply)
library(patchwork)
library(pdftools)
library(MetaboAnalystR)
library(jhtools)
library(imager)
library(magick)
install.packages('magick',lib = '/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
devtools::install_github('ropensci/magick',lib='/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
############################### function #######################################
subset_images <- function(sample = 'A',
                       np = 'neg',
                       mz = '100.00293',
                       fig = 'pdf',
                       merge = FALSE,
                       output_dir,
                       figname){
  output_dir %>% jhtools::checkdir()
  all_images_dir <- list.files(glue("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/images_smooth/{sample}/{np}/"),full.names = TRUE)
  print(all_images_dir)
  exists_images_dir <-  all_images_dir[grep(glue("{mz}.{fig}"),all_images_dir)]  
  if(length(exists_images_names)>0){
    if(merge==TRUE){
      pdf_combine(input = exists_images_dir,
                  output = glue("{output_dir}/{figname}merged.pdf"))
    }else{
      file.copy(exists_images_dir,output_dir,overwrite = TRUE)
    }
  }else{
    stop("no match")
  }
}


################################# transfer #####################################

select_all_sample_image <- function(mz,np = 'neg',output_dir){
  output_dir %>% jhtools::checkdir()
  samples <- c('A','A1','B','B1','C','C1','D','E')
  all_images_dir <- list.files(glue("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/images_smooth/{samples}/{np}/"),full.names = TRUE)
  exists_images_dir <-  all_images_dir[grep(glue("{mz}.jpg"),all_images_dir)] 
  print(exists_images_dir)
  
  imgA <- image_read(exists_images_dir[1]) 
  imgA1 <- image_read(exists_images_dir[2]) 
  imgB <- image_read(exists_images_dir[3]) 
  imgB1 <- image_read(exists_images_dir[4]) 
  imgC <- image_read(exists_images_dir[5]) 
  imgC1 <- image_read(exists_images_dir[6]) 
  imgD <- image_read(exists_images_dir[7]) 
  imgE <- image_read(exists_images_dir[8]) 
  
  img1 <- image_append(c(imgC1,imgC),stack = FALSE)#patient1
  img2 <- image_append(c(imgA,imgA1),stack = FALSE)#patient2
  img3 <- image_append(c(imgB1,imgB),stack = FALSE)#patient3
  img4 <- image_append(c(imgD,imgE),stack = FALSE)#patient4,5
  image_append(c(img1,img2,img3,img4),stack = TRUE) %>% image_write(path = glue("{output_dir}/{mz}.jpg"))
}


  samples <- c('A','A1','B','B1','C','C1','D','E')
  all_sample_sscc <- list.files(glue("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/1.成像图/SSCC/{samples}"),full.names = TRUE)
  
  imgA <- image_read(all_sample_sscc[grep(glue("Cluster-sscc-A-neg.jpg"),all_sample_sscc)]) 
  imgA1 <- image_read(all_sample_sscc[grep(glue("Cluster-sscc-A1-neg.jpg"),all_sample_sscc)]) 
  imgB <- image_read(all_sample_sscc[grep(glue("Cluster-sscc-B-neg.jpg"),all_sample_sscc)]) 
  imgB1 <- image_read(all_sample_sscc[grep(glue("Cluster-sscc-B1-neg.jpg"),all_sample_sscc)]) 
  imgC <- image_read(all_sample_sscc[grep(glue("Cluster-sscc-C-neg.jpg"),all_sample_sscc)]) 
  imgC1 <- image_read(all_sample_sscc[grep(glue("Cluster-sscc-C1-neg.jpg"),all_sample_sscc)]) 
  imgD <- image_read(all_sample_sscc[grep(glue("Cluster-sscc-D-neg.jpg"),all_sample_sscc)]) 
  imgE <- image_read(all_sample_sscc[grep(glue("Cluster-sscc-E-neg.jpg"),all_sample_sscc)]) 
  sscc_dir <- checkdir("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/SSCC")
  img1 <- image_append(c(imgC1,imgC),stack = FALSE)#patient1
  img2 <- image_append(c(imgA,imgA1),stack = FALSE)#patient2
  img3 <- image_append(c(imgB1,imgB),stack = FALSE)#patient3
  img4 <- image_append(c(imgD,imgE),stack = FALSE)#patient4,5
  image_append(c(img1,img2,img3,img4),stack = TRUE) %>% 
    image_write(path = glue("{sscc_dir}/allsample.jpg"))

#################################### run #######################################
  #load data
  project_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/'
  #exp
  sample_exp <- list.files('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/counts/') 
  sample_exp <- sample_exp[-grep('mean',sample_exp)]
  sample_exp_dir <- glue("{project_dir}/counts/{sample_exp}")
  #rawdata
  sample_rawdata <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/rawdata/imzML/'
  sample_rawdata_dir <- glue("{project_dir}/rawdata/imzML/{sample_rawdata}")
  #成像图
  sample_images <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/images_smooth/'
  
 
  #subset
  output_dir <- file.path(project_dir,'test1')
  samples <- c('A','A1','B','B1','C','C1','D','E')
  names(samples) <- c('mhcc4T', 'mhcc4C','mhcc3T','mhcc3C','mhcc1T','mhcc1C','mhcc2T','mhcc5T')
  subset_images(sample = samples,merge = TRUE,output_dir = output_dir)
  #results
  Qualitative <- readxl::read_xlsx('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/2.定性结果/Qualitative.xlsx')
  Qualitative_Metabolites_filter <- Qualitative[!is.na(Qualitative$Metabolites),]
  Qualitative_mz_filter <- Qualitative[!is.na(Qualitative$mz),]
  Qualitative_Formula_filter<- Qualitative[!is.na(Qualitative$Formula),]
  Qualitative_kegg_filter<- Qualitative[!is.na(Qualitative$KEGG),]
  
  
  output_dir <- file.path(project_dir,'merge_mz_smooth')
  samples <- c('A','A1','B','B1','C','C1','D','E')
  pbmclapply(Qualitative_mz_filter$mz, function(x){
    subset_images(sample = samples,
                  mz = x,
                  merge = TRUE,
                  output_dir = output_dir,
                  figname = x)
  },mc.cores = 64)
  
  output_dir <- file.path(project_dir,'mz_all_smooth')
  
  pbmclapply(Qualitative_mz_filter$mz, function(x){
    subset_images(sample = 'all',
                  mz = x,
                  merge = FALSE,
                  output_dir = output_dir,
                  figname = x)
  },mc.cores = 64)
  
  
  output_dir <- file.path(project_dir,'merge_kegg_smooth')
  samples <- c('A','A1','B','B1','C','C1','D','E')
  pbmclapply(Qualitative_kegg_filter$mz, function(x){
    figname <- Qualitative_kegg_filter %>% 
                  filter(mz==x) %>% 
                  select(KEGG) %>% 
                  as.character()%>% 
                  str_replace_all(c(' '='',';'='_','\n'=''))
    print(figname)
    subset_images(sample = samples,
                  mz = x,
                  merge = TRUE,
                  output_dir = output_dir,
                  figname = figname)
  },mc.cores = 64)
  
  output_dir <- file.path(project_dir,'merge_Metabolites_smooth')
  samples <- c('A','A1','B','B1','C','C1','D','E')
  pbmclapply(Qualitative_Metabolites_filter$mz, function(x){
    figname <- Qualitative_Metabolites_filter %>% 
      filter(mz==x) %>% 
      select(Metabolites) %>% 
      as.character()%>% 
      str_replace_all(c(' '='',';'='','\n'='','-'=''),'('='',')'='',':'='','+'='','-'='',','='')
    print(figname)
    subset_images(sample = samples,
                  mz = x,
                  merge = TRUE,
                  output_dir = output_dir,
                  figname = figname)
  },mc.cores = 64)
  
  kegg_select <- c('hsa00310',
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
           'hsa04964')
  work_dir <- c('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/')
  setwd(work_dir)
  diff <- readxl::read_xlsx('空间代谢组报告/4.比较分析/HCC_control-mean/neg/4.差异代谢物/差异代谢物.xlsx')
  output_dir <- file.path(work_dir,'diff_images')
  
  
  
  #肿瘤癌旁分析
  matrix_dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/3.选区数据/数据矩阵/"
  anno_meta <- openxlsx::read.xlsx('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/2.定性结果/Qualitative.xlsx')
  anno_meta <- anno_meta[!is.na(anno_meta$Metabolites),]
  trans_exp <- function(x){
    matrix <- openxlsx::read.xlsx(x)
    matrix <- matrix[!is.na(matrix$Metabolites),]
    anno_meta <- matrix[,1:17]
    exp <- matrix[,-(1:17)] %>% apply(2,as.numeric)
    colnames(exp) <- colnames(matrix[,-(1:17)])
    row.names(exp) <- matrix$mz
    return(exp)
  }
  
  matrix_A_neg <- trans_exp(glue("{matrix_dir}/A-ALL-neg-all.xlsx"))
  matrix_A1_neg <- trans_exp(glue("{matrix_dir}/A1-ALL-neg-all.xlsx"))
  matrix_B_neg <- trans_exp(glue("{matrix_dir}/B-ALL-neg-all.xlsx"))
  matrix_B1_neg <- trans_exp(glue("{matrix_dir}/B1-ALL-neg-all.xlsx"))
  matrix_C_neg <- trans_exp(glue("{matrix_dir}/C-ALL-neg-all.xlsx"))
  matrix_C1_neg <- trans_exp(glue("{matrix_dir}/C1-ALL-neg-all.xlsx"))
  matrix_D_neg <- trans_exp(glue("{matrix_dir}/D-ALL-neg-all.xlsx"))
  matrix_E_neg <- trans_exp(glue("{matrix_dir}/E-ALL-neg-all.xlsx"))
  
  merge_mz <- Reduce(intersect,list(row.names(matrix_A_neg),
                   row.names(matrix_A1_neg),
                   row.names(matrix_B_neg),
                   row.names(matrix_B1_neg),
                   row.names(matrix_C_neg),
                   row.names(matrix_C1_neg),
                   row.names(matrix_D_neg),
                   row.names(matrix_E_neg)))
  
  #pca
  combine_N <- cbind(matrix_C1_neg[merge_mz,],
                     matrix_A_neg[merge_mz,],
                     matrix_B1_neg[merge_mz,]) 
  combine_T <- cbind(matrix_C_neg[merge_mz,],
                     matrix_A1_neg[merge_mz,],
                     matrix_B_neg[merge_mz,],
                     matrix_D_neg[merge_mz,],
                     matrix_E_neg[merge_mz,]) 
  
  combine_all <- cbind(combine_N,combine_T)
  
  averge <- aggregate(t(combine_all),by=list(group),mean)
  averge <- aggregate(t(combine_all),by=list(group),mean)
  
  
  group <- c(rep("normal",length(colnames(combine_N))),rep("tumor",length(colnames(combine_T))))
  sample <- c(rep("C1",length(colnames(matrix_C1_neg))),
              rep("A",length(colnames(matrix_A_neg))),
              rep("B1",length(colnames(matrix_B1_neg))),
              rep("C",length(colnames(matrix_C_neg))),
              rep("A1",length(colnames(matrix_A1_neg))),
              rep("B",length(colnames(matrix_B_neg))),
              rep("D",length(colnames(matrix_D_neg))),
              rep("E",length(colnames(matrix_E_neg))))
  dim(combine_all)
  combine_all_log[1:4,1:4]
  #normalize
  combine_all_log <- log_transform(combine_all)
  combine_all_scale <- auto_scale(combine_all)
  averge_group <- aggregate(t(combine_all_scale),by=list(group),mean)
  averge_sample <- aggregate(t(combine_all_scale),by=list(sample),mean)
  averge_sample <- t(averge_sample)
  colnames(averge_sample) <- averge_sample[1,]
  averge_sample <- averge_sample[-1,]
  
  pheatmap(averge_sample)
  #pca
  pca_data <- prcomp(t(averge_sample))
  pca_x <- data.frame(pca_data$x,group = group,sample=sample)
  ggplot(pca_x,aes(x = PC1,y = PC2,color = group)) + geom_point()
  ggplot(pca_x,aes(x = PC1,y = PC2,color = sample)) + geom_point()
  #
  
  #diff
  diff_limma <- function(exp,group){
    design <- model.matrix(~group)
    colnames(design) <- levels(group)
    rownames(design) <- colnames(exp)
    fit <- lmFit(exp, design)
    fit <- eBayes(fit, trend=TRUE)
    results_limma <- topTable(fit, coef=2,n=Inf)
    return(results_limma)
  }
  results_limma <- diff_limma(combine_all_log,group)
  select_mz <- results_limma %>% filter(abs(logFC)>1&P.Value<0.05)
  select_mz_KEGG<- anno_meta[anno_meta$mz%in%row.names(select_mz),]$KEGG
  
  
  #enrichment
  library(MetaboSignal,lib.loc = '/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
  #BiocManager::install("MetaboSignal",lib = '/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
  
  pca_data <- prcomp(t(combine_all))
  head(pca_data,10)
  cluster_km <- kmeans(pca_data$x[,1:30],10)
  pca_x <- data.frame(pca_data$x,cluster = as.character( cluster_km$cluster),group = group)
  
  ggplot(pca_x,aes(x = PC1,y = PC2,color = sample)) + geom_point()
  
  
  tsne_data <- tsne::tsne(pca_data$x[,1:30])
  tsne_ <- tsne_data$layout
  colnames(tsne_) <- c('tsne1','tsne2')
  ggplot(tsne_,aes(x = tsne1,y = tsne2,color = as.character(cluster_km$cluster))) + geom_point()
  
  umap_data <- umap::umap(pca_data$x[,1:10])
  head(umap_data$layout,10)
  
  umap_ <- umap_data$layout
  colnames(umap_) <- c('umap1','umap2')
  ggplot(umap_,aes(x = umap1,y = umap2,color = as.character( cluster_km$cluster))) + geom_point()
  #
  combine_A <- cbind(matrix_A1_neg,matrix_A_neg) 
  combine_A <- log_transform(combine_A)
  group <- c(rep("normal",length(colnames(matrix_A1_neg))),rep("tumor",length(colnames(matrix_A_neg))))
 
  pca_data <- prcomp(t(combine_A))
  head(pca_data,10)
  cluster_km <- kmeans(pca_data$x[,1:30],10)
  pca_x <- data.frame(pca_data$x,cluster = as.character( cluster_km$cluster),group = group)
  
  ggplot(pca_x,aes(x = PC1,y = PC2,color = group)) + geom_point()
  
  results_limma <- diff_limma(combine_A,group)
  select_mz <- results_limma %>% filter(abs(logFC)>1&P.Value<0.05)
  select_mz_meta <- anno_meta[anno_meta$mz%in%row.names(select_mz),]$Metabolites
  
  
  #enrichment
 
  
  #plot
  
  
  #划区分析
  matrix <- readxl::read_xlsx('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/空间代谢组报告/3.选区数据/数据矩阵/E-ALL-neg-all.xlsx')
 
  matrix <- matrix[!is.na(matrix$Metabolites),]
  anno_meta <- matrix[,1:17]
  exp <- matrix[,-(1:17)] %>% apply(2,as.numeric)
  colnames(exp) <- colnames(matrix[,-(1:17)])
  row.names(exp) <- matrix$mz
  
  location <- read_table('/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/location/Esub1.txt',
                         col_names = c('number','x','y','np'))
  pre <- colnames(exp)[1] %>% str_split('-') %>% unlist() 
  pre <- pre[1:3] %>% str_c(collapse = '-')
  coord <- glue("{pre}-{location$x}-{location$y}")
  group <- rep('control',length(colnames(exp)))
  group[!is.na(match(colnames(exp),coord))] <- 'select'
  
  table(group)
  
  
  
  #T test
  design <- model.matrix(~group)
  colnames(design) <- levels(group)#normal tumor
  rownames(design) <- colnames(exp)
  fit <- lmFit(exp, design)
  fit <- eBayes(fit, trend=TRUE)
  results_limma <- topTable(fit, coef=2,n=Inf)
  select_mz_up <- results_limma %>% filter(logFC > 1.5&P.Value<0.05)
  select_mz_down <- results_limma %>% filter(logFC < (-1.5)&P.Value<0.05)
  select_mz_compund_up <- anno_meta[anno_meta$mz%in%row.names(select_mz_up),]$KEGG
  select_mz_compund_down <- anno_meta[anno_meta$mz%in%row.names(select_mz_down),]$KEGG
  
  compound_up <- select_mz_compund_up %>% str_replace_all(c(' '='','\n'='')) %>% str_split(';') %>% unlist() %>% unique()
  
  compound_down <- select_mz_compund_down %>% str_replace_all(c(' '='','\n'='')) %>% str_split(';') %>% unlist() %>% unique()
  
  enrich_compound_down <- clusterProfiler::enricher(gene = compound_down,
                                               TERM2GENE = TERM2GENE,
                                               TERM2NAME = TERM2NAME)
  
  tmp <- enrich_compound_down@result$geneID %>% str_split('/')
  
  enrich_compound_up <- clusterProfiler::enricher(gene = compound_up,
                                                    TERM2GENE = TERM2GENE,
                                                    TERM2NAME = TERM2NAME)
  enrich_group_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/enrichment/group/'
  barplot(enrich_compound_down,showCategory=20,x="GeneRatio",font.size=13) %>% 
    ggsave(filename = glue("{enrich_group_dir}/down_barplot.png"),width = 8,height = 10)
  barplot(enrich_compound_up,showCategory=20,x="GeneRatio",font.size=13) %>% 
    ggsave(filename = glue("{enrich_group_dir}/up_barplot.png"),width = 6,height = 4)
  enrich_compound_down@result %>% write_csv(file = glue("{enrich_group_dir}/down.csv"))
  enrich_compound_up@result %>% write_csv(file = glue("{enrich_group_dir}/up.csv"))
  
  diff_mz_up <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/select_mz/diff_up'
  diff_mz_down <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/select_mz/diff_down'
  tumor_associated_dir <- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/select_mz/associated'
  all_mz_dir<- '/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/select_mz/all_mz'
  select_mz_up <- results_limma %>% filter(logFC>7)
  select_mz_down <- results_limma %>% filter(logFC<(-7))
  
  tumor_associated <- c('89.02304','182.00884','87.00758','89.02304'
                       )
  lapply(row.names(select_mz_up), function(x){
    select_all_sample_image(mz = x,output_dir = diff_mz_up)
  })
  lapply(row.names(select_mz_down), function(x){
    select_all_sample_image(mz = x,output_dir = diff_mz_down)
  })
  lapply(tumor_associated, function(x){
    select_all_sample_image(mz = x,output_dir = tumor_associated_dir)
  })
  lapply(row.names(combine_all), function(x){
    select_all_sample_image(mz = x,output_dir = all_mz_dir)
  })
################################################################################
  library(ropls,lib.loc='/cluster/home/yzy_jh/sbin/R/library/R5.19/')
  
  group <- c(rep("normal",length(colnames(combine_N))),rep("tumor",length(colnames(combine_T))))
  sample <- c(rep("C1",length(colnames(matrix_C1_neg))),
              rep("A",length(colnames(matrix_A_neg))),
              rep("B1",length(colnames(matrix_B1_neg))),
              rep("C",length(colnames(matrix_C_neg))),
              rep("A1",length(colnames(matrix_A1_neg))),
              rep("B",length(colnames(matrix_B_neg))),
              rep("D",length(colnames(matrix_D_neg))),
              rep("E",length(colnames(matrix_E_neg))))
  meta_matrix <- data.frame(loc=colnames(combine_all),group=group,sample=sample)
  meta_matrix_select <- meta_matrix[sample(1:length(colnames(combine_all)),8000),]
  #PCA分析
  select_matrix <- combine_all[,meta_matrix_select$loc]
  dataMatrix <- t(select_matrix)
  pca_data <- opls(dataMatrix)
  
  #根据group分类变量标记
  group_NT <- meta_matrix_select$group %>% factor()
  group_sample <- meta_matrix_select$sample %>% factor()
  plot(pca_data,
       typeVc = "x-score",
       parAsColFcVn = group_NT)
  
  oplsda_NT<- opls(dataMatrix, group_NT)
  oplsda_sample <- opls(dataMatrix, group_sample)
  oplsda_sample
  save(oplsda_NT,oplsda_sample,file = "opls.Rdata")
  #检查训练子集的预测
  trainVi <- getSubsetVi(oplsda)
  table(group[trainVi], fitted(oplsda))
  #计算测试子集的性能
  table(group[-trainVi],
        predict(plsda_data, dataMatrix[-trainVi, ]))
  library(ggplot2)
  library(ggsci)
  library(tidyverse)
  #提取样本在 OPLS-DA 轴上的位置
  sample.score = oplsda@scoreMN %>%  #得分矩阵
    as.data.frame() %>%
    mutate(group = group,
           o1 = oplsda@orthoScoreMN[,1]) #正交矩阵
  head(sample.score)#查看
  p <- ggplot(sample.score, aes(p1, o1, color = gender)) +
    geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) + #横向虚线
    geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
    geom_point() +
    #geom_point(aes(-10,-10), color = 'white') +
    labs(x = 'P1(5.0%)',y = 'to1') +
    stat_ellipse(level = 0.95, linetype = 'solid',
                 size = 1, show.legend = FALSE) + #添加置信区间
    scale_color_manual(values = c('#008000','#FFA74F')) +
    theme_bw() +
    theme(legend.position = c(0.1,0.85),
          legend.title = element_blank(),
          legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
          axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
          axis.ticks = element_line(color = 'black'))
  tmp <- data.frame(dataMatrix)
  tmp1 <- t.test(tmp~group_NT,tmp)$p.value
  
  tmp1 <- apply(array, margin, ...)
  #6-1-展示VIP值分布情况
  ageVn <- group_NT
  pvaVn <- apply(dataMatrix, 2,
                 function(feaVn) cor.test(ageVn, feaVn)[["p.value"]])
  
  
  
  
  vipVn <- getVipVn(opls(dataMatrix, ageVn,
                         predI = 1, orthoI = NA,
                         fig.pdfC = "none"))
  quantVn <- qnorm(1 - pvaVn / 2)
  rmsQuantN <- sqrt(mean(quantVn^2))
  opar <- par(font = 2, font.axis = 2, font.lab = 2,
              las = 1,
              mar = c(5.1, 4.6, 4.1, 2.1),
              lwd = 2, pch = 16)
  plot(pvaVn, vipVn,
       col = "red",
       pch = 16,
       xlab = "p-value", ylab = "VIP", xaxs = "i", yaxs = "i")
  box(lwd = 2)
  curve(qnorm(1 - x / 2) / rmsQuantN, 0, 1, add = TRUE, col = "red", lwd = 3)
  abline(h = 1, col = "blue")
  abline(v = 0.05, col = "blue")
  par(opar)
  
  #6-2-绘制棒棒糖图VIP值筛选差异代谢物
  #VIP值帮助寻找重要的代谢物
  vip <- getVipVn(oplsda)
  vip_select <- vip[vip > 1]    #通常以VIP值>1作为筛选标准
  head(vip_select)
  
  vip_select <- cbind(oplsda[names(vip_select), ], vip_select)
  names(vip_select)[4] <- 'VIP'
  vip_select <- vip_select[order(vip_select$VIP, decreasing = TRUE), ]
  head(vip_select)    #带注释的代谢物，VIP>1 筛选后，并按 VIP 降序排序
  #对差异代谢物进行棒棒糖图可视化
  #代谢物名字太长进行转换
  vip_select$cat = paste('A',1:nrow(vip_select), sep = '')
  p2 <- ggplot(vip_select, aes(cat, VIP)) +
    geom_segment(aes(x = cat, xend = cat,
                     y = 0, yend = VIP)) +
    geom_point(shape = 21, size = 5, color = '#008000' ,fill = '#008000') +
    geom_point(aes(1,2.5), color = 'white') +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = '', y = 'VIP value') +
    theme_bw() +
    theme(legend.position = 'none',
          legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
          axis.text.x = element_text(angle = 90),
          axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
          axis.ticks = element_line(color = 'black'),
          axis.ticks.x = element_blank())
  p2
  
  