#function




cell_proportion<-function(object,output_dir,group=NULL,sample=NULL){
  colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", 
              "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", 
              colorRampPalette(brewer.pal(8,'Dark2'), alpha = 1)(8))
  meta.data <- object@meta.data
  if(is.null(meta.data$sample)){
    meta.data$sample <- meta.data$orig.ident
  }
  #plot 
  if(!is.null(sample)){
    sample.color <- colors[1:length(unique(meta.data$cell_type))]
    
    p1 <- ggplot(data = meta.data,mapping = aes(x = sample, fill=cell_type))+
      geom_bar(stat = "count",width = 0.7, position = 'fill') +
      scale_fill_manual(values = sample.color) +
      labs(y = "proportion of cells", x = "") +
      theme(text = element_text(family = "serif"),
            axis.text = element_text(family = "serif", size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
            panel.border = element_blank(),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line.x = element_line (colour = "black"), axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = NA, color = NA), legend.position = "left") 
    ggsave(p1,filename = glue("{output_dir}/sample_cell_type_stat_bar.png"), width = 8, height = 6)
    
  }
  if(!is.null(group)){
    group.color <- colors[1:length(unique(meta.data$cell_type))]
    
    p2 <- ggplot(data = meta.data, mapping = aes(x = group, fill= cell_type))+
      geom_bar(stat = "count", width = 0.7, position = 'fill') +
      scale_fill_manual(values = group.color) +
      labs(y = "proportion of cells", x = "") +
      theme(text = element_text(family = "serif"),
            axis.text = element_text(family = "serif", size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
            panel.border = element_blank(),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line.x = element_line (colour = "black"), axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = NA, color = NA), legend.position = "left") 
    ggsave(p2,filename = glue("{output_dir}/group_cell_type_stat_bar.png"), width = 4, height = 6)
    
  }
}


#' @export run_seurat
run_seurat <- function(object,...){
  UseMethod(generic = 'run_seurat',object = object)
}

#' @method run_seurat
#' @export  
run_seurat.monocle <- function(object,...){
  
}  

methods(run_seurat)


  
  
  
library(Seurat)
library(SeuratData)
library(ggplot2)
library(ggpubr)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat


## extract the max value of the y axis



## main function
StackedVlnPlot<- function(object, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  modify_vlnplot<- function(object, 
                            feature, 
                            pt.size = 0, 
                            plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                            ...) {
    p<- VlnPlot(object, features = feature, pt.size = pt.size, ... )  + 
      xlab("") + ylab(feature) + ggtitle("") + 
      theme(legend.position = "none", 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.title.y = element_text(size = rel(1), angle = 0), 
            axis.text.y = element_text(size = rel(1)), 
            plot.margin = plot.margin ) #+ rotate_x_text() 
    return(p)
  }
  extract_max<- function(p){
    ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
  }
  plot_list<- purrr::map(features, function(x) modify_vlnplot(object = object,feature = x,...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  # rotate_x_text() 转置横坐标的字体，ggpubr 包，在这里添加才会只出现一次横坐标，不会出现多个横坐标
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line()) + rotate_x_text()
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# call function
features<- c("CD8A","CD14")

vlnplt <- StackedVlnPlot(object = object, features = features,pt.size = 0)
vlnplt
ggsave(filename = "./vlnplot.png", height = 18, width = 14, plot = vlnplt)


write2cellphonedb <- function(object){
  output_dir='/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/cellphonedb/'
  write.table(as.matrix(object@assays$RNA@data), glue("{output_dir}/cellphonedb_count.txt"), sep='\t', quote=F)
  meta_data <- cbind(rownames(object@meta.data), object@meta.data[,'cell_type', drop=F])  
  meta_data <- as.matrix(meta_data)
  meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
  write.table(meta_data, glue("{output_dir}/cellphonedb_meta.txt"), sep='\t', quote=F, row.names=F)
}


library(parallel)
library(foreach)
single_parallel <- function(func,iterable,...){
  #1 
  "
  :param func:被并行函数
  :param iteralbe:func的1个动态参数(vector、list)
  :param ...:func的静态参数
  :return list,与iterable等长
  "
  #2.计算计算机内核数
  cores <- parallel::detectCores(logical = FALSE)
  #3.打开并行计算
  cl <- parallel::makeCluster(cores/2)
  #4.给每个单独内核传递变量，函数等
  parallel::clusterExport(cl,deparse(substitute(func)))
  #5.开始并行计算
  result <- parallel::parLapply(cl,iterable,func,...)
  #6.关闭并行计算
  stopCluster(cl)
  return(result)
}

sys.time({
  single_parallel()
})

# 多变量并行计算
multi_parallel <- function(func,...,MoreArgs=NULL){
  "
  :param func:被并行函数
  :param ...:func的多个动态参数
  :param MoreArgs:func的静态参数(list)
  :return list
  "
  # 加载包
  library(foreach)
  library(doParallel)
  li
  # 内核数
  cores <- parallel::detectCores(logical = FALSE)
  # 打开
  cl <- parallel::makeCluster(cores/2)
  # 注册
  doParallel::registerDoParallel(cl)
  # 并行计算
  dots <- list(...)  # 动态参数list
  result <- foreach(i=seq(length(dots))) %dopar% 
    do.call(func,c(lapply(dots,`[`,i),MoreArgs))  # 数据与参数组成list传入函数
  # 关闭
  stopCluster(cl)
  return(result)
}