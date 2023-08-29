
suppressMessages(library(Seurat))
suppressMessages(library(jhtools))
suppressMessages(library(copykat))

run_copkat <- function(object,
                       output_dir,
                       tumor_number,
                       ...)
output_dir %>% checkdir()
if(file.exists(glue("{output_dir}/copykat_result.rds"))){
  copykat_result <- read_rds(glue("{output_dir}/copykat_result.rds"))
}else{
  exp_rawdata <- as.matrix(object@assays$RNA@counts)
  copykat_result <- copykat(rawmat=exp_rawdata, 
                            id.type="S", 
                            ngene.chr=5, 
                            win.size=25, 
                            KS.cut=0.1, 
                            sam.name="test", 
                            distance="euclidean", 
                            norm.cell.names="",
                            output.seg="FLASE",
                            n.cores = 200,
                            ...)
  write_rds(copykat.test, file=glue("{output_dir}/copykat_result.rds"))
}
############################# copykat plot #####################################
copykat_plot <- function(copykat_result,
                         tumor_number = 2,
                         ...){
  pred.test <- data.frame(copykat_result$prediction)
  pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
  CNA.test <- data.frame(copykat_result$CNAmat)
  
  tumor.cells <- pred.test$cell.names[which(copykat_result$copykat.pred=="aneuploid")] %>% 
    gsub(":",".",.) %>% 
    substr(1,28) %>% 
    paste(".1",sep = "")
  tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
  hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
  hc.umap <- cutree(hcc,tumor_number)#define number of tumor
  
  sce <- sce[,pred.test$cell.names]
  sce@meta.data$copykat.pred <- pred.test$copykat.pred
  sce@meta.data$copykat.tumor.pred <- rep("normal", nrow(sce@meta.data))
  sce@meta.data$copykat.tumor.pred[rownames(sce@meta.data) %in% names(hc.umap[hc.umap==1])] <- "tumor cluster 1"
  sce@meta.data$copykat.tumor.pred[rownames(sce@meta.data) %in% names(hc.umap[hc.umap==2])] <- "tumor cluster 2"
  
  DimPlot(sce, label = T) %>% ggsave(filename = glue("{output_dir}/"), width = 10, height = 10)
  DimPlot(sce, 
          split.by = "orig.ident",
          group.by = "copykat.pred") %>% ggsave(filename = glue("{output_dir}/"), width = 10, height = 10)
  DimPlot(sce, 
          group.by = "copykat.tumor.pred") %>% ggsave(filename = glue("{output_dir}/"), width = 10, height = 10)
  
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
  
  chr <- as.numeric(CNA.test$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)
  
  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  com.preN <- pred.test$copykat.pred
  pred <- rbPal5(2)[as.numeric(factor(com.preN))]
  
  cells <- rbind(pred,pred)
  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
  
  heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
  
  legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
  
  
  rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
  subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
  cells <- rbind(subpop,subpop)
  
  heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
  legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')
  
  
}
################################# run ###########################################
work_dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/"
project_dir <- file.path(work_dir,"copykat")
setwd(project_dir)
copykat_result <- read_rds('data/copykat.rds')
