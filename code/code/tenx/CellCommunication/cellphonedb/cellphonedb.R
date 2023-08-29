setwd("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/cellphonedb/")
write.table(as.matrix(sce@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sce@meta.data), sce@meta.data[,'cell_type', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
