setwd("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/cellphonedb/")
write.table(as.matrix(sce@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sce@meta.data), sce@meta.data[,'cell_type', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

setwd("/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/cellphonedb/out/")
mypvals <- read.delim("pvalues.txt", check.names = FALSE)
mymeans <- read.delim("means.txt", check.names = FALSE)



# 这些基因list很有意思啊，建议保存
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair,value = T)
chemokines <- grep("^CXC|CCL|CCR|CX3|XCL|XCR", mymeans$interacting_pair,value = T)
th1 <- grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4", 
            mymeans$interacting_pair,value = T)
th2 <- grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", 
            mymeans$interacting_pair,value = T)
th17 <- grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB", 
             mymeans$interacting_pair,value = T)
treg <- grep("IL35|IL10|FOXP3|IL2RA|TGFB", mymeans$interacting_pair,value = T)
costimulatory <- grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]", 
                      mymeans$interacting_pair,value = T)
coinhibitory <- grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR", 
                     mymeans$interacting_pair,value = T)
niche <- grep("CSF", mymeans$interacting_pair,value = T)


mymeans %>% dplyr::filter(interacting_pair %in% costimulatory)%>%
  dplyr::select("interacting_pair",starts_with("NK"),ends_with("NK"))  %>%  
  reshape2::melt() -> meansdf

colnames(meansdf)<- c("interacting_pair","CC","means")

mypvals %>% dplyr::filter(interacting_pair %in% costimulatory)%>%
  dplyr::select("interacting_pair",starts_with("NK"),ends_with("NK"))%>%  
  reshape2::melt()-> pvalsdf

colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >1))$means)

 
p  <- ggplot(pldf%>% filter(means >1),aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 25  )+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 


  ggsave(p,filename = "interacting_pair.png",width = 10,height = 10)
