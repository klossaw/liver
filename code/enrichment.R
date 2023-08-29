#BiocManager::install("KEGGREST",lib='/cluster/home/yzy_jh/sbin/R/library/R5.19/',force = T)
#BiocManager::install('pathview',lib = '/cluster/home/yzy_jh/sbin/R/library/4.1.1/')

#### get some spcific pathway information ####
library(KEGGREST,lib.loc = '/cluster/home/yzy_jh/sbin/R/library/R5.19/')
library(pathview,lib.loc = '/cluster/home/yzy_jh/sbin/R/library/4.1.1/')
library(jhtools)

metabolism_enrich <-
  c(
    '/cluster/home/danyang_jh/projects/liver/analysis/qzhang/human/metabolism/paired_test&PLS/enrichment/'
  )
# check the output directory or create it
checkdir(metabolism_enrich)
setwd(metabolism_enrich)


kegg_ref <- keggList("pathway","hsa") %>% data.frame()
kegg_ref_df <-
  data.frame(pathID = str_split(rownames(kegg_ref),
                                pattern = ':', simplify = T)[, 2],
             name = kegg_ref[, 1])

# check some specific pathways
path <- keggGet("hsa00010")
head(path)
gene.info <- path[[1]]$GENE
# extract gene symbol and Entrez ID
genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))
gene.symbol <- genes[1:length(genes)%%3 == 2]
gene.id <- genes[1:length(genes)%%3 == 1]

# get metabolites in one pathway
# extract all ID and compounds ID in this pathway
match.df <- vector()
for (i in 1:nrow(hsa_path)) {
  hsa_info <- keggGet(hsa_path[i,"pathID"])
  hsa_compound <- hsa_info[[1]]$COMPOUND
  path_name <- hsa_info[[1]]$NAME
  if(length(hsa_compound)>0)
  {
    cpd <- names(hsa_compound[1])
    cpd_name <- as.character(hsa_compound[1])
    for (j in 1:(length(hsa_compound)-1)) {
      cpd <- paste(cpd,names(hsa_compound)[j+1],sep = ";")
      cpd_name <- paste(cpd_name,as.character(hsa_compound)[j+1],sep = ";")
    }
    match.ver <- c(hsa_path[i,"pathID"],path_name,cpd,cpd_name)
    match.df <- rbind(match.df,match.ver)
  }
  if(length(hsa_compound)==0){
    match.df <- rbind(match.df,c(hsa_path[i,"pathID"],path_name,"",""))
  }
}

rownames(match.df) <- match.df[,1]
colnames(match.df) <- c("Pathway_ID","Pathway_Name","Compound_ID","Metabolism_Name")


pathview(species = 'hsa',cpd.data = compound,pathway.id = '00010')



################################ kegg ref ######################################

hsa_pathway <- KEGGREST::keggList(database = 'pathway',organism ='hsa' ) # 获取KEGG数据库中所有人类通路
hsa_path <- data.frame(hsa_pathway) # 转成数据框,方便后续分析
hsa_path$pathID <- substr(rownames(hsa_path),6,nchar(rownames(hsa_path)[1])) # 提取pathway ID
ref_list <- lapply(hsa_path$pathID, function(x){
  hsa_info <- keggGet(x)
  hsa_compound <- hsa_info[[1]]$COMPOUND
  path_name <- hsa_info[[1]]$NAME %>% str_split('-') %>% unlist()
  path_name <- path_name[1]
  ref <- data.frame("Pathway_ID" = rep(x,length(hsa_compound)),
                    "Pathway_Name" = rep(path_name,length(hsa_compound)),
                    "Compound_ID" = names(hsa_compound),
                    "Metabolism_Name" = hsa_compound)
  return(ref)
})
ref <- do.call(rbind,ref_list)
ref %>% write_tsv(
  glue(
    '/cluster/home/danyang_jh/projects/liver/analysis/qzhang/human/metabolism/paired_test&PLS/enrichment/ref_KEGG.tsv'
  )
)
TERM2GENE <- data.frame(ref[,1],ref[,3])
TERM2NAME <- data.frame(ref[,1],ref[,2])



##################################### make gsva ref ############################

neg_anno <-
  openxlsx::read.xlsx(
    '/cluster/home/yzy_jh/projects/liver/data/qzhang/human/metabolism/空间代谢组报告/2.定性结果/Qualitative.xlsx',
    sheet = 1
  )
pos_anno <-
  openxlsx::read.xlsx(
    '/cluster/home/yzy_jh/projects/liver/data/qzhang/human/metabolism/空间代谢组报告/2.定性结果/Qualitative.xlsx',
    sheet = 3
  )

for (i in 1:3) {
  diff_pairs_neg[[i]]$sel_mz$mz <-
    sprintf(as.numeric(rownames(diff_pairs_neg[[i]]$sel_mz)), fmt = '%0.5f')
  diff_pairs_pos[[i]]$sel_mz$mz <-
    sprintf(as.numeric(rownames(diff_pairs_pos[[i]]$sel_mz)), fmt = '%0.5f')
}

neg_up <- lapply(X = 1:3, function(x) {
  df <- diff_pairs_neg[[x]]$sel_mz %>% filter(logFC > 0)
  mz <- df$mz
  return(mz)
})

pos_up <- lapply(X = 1:3, function(x) {
  df <- diff_pairs_pos[[x]]$sel_mz %>% filter(logFC > 0)
  mz <- df$mz
  return(mz)
})

pos_dw <- lapply(X = 1:3, function(x) {
  df <- diff_pairs_pos[[x]]$sel_mz %>% filter(logFC < 0)
  mz <- df$mz
  return(mz)
})

neg_dw <- lapply(X = 1:3, function(x) {
  df <- diff_pairs_neg[[x]]$sel_mz %>% filter(logFC < 0)
  mz <- df$mz
  return(mz)
})

par(mfrow = c(2, 2))
gplots::venn(list(
  `A vs A1` = neg_up[[1]],
  `B vs B1` = neg_up[[2]],
  `C vs C1` = neg_up[[3]]
))
gplots::venn(list(
  `A vs A1` = neg_dw[[1]],
  `B vs B1` = neg_dw[[2]],
  `C vs C1` = neg_dw[[3]]
))
gplots::venn(list(
  `A vs A1` = pos_up[[1]],
  `B vs B1` = pos_up[[2]],
  `C vs C1` = pos_up[[3]]
))
gplots::venn(list(
  `A vs A1` = pos_dw[[1]],
  `B vs B1` = pos_dw[[2]],
  `C vs C1` = pos_dw[[3]]
))



select_mz_compund_up <-
  anno_meta[anno_meta$mz %in% row.names(select_mz_up), ]$KEGG
select_mz_compund_down <-
  anno_meta[anno_meta$mz %in% row.names(select_mz_down), ]$KEGG

compound_up <-
  select_mz_compund_up %>% str_replace_all(c(' ' = '', '\n' = '')) %>% str_split(';') %>% unlist() %>% unique()

compound_down <-
  select_mz_compund_down %>% str_replace_all(c(' ' = '', '\n' = '')) %>% str_split(';') %>% unlist() %>% unique()

enrich_compound_down <- clusterProfiler::enricher(gene = compound_down,
                                                  TERM2GENE = TERM2GENE,
                                                  TERM2NAME = TERM2NAME)

tmp <- enrich_compound_down@result$geneID %>% str_split('/')

enrich_compound_up <- clusterProfiler::enricher(gene = compound_up,
                                                TERM2GENE = TERM2GENE,
                                                TERM2NAME = TERM2NAME)



res_enrich <- lapply(1:3, function(x) {
  compound_up <-
    c(neg_anno[neg_anno$mz %in% neg_up[[x]],]$KEGG, pos_anno[pos_anno$mz %in%
                                                               pos_up[[x]],]$KEGG)
  compound_dw <-
    c(neg_anno[neg_anno$mz %in% neg_dw[[x]],]$KEGG, pos_anno[pos_anno$mz %in%
                                                               pos_dw[[x]],]$KEGG)
  compUp <-
    compound_up %>% str_replace_all(c(' ' = '', '\n' = '')) %>% str_split(';') %>% unlist() %>% unique()
  compDw <-
    compound_dw %>% str_replace_all(c(' ' = '', '\n' = '')) %>% str_split(';') %>% unlist() %>% unique()
  Up_enrich <- clusterProfiler::enricher(gene = compUp,
                                         TERM2GENE = TERM2GENE,
                                         TERM2NAME = TERM2NAME)
  Dw_enrich <- clusterProfiler::enricher(gene = compDw,
                                         TERM2GENE = TERM2GENE,
                                         TERM2NAME = TERM2NAME)
  return(list(Up_enrich, Dw_enrich))
})

up_res_df_enrich <- lapply(1:3,function(x){
  res_enrich[[x]][[1]]@result
})

dw_res_df_enrich <- lapply(1:3,function(x){
  res_enrich[[x]][[2]]@result
})

names(up_res_df_enrich) <- c('A_pair','B_pair','C_pair')
names(dw_res_df_enrich) <- names(up_res_df_enrich)


## GSEA with KEGG metabolites ####
cpd_anno_mz_neg <- openxlsx::read.xlsx(
  '/cluster/home/yzy_jh/projects/liver/data/qzhang/human/metabolism/空间代谢组报告/2.定性结果/Qualitative.xlsx',
  sheet = 2
)
cpd_anno_mz_pos <-
  openxlsx::read.xlsx(
    '/cluster/home/yzy_jh/projects/liver/data/qzhang/human/metabolism/空间代谢组报告/2.定性结果/Qualitative.xlsx',
    sheet = 4
  )
gsea_pair_res <- lapply(X = c(1:3), function(x) {
  mz_res_neg <-
    data.frame(diff_pairs_neg[[x]]$results, check.names = F)
  mz_res_pos <-
    data.frame(diff_pairs_pos[[x]]$results, check.names = F)
  
  rownames(mz_res_neg) <-
    sprintf(as.numeric(rownames(mz_res_neg)), fmt = '%0.5f')
  rownames(mz_res_pos) <-
    sprintf(as.numeric(rownames(mz_res_pos)), fmt = '%0.5f')
  
  logFC_mz_neg <- mz_res_neg$logFC
  names(logFC_mz_neg) <- rownames(mz_res_neg)
  logFC_mz_neg <- logFC_mz_neg[names(logFC_mz_neg) %in% neg_anno$mz]
  logFC_mz_pos <- mz_res_pos$logFC
  names(logFC_mz_pos) <- rownames(mz_res_pos)
  logFC_mz_pos <- logFC_mz_pos[names(logFC_mz_pos) %in% pos_anno$mz]
  
  for (i in 1:length(logFC_mz_neg)) {
    sel <- cpd_anno_mz_neg$mz %in% names(logFC_mz_neg[i])
    if (i == 1) {
      kegg_neg_logFC <- rep(logFC_mz_neg[i], sum(sel))
      names(kegg_neg_logFC) <- cpd_anno_mz_neg[sel,]$KEGG
    } else{
      temp_name <- names(kegg_neg_logFC)
      kegg_neg_logFC <-
        append(kegg_neg_logFC, rep(logFC_mz_neg[i], sum(sel)))
      names(kegg_neg_logFC) <-
        append(temp_name, cpd_anno_mz_neg[sel,]$KEGG)
    }
    kegg_neg_logFC <- kegg_neg_logFC[names(kegg_neg_logFC) != '']
  }
  
  for (i in 1:length(logFC_mz_pos)) {
    sel <- cpd_anno_mz_pos$mz %in% names(logFC_mz_pos[i])
    if (i == 1) {
      kegg_pos_logFC <- rep(logFC_mz_pos[i], sum(sel))
      names(kegg_pos_logFC) <- cpd_anno_mz_pos[sel,]$KEGG
    } else{
      temp_name <- names(kegg_pos_logFC)
      kegg_pos_logFC <-
        append(kegg_pos_logFC, rep(logFC_mz_pos[i], sum(sel)))
      names(kegg_pos_logFC) <-
        append(temp_name, cpd_anno_mz_pos[sel,]$KEGG)
    }
    kegg_pos_logFC <- kegg_pos_logFC[names(kegg_pos_logFC) != '']
  }
  
  kegg_logFC <- c(kegg_neg_logFC, kegg_pos_logFC)
  kegg_logFC <- kegg_logFC[order(kegg_logFC, decreasing = T)]
  
  gsea_mtb <-
    clusterProfiler::GSEA(geneList = kegg_logFC,
                          TERM2GENE = TERM2GENE,
                          TERM2NAME = TERM2NAME)@result
  
  return(gsea_mtb)
  
})

names(gsea_pair_res) <- c('A_pair','B_pair','C_pair')

##### plot #####

barplot(enrich_compound,showCategory=20,x="GeneRatio",font.size=13) 
dotplot(enrich_compound,showCategory=20,x="GeneRatio",font.size=13) 

metabiolism_enrich  <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/metabolism/enrichment"
metabiolism_enrich_group <-  checkdir(glue("{metabiolism_enrich}/group"))
barplot(enrich_compound,showCategory=20,x="GeneRatio",font.size=13) %>% 
  ggsave(filename = glue("{metabiolism_enrich_group}/bar_up"))

