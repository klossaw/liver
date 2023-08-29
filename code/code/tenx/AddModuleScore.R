work_dir <- "/cluster/home/yzy_jh/projects/liver/analysis/qzhang/human/tenx/"
project_dir <- file.path(work.dir,"AddModuleScore")

library(Seurat)

run_AddMduleScore <- function(object,
                              project_dir,
                              gene_list){
  project_dir %>% checkdir()
  if(is.list(gene_list)){
    lapply(gene_list, function(x){
    AddModuleScore(object = object,
                   features = x,
                   name = names(x))
    return(object)
    })  
  }else if(is.character(gene_list)){
    AddModuleScore(object = object,
                   features = gene_list,
                   name = name
                    )
    return(object)
  }
  return(object)
}


