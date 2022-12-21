# Classes, Generics and Methods

#### Generics ####

setClassUnion("NumNULL", members=c("numeric", "NULL"))
setClassUnion("ListNULL", members=c("list", "NULL"))
setClassUnion("DfNULL", members=c("data.frame", "NULL"))

setClass("readinRes", 
         slots=c(obj_ls="list",
                 meta_data="data.frame",
                 meta_data_ext="DfNULL",
                 obj_summ="data.frame"))

setClass("scrubletRes", 
         slots=c(obj_scrub_ls="ListNULL"))

setClass("qcfilterRes", 
         slots=c(sobj_flt_ls="list",
                 flt_summ="data.frame"))

setClass("mergeRes", 
         slots=c(sobj="Seurat"))

setClass("seuratRes", 
         slots=c(sobj="Seurat"))

setClass("fatlasRes", 
         slots=c(sobj="Seurat"))

setClass("signacxRes", 
         slots=c(sobj="Seurat"))

setClass("sargentRes", 
         slots=c(sobj="Seurat"))

setClass("signaturesRes", 
         slots=c(sobj="Seurat"))

#### Specific ####

`%nin%` = Negate(`%in%`)