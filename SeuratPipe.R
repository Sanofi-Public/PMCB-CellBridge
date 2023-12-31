# ===================================
# script seurat pipe
# ===================================
seuratPipe <- function(sobj, opt) {
  sobj <- NormalizeData(sobj, 
                        normalization.method=opt$seu_nrmlz_method, 
                        scale.factor=opt$seu_scale_factor) # Normalized values are stored in pbmc[["RNA"]]@data
  sobj <- FindVariableFeatures(sobj, 
                               selection.method="vst", 
                               nfeatures=opt$seu_n_hvg)
  sobj <- ScaleData(sobj, do.scale=TRUE, do.center=TRUE)
  sobj <- RunPCA(sobj, features=VariableFeatures(object=sobj), verbose=FALSE)
  # ===================================
  sobj_noharm <- NULL
  rdction <- "pca"
  if (opt$harmony != "none") {
    sobj_noharm <- sobj
    # run harmony
    vars <- strsplit(opt$harmony, split = ",")[[1]]
    vars <- gsub(" ", "", vars)
    if (!all(vars %in% names(sobj@meta.data))) {
      msg <- paste("variable column name(s) not exist. (RE: harmony)")
      stop(msg)
    }
    sobj <- RunHarmony(sobj, reduction.use=rdction, group.by.vars=vars)
    rdction <- "harmony"
  }
  # ===================================
  set.seed(12345)
  sobj <- FindNeighbors(sobj, reduction=rdction, 
                        dims=1:opt$seu_n_dim, 
                        k.param=opt$seu_k_param, 
                        verbose=TRUE)
  message(paste("*** Find clusters"))
  sobj <- FindClusters(sobj, resolution=opt$seu_cluster_res, 
                       verbose=TRUE)
  if (opt$tsne) {
    message(paste("*** Run TSNE"))
    sobj <- RunTSNE(sobj, reduction=rdction, dims=1:opt$seu_n_dim, 
                    check_duplicates=FALSE)  
  }
  message(paste("*** Run UMAP"))
  sobj <- RunUMAP(sobj, reduction=rdction, dims=1:opt$seu_n_dim) 
  # ===================================
  if (opt$harmony != "none") {
    rdction_noharm <- "pca"
    set.seed(12345)
    sobj_noharm <- FindNeighbors(sobj_noharm, 
                                 reduction=rdction_noharm, 
                                 dims=1:opt$seu_n_dim, 
                                 k.param=opt$seu_k_param, 
                                 verbose=TRUE)
    message(paste("*** Find clusters"))
    sobj_noharm <- FindClusters(sobj_noharm, resolution=opt$seu_cluster_res, 
                         verbose=TRUE)
    if (opt$tsne) {
      message(paste("*** Run TSNE"))
      sobj_noharm <- RunTSNE(sobj_noharm, reduction=rdction_noharm, dims=1:opt$seu_n_dim,
                             check_duplicates=FALSE)  
    }
    message(paste("*** Run UMAP"))
    sobj_noharm <- RunUMAP(sobj_noharm, reduction=rdction_noharm, dims=1:opt$seu_n_dim) 
  }
  # ===================================
  if (opt$adt){
    # https://github.com/satijalab/seurat/issues/3605
    # margin=1 means "perform CLR normalization within a feature"
    # margin=2 means "perform CLR normalization within a cell"
    sobj <- NormalizeData(sobj, normalization.method="CLR", margin=2, assay="ADT")
    sobj <- ScaleData(sobj, assay="ADT")
  }
  # ===================================
  # returns
  seurat.res <- new("seuratRes",
                    sobj_noharm=sobj_noharm,
                    sobj=sobj)
  return(seurat.res)
  # ===================================
}