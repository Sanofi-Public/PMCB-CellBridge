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
  rdction <- "pca"
  if (opt$harmony) {
    vars <- strsplit(opt$harmony_var, split = ",")[[1]]
    vars <- gsub(" ", "", vars)
    # run harmony
    sobj <- RunHarmony(sobj, reduction=rdction, group.by.vars=vars)
    rdction <- "harmony"
  }
  # ===================================
  set.seed(12345)
  sobj <- FindNeighbors(sobj, reduction=rdction, 
                        dims=1:opt$seu_n_dim, k.param=opt$seu_k_param, verbose=TRUE)
  message(paste("*** Find clusters"))
  sobj <- FindClusters(sobj, resolution=opt$seu_cluster_res, verbose=TRUE)
  if (opt$tsne) {
    message(paste("*** Run TSNE"))
    sobj <- RunTSNE(sobj, reduction=rdction, dims=1:opt$seu_n_dim)  
  }
  message(paste("*** Run UMAP"))
  sobj <- RunUMAP(sobj, reduction=rdction, dims=1:opt$seu_n_dim) 
  # ===================================
  # returns
  seurat.res <- new("seuratRes",
                    sobj=sobj)
  return(seurat.res)
  # ===================================
}