# ===================================
# Run sargent annotation 
# ===================================
Trajectory <- function(sobj, opt) {
  # ===================================
  sds_ls <- list()
  if (opt$trajectory != "none"){
    if (opt$trajectory == "umap" | opt$trajectory == "both") {
      # slingshot
      cat(paste("*** running slingshot on umap..."), sep = "\n")
      sds <- slingshot(data = Embeddings(sobj, "umap"), 
                       reducedDim = "umap",
                       clusterLabels = sobj$seurat_clusters)    
      sds_ls[[length(sds_ls) + 1]] <- sds
      names(sds_ls)[length(sds_ls)] <- paste("seurat_clusters", "umap", sep="-")
    }
    if ((opt$trajectory == "tsne" | opt$trajectory == "both") & opt$tsne) {
      # slingshot
      cat(paste("*** running slingshot on tsne ..."), sep = "\n")
      sds <- slingshot(data = Embeddings(sobj, "tsne"), 
                       reducedDim = "tsne",
                       clusterLabels = sobj$seurat_clusters)    
      sds_ls[[length(sds_ls) + 1]] <- sds
      names(sds_ls)[length(sds_ls)] <- paste("seurat_clusters", "tsne", sep="-")
    }
  } else {
    sds_ls <- NULL
  }
  # ===================================
  # returns
  traj.res <- new("trajectoryRes",
                  sds_ls=sds_ls)
  return(traj.res)
  # ===================================
}