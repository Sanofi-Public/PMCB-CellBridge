# ===================================
# Run Scrublet on each sample separately
# ===================================
forceAtlas <- function(sobj, opt) {
  # ===================================
  E <- GetAssayData(sobj, assay="RNA", slot="counts")
  E <- as(Matrix::t(E), "TsparseMatrix")
  gene_list <- rownames(sobj)
  # ===================================
  if (!is.null(opt$harmony)){
    rdction <- Embeddings(sobj, reduction="harmony")
    res <- make_spring_plot(E=E, 
                            gene_list=gene_list, 
                            precomputed_pca=rdction,
                            num_pc=opt$spr_n_dim)
  } else {
    res <- make_spring_plot(E=E, 
                            gene_list=gene_list,
                            num_pc=opt$spr_n_dim) 
  }
  # ===================================
  names(res) <- c("positions", "links", "info_dict") 
  # ===================================
  spring_coords <- res$positions
  stopifnot(nrow(spring_coords) == dim(sobj)[2])
  # head(spring_coords)
  colnames(spring_coords) <- paste0("spring_", 1:2)
  rownames(spring_coords) <- colnames(sobj)
  # head(spring_coords)
  sobj[["spring"]] <- CreateDimReducObject(embeddings=Matrix::as.matrix(spring_coords), 
                                           key="SPRING_", assay="RNA")
  # ===================================
  edges <- lapply(1:length(res$links), function(x){
    unlist(res$links[[x]])
  })
  edges <- data.frame(matrix(unlist(edges), nrow=length(edges), byrow=TRUE))
  sobj@misc$edges <- edges
  # ===================================
  sobj@misc$info_dict <- res$info_dict
  # ===================================
  # returns
  fatlas.res <- new("fatlasRes",
                    sobj=sobj)
  return(fatlas.res)
  # ===================================
}
