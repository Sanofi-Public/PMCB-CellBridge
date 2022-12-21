# ===================================
# Run Scrublet on each sample separately
# ===================================
scrubletR <- function(obj_ls, opt) {
  # ===================================
  if (is.null(opt$scr_th)) {
    scrub.res <- new("scrubletRes",
                     obj_scrub_ls=NULL)
    return(scrub.res)
  }
  # ===================================
  obj_scrub_ls <- list()
  for (i in 1:length(obj_ls)) {
    # ===================================
    obj <- obj_ls[[i]]
    smpl <- names(obj_ls)[i]
    # ===================================
    message(paste("***", smpl, "doublet scores"))
    # ===================================
    # Load the raw counts matrix as a sparse matrix with cells as rows and genes as columns.
    X <- as(Matrix::t(obj), "TsparseMatrix")
    i <- as.integer(X@i)
    j <- as.integer(X@j)
    val <- X@x
    dim <- as.integer(X@Dim)
    # ===================================
    scrublet_py_args <- c(list(i=i, j=j, val=val, dim=dim#,
                               # expected_doublet_rate=opt$scr_expected_doublet_rate,
                               # min_counts=opt$scr_min_counts, 
                               # min_cells=opt$scr_min_cells, 
                               # min_gene_variability_pctl=opt$scr_min_gene_variability_pctl,
                               # n_prin_comps=opt$scr_n_prin_comps,
                               # sim_doublet_ratio=opt$scr_sim_doublet_ratio,
                               # n_neighbors=if (is.null(opt$scr_n_neighbors)) { round(0.5*sqrt(nrow(X))) }
                               ))
    scrub_res <- do.call(scrublet_py, scrublet_py_args)
    names(scrub_res) <- c("doublet_scores_obs", "doublet_scores_sim") # "predicted_doublets", 
    # ===================================
    stopifnot(length(scrub_res$doublet_scores_obs) == dim(obj)[2])
    # stopifnot(length(scrub_res$predicted_doublets) == dim(obj)[2])
    # ===================================
    scrub_res$doublet_scores_obs <- setNames(scrub_res$doublet_scores_obs, colnames(obj))
    # scrub_res$predicted_doublets <- setNames(scrub_res$predicted_doublets, colnames(obj))
    # ===================================
    obj_scrub_ls[[length(obj_scrub_ls) + 1]] <- scrub_res
    names(obj_scrub_ls)[length(obj_scrub_ls)] <- smpl
  }
  # ===================================
  # returns
  scrub.res <- new("scrubletRes",
                   obj_scrub_ls=obj_scrub_ls)
  return(scrub.res)
  # ===================================
}