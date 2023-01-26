# ===================================
# script for quality control filter
# ===================================
qcFilter <- function(obj_ls, scrub_ls, opt) {
  # ===================================
  sobj_flt_ls <- list()
  fltr_summ <- data.frame()
  for (i in 1:length(obj_ls)) {
    # ===================================
    obj <- obj_ls[[i]]
    smpl <- names(obj_ls)[i]
    # ===================================
    message(paste("***", smpl, "QC"))
    message(paste("***", "ngene:", dim(obj)[1], " ncell:", dim(obj)[2]))
    # ===================================
    if (opt$scr_th != 0){
      dblts <- scrub_ls[[smpl]][["doublet_scores_obs"]]
      dblts <- names(dblts[dblts > opt$scr_th])
      obj <- obj[ , !(colnames(obj) %in% dblts)] 
    }
    # ===================================
    # filter low quality cells
    # itr <- 0
    if (is.null(opt$max_umi_per_cell)) { 
      max.umi.per.cell <- Inf 
    } else {
      max.umi.per.cell <- opt$max_umi_per_cell
    }
    
    if (is.null(opt$max_genes_per_cell)) { 
      max.genes.per.cell <- Inf 
    } else {
      max.genes.per.cell <- opt$max_genes_per_cell
    }
    # ===================================
    checks <- FALSE
    while (!checks) {
      # cat(paste(smpl, itr <- itr + 1), sep = "\n")
      # print(paste(dim(obj)))
      # Include genes detected in at least this many cells.
      obj <- obj[rowSums(obj != 0) >= opt$min_cell, ] 
      # print(paste(dim(obj)))
      # Include cells where at least this many genes are detected.
      obj <- obj[ , colSums(obj != 0) >= opt$min_genes_per_cell] 
      # print(paste(dim(obj)))
      # Include cells where at most this many genes are detected.
      obj <- obj[ , colSums(obj != 0) < max.genes.per.cell] 
      # print(paste(dim(obj)))
      # Include cells where at least this many transcripts are detected.
      obj <- obj[ , colSums(obj) >= opt$min_umi_per_cell] 
      # print(paste(dim(obj)))
      # Include cells where at most this many transcripts are detected.
      obj <- obj[ , colSums(obj) < max.umi.per.cell] 
      # print(paste(dim(obj)))
      # Creat a raw seurat object
      sobj <- CreateSeuratObject(counts=obj, project=smpl, 
                                 min.cells=0, min.features=0)
      # print(paste(dim(sobj)))
      # mitochondrial genes
      if (opt$species == "hs") {
        sobj[["percent_mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")  
      }
      if (opt$species == "mm") {
        sobj[["percent_mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")  
      }
      if (opt$species == "mf") {
        mito_genes <- c('ND5','COX3','ATP8','COX1','ND6','ND3','ND4L','COX2','ND1','CYTB','ATP6','ND4','ND2')
        sobj[["percent_mt"]] <- PercentageFeatureSet(sobj, pattern = mito_genes)  
      }
      # filter high percentage mitochondrial
      sobj <- subset(sobj, subset = percent_mt <= opt$max_mt_percent)
      # print(paste(dim(sobj)))
      # checks
      obj <- sobj@assays$RNA@counts
      checks <- all(
        c(sum(rowSums(obj != 0) <  opt$min_cell) == 0,
          sum(colSums(obj != 0) < opt$min_genes_per_cell) == 0,
          sum(colSums(obj != 0) > max.genes.per.cell) == 0,
          sum(colSums(obj) < opt$min_umi_per_cell) == 0,
          sum(colSums(obj) > max.umi.per.cell) == 0)
      )
    }
    # double check
    stopifnot(dim(sobj)[1] == sum(rowSums(sobj) != 0))
    stopifnot(dim(sobj)[2] == sum(colSums(sobj) != 0))
    # save counts
    fltr_summ <- fltr_summ %>%
      dplyr::bind_rows(data.frame(sample_id=smpl,
                                  type="post-qc",
                                  gene=dim(sobj)[1],
                                  cell=dim(sobj)[2]))
    # ===================================
    message(paste("***","ngene:", dim(obj)[1], " ncell:", dim(obj)[2]))
    # ===================================
    if (opt$scr_th != 0){
      y <- scrub_ls[[smpl]][["doublet_scores_obs"]]
      sobj <- AddMetaData(
        object = sobj,
        metadata = y[Cells(sobj)],
        col.name = "doublet_scores_obs"
      )
    }
    # ===================================
    sobj_flt_ls[[length(sobj_flt_ls) + 1]] <- sobj
    names(sobj_flt_ls)[length(sobj_flt_ls)] <- sobj@project.name
  }
  # ===================================
  # returns
  qcfilter.res <- new("qcfilterRes",
                      sobj_flt_ls=sobj_flt_ls,
                      flt_summ=fltr_summ)
  return(qcfilter.res)
  # ===================================
}