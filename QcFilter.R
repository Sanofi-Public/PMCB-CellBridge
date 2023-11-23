# ===================================
# script for quality control filter
# ===================================
qcFilter <- function(obj_ls, adt_ls, scrub_ls, opt) {
  # ===================================
  sobj_flt_ls <- list()
  rna_flt_summ <- data.frame()
  if (!opt$adt) {
    adt_flt_summ <- NULL  
  } else {  
    adt_flt_summ <- data.frame()
  }
  for (i in 1:length(obj_ls)) {
    # ===================================
    obj <- obj_ls[[i]]
    smpl <- names(obj_ls)[i]
    # ===================================
    message(paste("***", smpl, "QC"))
    # message(paste("***", "ngene:", dim(obj)[1], " ncell:", dim(obj)[2]))
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
      obj <- GetAssayData(object=sobj, assay="RNA", layer="counts")
      checks <- all(
        c(sum(rowSums(obj != 0) <  opt$min_cell) == 0,
          sum(colSums(obj != 0) < opt$min_genes_per_cell) == 0,
          sum(colSums(obj != 0) > max.genes.per.cell) == 0,
          sum(colSums(obj) < opt$min_umi_per_cell) == 0,
          sum(colSums(obj) > max.umi.per.cell) == 0)
      )
    }
    # double check
    stopifnot(dim(sobj)[1] == sum(rowSums(sobj@assays$RNA$counts) != 0))
    stopifnot(dim(sobj)[2] == sum(colSums(sobj@assays$RNA$counts) != 0))
    # ===================================
    # I saw a sample with "S6T_AAACCTGAGTACCGGA-1"
    # CreateSeuratObject take "S6T" as the orig.ident, instead of smpl,
    # which will be prolematic when we are adding the metadata. So, I enforce 
    # the orig.ident below
    sobj$orig.ident <- smpl
    # ===================================
    # save counts
    rna_flt_summ <- rna_flt_summ %>%
      dplyr::bind_rows(data.frame(sample_id=smpl,
                                  type="post-qc",
                                  gene=dim(sobj)[1],
                                  cell=dim(sobj)[2]))
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
    # It is recommended to remove cells with zero reads from both the RNA and 
    # ADT (CITE) data during the quality control step to ensure data integrity and 
    # consistency in subsequent analyses. This is because cells with zero reads 
    # in the RNA data are likely to be dead or dying, and they are unlikely to 
    # be contributing to the ADT (CITE) data.
    if (opt$adt) {
      adt <- adt_ls[[i]]
      sobj[["ADT"]] <- CreateAssayObject(counts=adt[, colnames(adt) %in% colnames(sobj)])
      DefaultAssay(sobj) <- "RNA"
      adt_flt_summ <- adt_flt_summ %>%
        dplyr::bind_rows(data.frame(sample_id=smpl,
                                    type="post-qc",
                                    gene=sum(rowSums(sobj@assays$ADT$counts) != 0),
                                    cell=sum(colSums(sobj@assays$ADT$counts) != 0)))
    } 
    # ===================================
    if (!opt$adt) {
      message(paste("*** (non-zeros)", "ngene:", dim(sobj)[1], " ncell:", dim(sobj)[2]))
    } else {
      message(paste("*** RNA (non-zeros) ->", "ngene:", dim(sobj)[1], " ncell:", dim(sobj)[2]))
      message(paste("*** ADT (non-zeros) ->", 
                    "ngene:", sum(rowSums(sobj@assays$ADT$counts) != 0), 
                    " ncell:", sum(colSums(sobj@assays$ADT$counts) != 0)))
    }
    # ===================================
    sobj_flt_ls[[length(sobj_flt_ls) + 1]] <- sobj
    names(sobj_flt_ls)[length(sobj_flt_ls)] <- sobj@project.name
  }
  # ===================================
  # returns
  qcfilter.res <- new("qcfilterRes",
                      sobj_flt_ls=sobj_flt_ls,
                      rna_flt_summ=rna_flt_summ,
                      adt_flt_summ=adt_flt_summ)
  return(qcfilter.res)
  # ===================================
}