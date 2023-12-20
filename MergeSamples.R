# ===================================
# script to merge all samples
# ===================================
mergeSamples <- function(sobj_ls, meta_data, meta_data_ext, opt) {
  # ===================================
  # rename cells using object names as prefix
  for (i in 1:length(sobj_ls)) {
    smpl <- names(sobj_ls)[i]
    sobj_ls[[i]] <- RenameCells(sobj_ls[[i]], add.cell.id=smpl)
  }
  # ===================================
  # merge all the objects in the list
  if (length(sobj_ls) > 1) {
    # sobj <- merge(sobj_ls[[1]], 
    #               y = sobj_ls[2:length(sobj_ls)], 
    #               project = opt$project)
    # sobj <- JoinLayers(sobj, assay="RNA")
    cnt_ls <- lapply(sobj_ls, function(x){
      GetAssayData(x, assay="RNA", layer="counts")
    })
    mtd_ls <- lapply(sobj_ls, function(x){
      x[[]]
    })
    mtd <- Reduce(rbind, mtd_ls)
    sobj <- CreateSeuratObject(counts=cnt_ls, meta.data=mtd, 
                               project=opt$project)
    sobj <- JoinLayers(sobj, assay="RNA")
    # ===================================
    # checks
    stopifnot(length(unique(unlist(lapply(sobj_ls, function(x) rownames(x))))) == dim(sobj)[1])
    stopifnot(length(unique(unlist(lapply(sobj_ls, function(x) colnames(x))))) == dim(sobj)[2])
    stopifnot(sum(unlist(lapply(sobj_ls, function(x) dim(x)[2]))) == dim(sobj)[2])
    stopifnot(all(colnames(sobj@assays$RNA) == colnames(sobj@assays$ADT)))
    # ===================================
    message(paste("***", "merged", length(sobj_ls), "objects."))
    # ===================================
  } else if (length(sobj_ls) == 1) {
    sobj <- sobj_ls[[1]]
    sobj@project.name <- opt$project
    # ===================================
    message(paste("***", "A single object with:"))
  } else {
    stop()
  }
  # ===================================
  DefaultAssay(sobj) <- "RNA"
  # ===================================
  if (!opt$adt) {
    message(paste("*** (non-zeros)","ngene:", dim(sobj)[1], " ncell:", dim(sobj)[2]))
  } else {
    message(paste("*** RNA (non-zeros) ->", "ngene:", dim(sobj)[1], " ncell:", dim(sobj)[2]))
    message(paste("*** ADT (non-zeros) ->", 
                  "ngene:", sum(rowSums(sobj@assays$ADT$counts) != 0), 
                  " ncell:", sum(colSums(sobj@assays$ADT$counts) != 0)))
  }
  # ===================================
  if (is.null(meta_data_ext)) {
    metas <- colnames(meta_data)
    for (meta in metas) {
      y <- setNames(meta_data[[meta]], as.character(meta_data$sample_id))
      y <- as.vector(y[sobj@meta.data$orig.ident])
      names(y) <- colnames(sobj)
      # print(head(y))
      sobj <- AddMetaData(
        object = sobj,
        metadata = y,
        col.name = meta
      )
      # print(dim(sobj@meta.data))
    } 
  } else {
    # meta_data_ext shoule be trimed since some cells have been removed in QC
    stopifnot(all(colnames(sobj) %in% meta_data_ext$cell))
    meta_data_ext <- meta_data_ext %>%
      dplyr::filter(cell %in% colnames(sobj))
    stopifnot(all(meta_data_ext$cell %in% colnames(sobj)))
    stopifnot(dim(sobj)[2] == dim(meta_data_ext)[1])
    # ===================================
    metas <- colnames(meta_data_ext)
    metas <- metas[-which(metas %in% c("cell"))]
    for (meta in metas) {
      y <- setNames(meta_data_ext[[meta]], as.character(meta_data_ext$cell))
      # print(head(y))
      sobj <- AddMetaData(
        object = sobj,
        metadata = y,
        col.name = meta
      )
      # print(dim(sobj@meta.data))
    } 
  }
  message(paste("***", "added metadata."))
  # ===================================
  # tbl <- table(sobj$sample, sobj$orig.ident)
  # tbl <- data.frame(tbl) %>%
  #   dplyr::rename(Sample = Var1,
  #                 Sample_id = Var2,
  #                 Cell = Freq) %>%
  #   dplyr::filter(Cell != 0) %>%
  #   dplyr::mutate(To = "->") %>%
  #   dplyr::select(Sample, To, Sample_id, Cell)
  # message(paste0(capture.output(tbl), collapse = "\n"))
  # ===================================
  # returns
  merge.res <- new("mergeRes",
                   sobj=sobj)
  return(merge.res)
  # ===================================
}
# ===================================
