# ===================================
# Run sargent annotation 
# ===================================
sargentCellstates <- function(sobj, opt) {
  # ===================================
  if (opt$genesets == 'none') {
    sobj@misc$sargent_genesets_pos <- "none"
    sobj@misc$sargent_genesets_neg <- "none"
    # returns
    sargent.res <- new("sargentRes",
                       sobj=sobj)
    return(sargent.res)
  }
  # ===================================
  if (opt$genesets %in% c('pbmc', 'cns', 'nasal')) {
    gpos <- GENESETS[[opt$genesets]]
    gneg <- GENESETS.neg[[opt$genesets]]
  } 
  if (opt$genesets == 'curated') {
    # readin genesets
    exls <- excel_sheets(file.path(project.path, "genesets.xlsx"))
    gpos <- read_excel(path=file.path(project.path, "genesets.xlsx"), 
                       sheet="positive", col_names=TRUE, col_types="text")
    gpos <- lapply(as.list(gpos), function(x){
      x <- x[!is.na(x)]
    })
    gpos[lengths(gpos) == 0] <- NULL
    if ("negative" %in% exls) {
      gneg <- read_excel(path=file.path(project.path, "genesets.xlsx"), 
                         sheet="negative", col_names=TRUE, col_types="text")
      gneg <- lapply(as.list(gneg), function(x){
        x <- x[!is.na(x)]
      })
      gneg[lengths(gneg) == 0] <- NULL
    } else { gneg <- NULL }
  }
  # ===================================
  res <- sargentAnnotation_helper(sobj=sobj,
                                  gene.sets=gpos, 
                                  gene.sets.neg=gneg)
  # ===================================
  # label cells
  if (opt$genesets %in% c('pbmc', 'cns', 'nasal')) {  
    cellstates <- LABEL.trans[[opt$genesets]][as.character(res$cellstates)]
    names(cellstates) <- names(res$cellstates)
  } else {
    cellstates <- res$cellstates
    names(cellstates) <- names(res$cellstates)
  }
  # ===================================
  # add new metadata
  sobj <- AddMetaData(
    object = sobj,
    metadata = cellstates[Cells(sobj)],
    col.name = "sargent_cellstates"
  )
  # ===================================
  # update miscellaneous
  sobj@misc$sargent_genesets_pos <- gpos
  if (!is.null(gneg)) { sobj@misc$sargent_genesets_neg <- "none" }
  # ===================================
  # returns
  sargent.res <- new("sargentRes",
                     sobj=sobj)
  return(sargent.res)
  # ===================================
}