# ===================================
# Run sargent annotation 
# ===================================
sargentCellstates <- function(sobj, opt) {
  # ===================================
  if (!opt$genesets) {
    gpos <- GENESETS[[opt$tissue]]
    gneg <- GENESETS.neg[[opt$tissue]]
  } else {
    # readin genesets
    exls <- excel_sheets(file.path(project.path, "genesets.xlsx"))
    if (!"positive" %in% exls) {
      msg <- paste("genesets.xlsx must contain a 'positive' sheet.")
      stop(msg)
    }
    gpos <- read_excel(path=file.path(project.path, "genesets.xlsx"), 
                       sheet="positive", col_names=TRUE, col_types="text")
    gpos <- lapply(as.list(gpos), function(x){
      x <- x[!is.na(x)]
    })
    if ("negative" %in% exls) {
      gneg <- read_excel(path=file.path(project.path, "genesets.xlsx"), 
                         sheet="negative", col_names=TRUE, col_types="text")
      gneg <- lapply(as.list(gneg), function(x){
        x <- x[!is.na(x)]
      })
    } else {
      gneg <- NULL
    }
  }
  # ===================================
  res <- sargentAnnotation_helper(sobj=sobj,
                                  gene.sets=gpos, 
                                  gene.sets.neg=gneg)
  # ===================================
  if (!opt$genesets) {  
    cellstates <- LABEL.trans[[opt$tissue]][as.character(res$cellstates)]
    names(cellstates) <- names(res$cellstates)
  } else {
    cellstates <- res$cellstates
    names(cellstates) <- names(res$cellstates)
  }
  # ===================================
  sobj <- AddMetaData(
    object = sobj,
    metadata = cellstates[Cells(sobj)],
    col.name = "sargent_cellstates"
  )
  # ===================================
  sobj@misc$sargent_genesets_pos <- gpos
  if (!is.null(gneg)) { sobj@misc$sargent_genesets_neg <- gneg }
  # ===================================
  # returns
  sargent.res <- new("sargentRes",
                     sobj=sobj)
  return(sargent.res)
  # ===================================
}