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
    alias <- FALSE
  } 
  # ===================================
  if (opt$genesets == 'curated') {
    # readin genesets
    exls <- excel_sheets(file.path(project.path, opt$sargent))
    # positive genes
    gpos <- read_excel(path=file.path(project.path, opt$sargent), 
                       sheet="positive", col_names=TRUE, col_types="text")
    gpos <- lapply(as.list(gpos), function(x){
      x <- x[!is.na(x)]
    })
    gpos[lengths(gpos) == 0] <- NULL
    names(gpos) <- toupper(names(gpos))
    # negative genes
    if ("negative" %in% exls) {
      gneg <- read_excel(path=file.path(project.path, opt$sargent), 
                         sheet="negative", col_names=TRUE, col_types="text")
      gneg <- lapply(as.list(gneg), function(x){
        x <- x[!is.na(x)]
      })
      gneg[lengths(gneg) == 0] <- NULL
      names(gneg) <- toupper(names(gneg))
    } else { gneg <- NULL }
    # alias
    alias <- "alias" %in% exls
  }
  # ===================================
  res <- sargentAnnotation_helper(sobj=sobj,
                                  gene.sets=gpos, 
                                  gene.sets.neg=gneg)
  # ===================================
  # label cells
  if (opt$genesets == 'curated' & alias) {
    cellstates_onto <- res$cellstates
    # add new metadata
    sobj <- AddMetaData(
      object = sobj,
      metadata = cellstates_onto[Cells(sobj)],
      col.name = "sargent_onto"
    )
    # relabeling by aliases
    message(paste("*** relabeling based on the aliases"))
    aliases <- read_excel(path=file.path(project.path, opt$sargent), 
                          sheet="alias", col_names=FALSE, col_types="text")
    colnames(aliases) <- c("old_label", "new_label") 
    aliases$old_label <- toupper(aliases$old_label)
    aliases <- setNames(aliases$new_label, aliases$old_label)
    cellstates <- res$cellstates
    cellstates_new <- ifelse(cellstates %in% names(aliases), aliases[cellstates], cellstates)
    cellstates <- setNames(cellstates_new, names(cellstates))
    sobj <- AddMetaData(
      object = sobj,
      metadata = cellstates[Cells(sobj)],
      col.name = "sargent_cellstates"
    )
  } else {
    cellstates <- res$cellstates
    # add new metadata
    sobj <- AddMetaData(
      object = sobj,
      metadata = cellstates[Cells(sobj)],
      col.name = "sargent_cellstates"
    )
    sobj$sargent_onto <- sobj$sargent_cellstates
  }
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