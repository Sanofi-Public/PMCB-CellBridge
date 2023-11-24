# ===================================
# script to check input arguments
# ===================================
checkOpts <- function(project.path, opt) {
  # ===================================
  if (is.null(opt$project)){
    stop("Please consider a project name.", call.=FALSE)
  }
  # ===================================
  if (is.null(opt$species)){
    stop("Please specify the species ('hs' or 'mm').", call.=FALSE)
  }
  # ===================================
  if (is.null(opt$tissue)){
    stop("Please specify the tissue.", call.=FALSE)
  }
  # ===================================
  if (is.null(opt$metadata)){
    stop("Please specify the type of metadata.", call.=FALSE)
  }
  # ===================================
  keylog <- opt$metadata %in% c("cell_based", "sample_based")
  if (!keylog) {
    msg <- paste("'metadata' must be one of 'cell_based' or 'sample_based'.")
    stop(msg)
  }
  # ===================================
  keylog <- opt$species %in% c("hs", "mm", "mf")
  if (!keylog) {
    msg <- paste("'species' must be one of 'hs', 'mm', or 'mf'.")
    stop(msg)
  }
  # ===================================
  keylog <- file.exists(file.path(project.path, "metadata.csv"))
  if (!keylog) {
    msg <- paste("'metadata.csv' not found.")
    stop(msg)
  }
  # ===================================
  keylog <- read.csv(file=file.path(project.path, "metadata.csv"), 
                     header=TRUE, nrows=1) %>%
    dplyr::rename_all(tolower)
  if ("sample_id" %in% names(keylog)) {
    msg <- paste("'sample_id' column name is reserved. Please remove or rename.")
    stop(msg)
  }
  if ("sample" %nin% names(keylog)) {
    msg <- paste("'sample' column must be included. Please add to metadata.")
    stop(msg)
  }
  if (opt$metadata == "cell_based" & "cell" %nin% names(keylog)) {
      msg <- paste("'cell' column must be included. Please add to metadata.")
      stop(msg)
  }
  # ===================================
  keylog <- opt$genesets %in% c("pbmc", "cns", "nasal", "curated", "none")
  if (!keylog) {
    msg <- paste("'genesets' must be one of 'pbmc', 'cns', 'nasal', 'curated', or 'none'.")
    stop(msg)
  }
  # ===================================
  if (opt$genesets == 'curated') {
    exls <- dir(file.path(project.path), pattern = "^genesets*", full.names = F)
    exls <- exls[grepl(".xlsx$", exls)]
    if (length(exls) > 1){
      message("Error: found more than 1 genesets file. Please keep one.")
      for (x in seq_along(exls)){
        message(paste0("*** ", x, ": ", exls[x]))
      }
      stop()      
    } else if (length(exls) == 0) {
      msg <- paste("no genesets file is found.")
      stop(msg)
    } else {
      opt$sargent <- exls
      message(paste("***", opt$sargent, "will be used for celltype annotation with Sargent."))
    }
    # ===================================
    exls <- excel_sheets(file.path(project.path, opt$sargent))
    if ("positive" %nin% exls) {
      msg <- paste("genesets must contain a sheet named 'positive'.")
      stop(msg)
    }
    # ===================================
    gpos <- read_excel(path=file.path(project.path, opt$sargent), 
                       sheet="positive", col_names=TRUE, col_types="text")
    gpos <- lapply(as.list(gpos), function(x){
      x <- x[!is.na(x)]
    })
    gpos[lengths(gpos) == 0] <- NULL
    names(gpos) <- toupper(names(gpos))
    df <- data.frame(parent=sub('[.][^.]+$', '', names(gpos)), 
                     name=names(gpos), 
                     stringsAsFactors=FALSE)
    g <- graph.data.frame(df)
    g <- igraph::simplify(g, remove.multiple=TRUE, remove.loops=TRUE, 
                          edge.attr.comb=igraph_opt("edge.attr.comb"))
    edges <- get.data.frame(g, what="edges")[1:2]
    nodes <- data.frame(id = V(g)$name, title = V(g)$name)
    for (x in nodes$title) {
      if (length(which(names(gpos) == x)) == 0) {
        msg <- paste0(x, ": parent/child incompatible. Please check the celltypes in the genesets.")
        stop(msg) 
      }
    }
    # ===================================
    if ("negative" %in% exls) {
      gneg <- read_excel(path=file.path(project.path, opt$sargent), 
                         sheet="negative", col_names=TRUE, col_types="text")
      gneg <- lapply(as.list(gneg), function(x){
        x <- x[!is.na(x)]
      })
      gneg[lengths(gneg) == 0] <- NULL
      names(gneg) <- toupper(names(gneg))
      if (!all(names(gneg) %in% names(gpos))) {
        msg <- paste0(names(gneg)[names(gneg) %nin% names(gpos)], ": all negative celltypes must exist in positive celltypes.")
        stop(msg) 
      }
    }
  }
  # ===================================
  message(paste0("** passed all checks"))
  return(opt)
  # ===================================
}
