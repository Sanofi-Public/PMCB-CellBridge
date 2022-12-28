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
    exls <- file.exists(file.path(project.path, "genesets.xlsx"))
    if (!exls) {
      msg <- paste("'genesets.xlsx' file not found.")
      stop(msg)
    }
    exls <- excel_sheets(file.path(project.path, "genesets.xlsx"))
    if ("positive" %nin% exls) {
      msg <- paste("'genesets.xlsx' must contain a sheet named 'positive'.")
      stop(msg)
    }
  }
  # ===================================
  message(paste0("** passed all checks"))
  # ===================================
}
