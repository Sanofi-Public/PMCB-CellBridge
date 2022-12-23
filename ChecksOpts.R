# ===================================
# script to check input arguments
# ===================================
checkOpts <- function(project.path, opt) {
  if (is.null(opt$project)){
    stop("Please consider a project name.", call.=FALSE)
  }
  if (is.null(opt$species)){
    stop("Please specify the species ('hs' or 'mm').", call.=FALSE)
  }
  if (is.null(opt$tissue)){
    stop("Please specify the tissue.", call.=FALSE)
  }
  if (is.null(opt$metadata)){
    stop("Please specify the type of metadata.", call.=FALSE)
  }
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
  message(paste0("** passed all checks"))
}
