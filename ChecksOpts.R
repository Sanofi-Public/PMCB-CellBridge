# ===================================
# script to check input arguments
# ===================================
checkOpts <- function(opt) {
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
  message(paste0("** passed all checks"))
}
