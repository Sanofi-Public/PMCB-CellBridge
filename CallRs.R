# ===================================
# script to call R scripts
# ===================================
callRs <- function(package.path) {
  funcs <- c("ChecksOpts.R", "Classes.R", "Helper.R", "ReadIn.R", 
             "Scrublet.R", "QcFilter.R", "MergeSamples.R", 
             "SeuratPipe.R", "ForceAtlas.R", "SignacxCellstates.R", 
             "GeneSets.R", "SargentCellstates.R", "FindSignatures.R",
             "CallDocker.R", "ControlPipe.R", "ControlQC.R", "Trajectory.R")
  for(func in funcs) {
    source(file.path(package.path, func))
  }
  message(paste0("** sourced R scripts"))
}
# ===================================