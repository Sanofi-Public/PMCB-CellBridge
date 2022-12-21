# ===================================
# script to call platform
# ===================================
callDocker <- function(docker) {
  if (!docker) {
    ### Run python on MAGELLAN
    # https://kb-am1.sanofi.com/display/MP/Centralized+python+environment%3A+Usage
    # Make sure to update the LD_LIBRARY_PATH for RStudio before running the script:
    # Open Magellan SSH terminal
    # Run the following command:
    # /opt/R/site_lib/RLib_Common/r-reload "/opt/py/conda/PyLib_Common/envs/scrna_doublet/lib"
    # Then re-open RStudio and test
    use_condaenv(condaenv="scrna_doublet",
                 conda="/opt/py/conda/PyLib_Common/bin/conda")
    py_config()
    use_python(system("which python", intern = TRUE))
    message(paste("** using", system("which python", intern = TRUE)))
  } else {
    ### Run python on DOCKER
    use_python(system("which python3", intern = TRUE))
    py_config()
    message(paste("** using", system("which python3", intern = TRUE)))
  }
}
# ===================================