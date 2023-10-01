pth <- "/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_PMCB/NOURI.Nima/work/repos/cellbridge_space/cellbridge"
rmarkdown::render(input=file.path(pth, "cellbridge_flags.Rmd"),
                  output_file="cellbridge_flags.html",
                  output_dir=pth, 
                  params=list(opt_parser=opt_parser),
                  quiet=TRUE)
