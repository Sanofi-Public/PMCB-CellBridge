#!/usr/bin/env Rscript
# ===================================
### Call libraries
libs <- list("optparse", "dplyr", "Seurat", "ggplot2", "ggrepel", "pheatmap", 
             "grid", "gridExtra", "cowplot", "RColorBrewer", "reticulate", 
             "purrr", "kableExtra", "harmony", 'SignacX', 'sargent', 'igraph', 
             'gridtext', 'gplots', 'gtools', 'readxl', 'DT', 'data.tree', 'plotly')
shh <- suppressPackageStartupMessages
loads <- sapply(libs, function(x){
  shh(require(x, character.only=TRUE))
})
if (!all(loads)) {
  msg <- paste(libs[which(!loads)], "load failed.")
  stop(msg)
} else {
  message(paste("** loaded all packages"))  
}
# ===================================
option_list <- list(
  # General args
  optparse::make_option(c("--input"), type="character", default=NULL, 
                        help="Path to data directory containing one folder per sample 
                        (default=current directory)", 
                        metavar="character"),
  optparse::make_option(c("--project"), type="character", default=NULL, 
                        help="name of project 
                        (required arg)", 
                        metavar="character"),
  optparse::make_option(c("--species"), type="character", default=NULL, 
                        help="name of species (required arg) currently:
                        'hs' for humman and 'mm' for mouse", 
                        metavar="character"),
  optparse::make_option(c("--tissue"), type="character", default=NULL, 
                        help="name of tissue (required arg) currently:
                        'pbmc', 'cns', 'nasal'", 
                        metavar="character"),
  optparse::make_option(c("--metadata"), type="character", default=NULL, 
                        help="metadata type (required arg) currently:
                        'sample_based' or 'cell_based'. 
                        NOTE: column with sample info must be named 'sample',
                        NOTE: column with cell id must be named 'cell',
                        NOTE: column name 'sample_id' is reserved.", 
                        metavar="character"),
  optparse::make_option(c("--genesets"), type="logical", default=FALSE, 
                        help="curated genesets for celltype annotation with sargent,
                        (default: 'false')", 
                        metavar="logical"),
  optparse::make_option(c("--docker"), type="logical", default=TRUE, 
                        help="defines to switch setup from docker to local,
                        (default: 'true')", 
                        metavar="logical"),
  optparse::make_option(c("--only_qc"), type="logical", default=FALSE, 
                        help="set true if an overview of the quality metrics of 
                        the data is needed before running the main pipeline,
                        (default: 'false')", 
                        metavar="logical"),
  # Quality control args
  optparse::make_option(c("--min_umi_per_cell"), type="integer", default=750, 
                        help="Minimum UMI counts per cell 
                        (default=750)", 
                        metavar="integer"),
  optparse::make_option(c("--max_umi_per_cell"), type="integer", default=NULL, 
                        help="Maximum UMI counts per cell 
                        (default=NULL). If 'NULL' no upper limit will be applied", 
                        metavar="integer"),
  optparse::make_option(c("--max_mt_percent"), type="double", default=10, 
                        help="Maximum percent of mitochondrial genes 
                        (default=15%)", 
                        metavar="double"),
  optparse::make_option(c("--min_genes_per_cell"), type="integer", default=250, 
                        help="Minimum genes per cell 
                        (default=250)", 
                        metavar="integer"),
  optparse::make_option(c("--max_genes_per_cell"), type="integer", default=NULL, 
                        help="Maximum genes per cell 
                        (default=NULL). If 'NULL' no upper limit will be applied)", 
                        metavar="integer"),
  optparse::make_option(c("--min_cell"), type="integer", default=3, 
                        help="Genes detected in at least this many cells 
                        (default=3)", 
                        metavar="integer"),
  # Scrublet args
  optparse::make_option(c("--scr_th"), type="double", default=NULL, 
                        help="Score threshold for calling a transcriptome a doublet 
                        (default: null). If 'null', no doublt score calculation", 
                        metavar="double"),
  # Seurat args
  optparse::make_option(c("--seu_nrmlz_method"), type="character", default="LogNormalize", 
                        help="Method for normalization 
                        (default='LogNormalize')", 
                        metavar="character"),
  optparse::make_option(c("--seu_scale_factor"), type="double", default=1e6, 
                        help="Sets the scale factor for cell-level normalization
                        (default: 1e6)", 
                        metavar="double"),
  optparse::make_option(c("--seu_n_hvg"), type="integer", default=2000, 
                        help="Number of features to select as top variable features 
                        (default: 2000)", 
                        metavar="integer"),
  optparse::make_option(c("--seu_n_dim"), type="integer", default=30, 
                        help="Dimensions of reduction to use as input  
                        (default: 30)", 
                        metavar="integer"),
  optparse::make_option(c("--seu_k_param"), type="integer", default=20, 
                        help="Defines k for the k-nearest neighbor algorithm
                        (default: 20)", 
                        metavar="integer"),
  optparse::make_option(c("--seu_cluster_res"), type="double", default=0.7, 
                        help="Value of the resolution parameter,
                        (default: 0.7)", 
                        metavar="double"),
  optparse::make_option(c("--harmony"), type="logical", default=FALSE, 
                        help="Run harmony,
                        (default: 'false')", 
                        metavar="logical"),
  optparse::make_option(c("--harmony_var"), type="character", default="sample", 
                        help="Which variable(s) to remove. It can be multiple
                        characters separated by 'comma'),
                        (default: 'sample')", 
                        metavar="character"),
  optparse::make_option(c("--tsne"), type="logical", default=FALSE, 
                        help="Run TSNE. Time consuming for large datasets.
                        (default: 'false')", 
                        metavar="logical"),
  # cellbridge args
  optparse::make_option(c("--spr_n_dim"), type="integer", default=30, 
                        help="Dimensions of reduction to use as input  
                        (default: 30)", 
                        metavar="integer"),
  # signitures args
  optparse::make_option(c("--mrk_logfc"), type="double", default=0.25, 
                        help="Limit testing to genes which show, on average, at 
                        least X-fold difference (log-scale) between the two 
                        groups of cells.  
                        (default: 0.25)", 
                        metavar="double"),
  optparse::make_option(c("--mrk_min_pct"), type="double", default=0.50, 
                        help="Only test genes that are detected in a minimum 
                        fraction of min.pct cells in either of the two populations  
                        (default: 0.25)", 
                        metavar="double"),
  optparse::make_option(c("--mrk_only_pos"), type="logical", default=TRUE, 
                        help="Only return positive markers  
                        (default: true)", 
                        metavar="logical"),
  optparse::make_option(c("--mrk_test"), type="character", default="wilcox", 
                        help="test to use  
                        (default: wilcox)", 
                        metavar="character")
) 
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)
# print_help(opt_parser)
# ===================================
if (FALSE) {
  opt$project <- "project"
  opt$species <- "hs"
  opt$tissue <- "pbmc"
  opt$harmony <- TRUE
  opt$tsne <- TRUE
  # opt$scr_th <- 0.25
  opt$metadata <- "sample_based"
  opt$docker <- FALSE
}
# ===================================
### Rscript --vanilla read_in_samples.R /cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_PMCB/NOURI.Nima/work/repos/RP/data_for_spring > project.out 2> project.err
# ===================================
if (!opt$docker) { 
  package.path <- "/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_PMCB/NOURI.Nima/work/repos/cellbridge_space/cellbridge"
} else {
  package.path <- "/opt/cellbridge_1.0.0"
}
# ===================================
### source R functions
source(file.path(package.path, "CallRs.R"))
callRs(package.path)
# ===================================
### Check opts
checkOpts(opt)
# ===================================
### call docker
callDocker(opt$docker)
# ===================================
### source Python functions
### Note: to use source_python, the scripy has to be sourced directly, and NOT as a function
source(file.path(package.path, "CallPys.R"))
# ===================================
# project.path <- "/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_PMCB/NOURI.Nima/work/repos/cellbridge_space/cellbridge_example_proj/3gz"
# project.path <- "/cloud-data/its-cmo-darwin-bgi-virginia/Downloads/Public_Datasets/GSE174332/pipeline_input"
if (is.null(opt$input)) {
  project.path <- getwd()
} else {
  project.path <- opt$input
}
# ===================================
if (!opt$only_qc) {
  controlPipe(package.path=package.path, 
              project.path=project.path, 
              opt=opt)  
} else {
  controlQC(package.path=package.path, 
            project.path=project.path, 
            opt=opt)
}
# ===================================