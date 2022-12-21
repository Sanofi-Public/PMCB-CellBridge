# ===================================
# script to call Python scripts
# ===================================
funcs <- c("Scrublet.py", "ForceAtlas.py")
for(func in funcs) {
  source_python(file.path(package.path, func))
  }
### check Py functions
py_funs <- c("scrublet_py", "tot_counts_norm_sparse", "sparse_var", 
             "sparse_multiply", "runningquantile", "get_PCA_sparseInput",
             "get_vscores_sparse", "get_knn_graph2", "make_spring_plot")
stopifnot(all(sapply(py_funs, function(x) {
  exists(x)
  })))
stopifnot(py_available("scrublet"))
message(paste0("** sourced python scripts"))
# ===================================