# ===================================
# script to read in count matrices
# ===================================
readIn <- function(project.path, opt) {
  # ===================================
  if (opt$meta == "sample_based") {
    # check sample column
    meta_data_ext <- NULL
    meta_data <- read.csv(file=file.path(project.path, "metadata.csv"), 
                          header=TRUE) %>%
      dplyr::rename_all(tolower)
    if ("sample_id" %in% names(meta_data)) {
      msg <- paste("'sample_id' column name is reserved.", 
                   "Please remove or rename.", sep = "\n")
      stop(msg)
    }
    if ("sample" %nin% names(meta_data)) {
      msg <- paste("'sample' column must be included.", 
                   "Please add to metadata.", sep = "\n")
      stop(msg)
    }
    meta_data <- meta_data %>%
      dplyr::mutate(sample_id = paste0("S", 1:n()))  
  }
  # ===================================
  if (opt$meta == "cell_based") {
    meta_data_ext <- read.csv(file=file.path(project.path, "metadata.csv"), 
                              header=TRUE) %>%
      dplyr::rename_all(tolower)
    if ("sample_id" %in% names(meta_data_ext)) {
      msg <- paste("'sample_id' column name is reserved.", 
                   "Please remove or rename.", sep = "\n")
      stop(msg)
    }
    if ("sample" %nin% names(meta_data_ext)) {
      msg <- paste("'sample' column must be included.", 
                   "Please add to metadata.", sep = "\n")
      stop(msg)
    }
    if ("cell" %nin% names(meta_data_ext)) {
      msg <- paste("'cell' column must be included.", 
                   "Please add to metadata.", sep = "\n")
      stop(msg)
    }
    meta_data_ext <- meta_data_ext %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(sample_id = paste0("S", cur_group_id())) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cell = paste(sample_id, cell, sep = "_"))
    
    meta_data <- meta_data_ext
    meta_data$cell <- NULL
    meta_data <- meta_data %>%
      dplyr::group_by(sample, sample_id) %>%
      dplyr::summarise(across(everything(), .fns = function(x){
        y <- paste(unique(x), collapse = ",")
        ifelse(nchar(y) > 50, "n_char > 50", y)
      }), .groups = "drop")
    meta_data[, apply(meta_data, 2, function(x) all(x == "n_char > 50"))] <- NULL
  }
  # ===================================
  samples <- unique(meta_data$sample)
  obj_ls <- list()
  obj_summ <- data.frame()
  for (smpl in samples) {
    id <- meta_data$sample_id[meta_data$sample == smpl]
    message(paste("***", smpl, "->", id))
    filein <- list.files(file.path(project.path, smpl), recursive=FALSE, full.names=FALSE)
    sfx <- tools::file_ext(filein)
    if (unique(sfx) == "gz") {
      # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
      if (length(sfx) != 3) {
        msg <- paste("Directory should contain 3 barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz files")
        stop(msg)
      }
      ind.data <- Read10X(data.dir=file.path(project.path, smpl)) 
    }
    # ===================================
    if (unique(sfx) == "h5") {
      # Should show filtered_feature_bc_matrix.h5
      if (length(sfx) > 1) {
        msg <- paste("Directory should contain one *.h5 file")
        stop(msg)
      }
      ind.data <- Read10X_h5(filename=file.path(project.path, smpl, filein), 
                             use.names=TRUE, unique.features=TRUE)
    }
    # ===================================
    message(paste("***", "ngene:", dim(ind.data)[1], " ncell:", dim(ind.data)[2]))
    obj_summ <- obj_summ %>%
      dplyr::bind_rows(data.frame(sample_id = id,
                                  type="pre-qc",
                                  gene=sum(rowSums(ind.data) != 0),   # genes: number of rows with non-zero values
                                  cell=sum(colSums(ind.data) != 0)))  # cells: number of cols with non-zero values
    obj_ls[[length(obj_ls) + 1]] <- ind.data
    names(obj_ls)[length(obj_ls)] <- id
  }
  # ===================================
  # returns
  readin.res <- new("readinRes",
                    obj_ls=obj_ls,
                    meta_data=meta_data,
                    meta_data_ext=meta_data_ext,
                    obj_summ=obj_summ)
  return(readin.res)
  # ===================================
}