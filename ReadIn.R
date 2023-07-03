# ===================================
# script to read in count matrices
# ===================================
readIn <- function(project.path, package.path, opt) {
  # ===================================
  meta_data_ext <- NULL
  # ===================================
  if (opt$metadata == "sample_based") {
    # readin sample based metadata
    meta_data <- read.csv(file=file.path(project.path, "metadata.csv"), header=TRUE) %>%
      dplyr::rename_all(tolower) %>%
      dplyr::mutate(sample_id = paste0("S", 1:n()))  
  }
  # ===================================
  if (opt$metadata == "cell_based") {
    # readin cell based metadata
    meta_data_ext <- read.csv(file=file.path(project.path, "metadata.csv"), header=TRUE) %>%
      dplyr::rename_all(tolower) %>%
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
  if (opt$genetype == "hgnc_symbol" & opt$species == "hs") {
    xxTOhs <- read.table(file.path(package.path, "Hs_ensembl.txt"))
  }
  # ===================================
  samples <- unique(meta_data$sample)
  obj_ls <- list()
  rna_summ <- data.frame()
  # ===================================
  if (!opt$adt) {
    adt_ls <- NULL  
    adt_summ <- NULL
  } else {
    adt_ls <- list()  
    adt_summ <- data.frame()
  }
  # ===================================
  for (smpl in samples) {
    id <- meta_data$sample_id[meta_data$sample == smpl]
    message(paste("***", smpl, "->", id))
    filein <- list.files(file.path(project.path, smpl), recursive=FALSE, full.names=FALSE)
    # sfx <- tools::file_ext(filein)
    nm <- tools::file_path_sans_ext(filein, compression = TRUE)
    sfx <- sapply(1:length(nm), function(x) { gsub(paste0(nm[x], "."), "", filein[x]) })
    # sfx <- gsub(paste0(nm, "."), "", filein)
    # ===================================
    if (all(unique(sfx) %in%  c("tsv.gz", "mtx.gz"))) {
      # ===================================
      if (length(sfx) != 3) {
        msg <- paste("Directory should contain 3 barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz files")
        stop(msg)
      }
      ind.data <- Read10X(data.dir=file.path(project.path, smpl)) 
      # ===================================
    } else if (unique(sfx) == "h5") {
      if (length(sfx) > 1) {
        msg <- paste("Directory should contain one *.h5 file")
        stop(msg)
      }
      ind.data <- Read10X_h5(filename=file.path(project.path, smpl, filein), 
                             use.names=TRUE, unique.features=TRUE)
      # ===================================
    } else if (unique(sfx) == "txt.gz") {
      if (length(sfx) > 1) {
        msg <- paste("Directory should contain one *.txt.gz file")
        stop(msg)
      }
      ind.data <- data.table::fread(file = file.path(project.path, smpl, filein)) %>%
        tibble::column_to_rownames("V1")
      # ===================================
    } else if (unique(sfx) == "rds") {
      sobj <- readRDS(file.path(project.path, smpl, filein))
      stopifnot(class(sobj) == "Seurat")
      stopifnot("sample" %in% names(sobj@meta.data))
      opt.res <- opt_in_sobj(sobj=sobj, opt=opt, id=id)
      meta_data <- opt.res$meta_data
      meta_data_ext <- opt.res$meta_data_ext
      ind.data <- GetAssayData(sobj, assay="RNA", slot="counts")
      # ===================================
    } else {
      msg <- paste("unidetified input file format.")
      stop(msg)
    }
    # ===================================
    if (opt$genetype == "hgnc_symbol" & opt$species == "hs") {
      message(paste("***", "converting ensembl gene ids"))
      ind.data <- EnsemblToGnHs(ind.data, convr = xxTOhs)  
    }
    # ===================================
    # check assay
    if (!opt$adt) {
      message(paste("*** (non-zeros)", "ngene:", sum(rowSums(ind.data) != 0), " ncell:", sum(colSums(ind.data) != 0)))
      rna_summ <- rna_summ %>%
        dplyr::bind_rows(data.frame(sample_id = id,
                                    type="pre-qc",
                                    gene=sum(rowSums(ind.data) != 0),   # genes: number of rows with non-zero values
                                    cell=sum(colSums(ind.data) != 0)))  # cells: number of cols with non-zero values
      obj_ls[[length(obj_ls) + 1]] <- ind.data
      names(obj_ls)[length(obj_ls)] <- id
    } else {
      rna <- ind.data[[1]]
      adt <- ind.data[[2]] # ADT stands for Antibody-Derived Tag
      message(paste("*** RNA (non-zeros) ->", "ngene:", sum(rowSums(rna) != 0), " ncell:", sum(colSums(rna) != 0)))
      message(paste("*** ADT (non-zeros) ->", "ngene:", sum(rowSums(adt) != 0), " ncell:", sum(colSums(adt) != 0)))
      
      rna_summ <- rna_summ %>%
        dplyr::bind_rows(data.frame(sample_id = id,
                                    type="pre-qc",
                                    gene=sum(rowSums(rna) != 0),   # genes: number of rows with non-zero values
                                    cell=sum(colSums(rna) != 0)))  # cells: number of cols with non-zero values
      adt_summ <- adt_summ %>%
        dplyr::bind_rows(data.frame(sample_id = id,
                                    type="pre-qc",
                                    gene=sum(rowSums(adt) != 0),   # genes: number of rows with non-zero values
                                    cell=sum(colSums(adt) != 0)))  # cells: number of cols with non-zero values
      
      obj_ls[[length(obj_ls) + 1]] <- rna
      names(obj_ls)[length(obj_ls)] <- id
      
      adt_ls[[length(adt_ls) + 1]] <- adt
      names(adt_ls)[length(adt_ls)] <- id
    }
    # ===================================
  }
  # ===================================
  # returns
  readin.res <- new("readinRes",
                    obj_ls=obj_ls,
                    adt_ls=adt_ls,
                    meta_data=meta_data,
                    meta_data_ext=meta_data_ext,
                    rna_summ=rna_summ,
                    adt_summ=adt_summ)
  return(readin.res)
  # ===================================
}