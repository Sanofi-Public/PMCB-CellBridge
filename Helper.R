# ===================================
# ===================================
date <- format(Sys.time(), "%Y-%m-%d")
# ===================================
# ===================================
`%nin%` = Negate(`%in%`)
# ===================================
# =================================== 
color_helper <- function(set.names, panel="Set1") {
  x <- brewer.pal.info[brewer.pal.info$category == "qual",]
  panel.size <- x$maxcolors[rownames(x) == panel]
  getPalette <- colorRampPalette(brewer.pal(panel.size, panel))  
  colrs <- setNames(getPalette(length(set.names)), 
                    set.names)
  return(colrs)
}
# ===================================
# =================================== 
getBaseTheme <- function() {
  # Define universal plot settings
  base_theme <- theme_bw() + 
    theme(text=element_text(size=10),
          plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), units="line")) +
    theme(plot.title=element_text(size=12, hjust=0.5)) +
    theme(panel.border = element_rect(size = 1.25)) +
    theme(plot.background=element_blank()#,
          # panel.grid.major=element_blank(),
          # panel.grid.minor=element_blank()
    ) +
    #theme(strip.background=element_rect(fill='white')) + 
    theme(strip.background=element_blank(),
          strip.text=element_text(size=11, face = "bold")) +
    theme(axis.title=element_text(size=11, vjust=0.25, face = "bold"),
          axis.text=element_text(size=10, face = "bold")) +
    theme(legend.margin=margin(t=-0.1, r=0, b=0, l=0, unit="line"),
          legend.position='bottom',
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          legend.key.height=grid::unit(10, "points"), 
          legend.key.width=grid::unit(15, "points"),
          legend.spacing=grid::unit(2, "points"))
  
  return(base_theme)
}
# ===================================
# =================================== 
sargentAnnotation_helper <- function(sobj, gene.sets, gene.sets.neg,
                                     nrmlz_method="LogNormalize", 
                                     scale_factor=1e6, n_hvg=500, 
                                     n_pc=30, k_param=20) {
  # ===================================
  df <- data.frame(parent=sub('[.][^.]+$', '', names(gene.sets)), 
                   name=names(gene.sets), 
                   stringsAsFactors=FALSE)
  g <- graph.data.frame(df)
  g <- igraph::simplify(g, remove.multiple=TRUE, remove.loops=TRUE, 
                        edge.attr.comb=igraph_opt("edge.attr.comb"))
  # plot(g, vertex.size=0)
  # ===================================
  # find the root nodes
  root_cltypes <- names(which(sapply(sapply(V(g), 
                                          function(x) { 
                                            neighbors(g, x, mode="in") 
                                          }), 
                                   length) == 0))
  # ===================================
  gsets <- gene.sets[names(gene.sets) %in% root_cltypes]
  gsets.ng <- gene.sets.neg[names(gene.sets.neg) %in% root_cltypes]
  if (length(gsets.ng) == 0) { gsets.ng <- NULL }
  # ===================================
  gex <- GetAssayData(sobj, assay="RNA", layer="counts")
  adj.mtx <- attr(sobj, which="graphs")[["RNA_nn"]]
  # ===================================
  sargent_anot <- sargentAnnotation(gex=gex, 
                                    gene.sets=gsets, 
                                    gene.sets.neg=gsets.ng,
                                    adjacent.mtx=adj.mtx)
  print(sargent_anot)
  # ===================================
  sargent_ls <- list()
  sargent_ls[[length(sargent_ls) + 1]] <- sargent_anot
  names(sargent_ls)[length(sargent_ls)] <- "root"
  # ===================================
  celltypes_fl <- sargent_anot@cells_type
  for (cltype in sort(names(gene.sets))) {
    # ===================================
    childs <- setdiff(names(neighborhood(g, nodes=cltype, mode="out")[[1]]), cltype)
    if (length(childs) == 0) { next }
    # cat(paste0(cltype, ": ", paste(childs, sep=", ")), sep = "\n")
    # ===================================
    cells <- celltypes_fl[[cltype]]
    if (length(cells) == 0) { next }
    # print(length(cells))
    # ===================================
    if (length(cells) > 100){
      adj.mtx <- attr(CreateSeuratObject(counts=gex) %>%
                        subset(., cell=cells) %>%
                        NormalizeData(., normalization.method=nrmlz_method, 
                                      scale.factor=scale_factor, verbose=FALSE) %>%
                        FindVariableFeatures(., selection.method="vst", nfeatures=n_hvg, 
                                             verbose=FALSE) %>%
                        ScaleData(., do.scale=TRUE, do.center=TRUE, verbose=FALSE) %>%
                        RunPCA(., features=VariableFeatures(.), verbose=FALSE) %>%
                        FindNeighbors(., reduction="pca", dims=1:n_pc, k.param=k_param, 
                                      verbose=FALSE), 
                      which="graphs")[["RNA_nn"]]  
    } else {
      adj.mtx <- NULL
    }
    # ===================================
    gsets <- gene.sets[names(gene.sets) %in% childs]
    gsets.ng <- gene.sets.neg[names(gene.sets.neg) %in% childs]
    if (length(gsets.ng) == 0) { gsets.ng <- NULL }
    # ===================================
    sargent_anot <- sargentAnnotation(gex=gex, cells=cells, 
                                      gene.sets=gsets, 
                                      gene.sets.neg=gsets.ng,
                                      adjacent.mtx=adj.mtx)
    print(sargent_anot)
    # ===================================
    names(sargent_anot@cells_type)[names(sargent_anot@cells_type) == "unclassified"] <- paste(cltype, "unspecified", sep = ".")
    # ===================================
    sargent_ls[[length(sargent_ls) + 1]] <- sargent_anot
    names(sargent_ls)[length(sargent_ls)] <- cltype
    # ===================================
    # remove parent celltype
    celltypes_fl[[cltype]] <- NULL 
    # add  progeny celltypes
    celltypes_fl <- c(celltypes_fl, sargent_anot@cells_type) 
    # merge "TBD" cells
    celltypes_fl <- base::with(utils::stack(celltypes_fl), split(values, ind)) 
    # ===================================
    # print(lengths(celltypes_fl))
  }
  # ===================================
  # print(sort(lengths(celltypes_fl)))
  cellstates <- utils::stack(celltypes_fl) %>%
    dplyr::rename(cell = values,
                  label = ind)
  cellstates <- setNames(as.character(cellstates$label), 
                         as.character(cellstates$cell))
  # print(sort(table(cellstates)))
  # ===================================
  # return
  return(list("sargent_ls" = sargent_ls,
              "cellstates" = cellstates))
}
# ===================================
# =================================== 
# if mouse, convert gene names
converta2 <- function(gene, conv) {
  ifelse(gene %in% rownames(conv), conv[gene,'HumanSymbol'], gene)
}
# ===================================
# =================================== 
# Replace gene names in different slots of a Seurat object. 
# Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
RenameGenesSeurat <- function(obj, newnames) { 
  RNA <- obj@assays$RNA
  if (nrow(RNA) == length(newnames)) {
    RNA$counts@Dimnames[[1]] <- newnames
    RNA$data@Dimnames[[1]] <- newnames
    # RNA@scale.data@Dimnames[[1]] <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
# ===================================
# =================================== 
create_unique_ids <- function(n, char_len = 5, seed_no = NULL){
  if (!is.null(seed_no)) {
    set.seed(seed_no)  
  }
  pool <- c(letters, LETTERS, 0:9)
  res <- character(n) # pre-allocating vector is much faster than growing it
  for(i in seq(n)){
    this_res <- paste0(sample(pool, char_len, replace = TRUE), collapse = "")
    while(this_res %in% res){ # if there was a duplicate, redo
      this_res <- paste0(sample(pool, char_len, replace = TRUE), collapse = "")
    }
    res[i] <- this_res
  }
  return(res)
}
# ===================================
# =================================== 
replacegrepl <- function(vec, from, to) {
  vec[grepl(paste0("^", from, "."), vec)] <- to
  return(vec)
}
# ===================================
# =================================== 
DotPlot_helper <- function(sobj, ftrs) {
  clst.idnts <- ifelse(sum(ftrs %in% rownames(sobj)) > 2 & length(unique(sobj$temp)) > 2, TRUE, FALSE)
  p <- DotPlot(sobj, group.by = "temp", features = ftrs,
               cluster.idents = clst.idnts, dot.min = 0.1) + 
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size = ifelse(length(ftrs) < 30, 10, 8)),
          axis.text.x = element_text(size = 8, angle = 25, hjust = 1, vjust = 1),
          legend.position = "none") +
    coord_flip()
  if (length(ftrs) > 30) { p <- p + scale_x_discrete(guide = guide_axis(n.dodge=2))}
  if (length(unique(sobj@meta.data[["temp"]])) > 10) { p <- p + scale_y_discrete(guide = guide_axis(n.dodge=2))}
  return(p)
}
# ===================================
# =================================== 
EnsemblToGnHs <- function(obj, convr){
  conv_gene <- sapply(rownames(obj), function(x){
    converta2(gene=x, conv=convr)
  })
  stopifnot(dim(obj)[1] == length(conv_gene))
  
  dups <- which(duplicated(conv_gene))
  dups <- which(as.vector(conv_gene) %in% as.vector(conv_gene)[dups]) 
  
  if (length(dups) > 0) {
    message(paste("*** removed", length(dups)/2, "duplicated genes"))
    conv_gene <- conv_gene[-dups]
    obj <- obj[-dups, ]
    stopifnot(dim(obj)[1] == length(conv_gene))
  }
  
  rownames(obj) <- as.vector(conv_gene)
  return(obj)
}
# ===================================
# =================================== 
update_metadata <- function(meta_data, rds_meta_data){
  # ===================================
  meta_data_ext <- lapply(rds_meta_data, function(x) {
    # x <- rds_meta_data[[1]]
    y <- filter(meta_data, sample_id == unique(x$sample_id))
    # ===================================
    logik <- which(names(x) != "sample_id")
    names(x)[logik] <- paste0("orig.", names(x)[logik])
    # ===================================
    stopifnot(intersect(names(x), names(y)) == "sample_id")
    # ===================================
    x %>%
      tibble::rownames_to_column(var = "cell") %>%
      dplyr::rename_all(tolower) %>%
      base::merge(., y, by="sample_id", all=TRUE) %>%
      dplyr::mutate(cell=paste(sample_id, cell, sep="_")) # I used orig.sample_id because in the merge section cells will be renamed the same way  
  })
  meta_data_ext <- do.call(rbind, meta_data_ext)
  # ===================================
  meta_data <- meta_data_ext
  meta_data$cell <- NULL
  meta_data <- meta_data %>%
    dplyr::group_by(sample, sample_id) %>%
    dplyr::summarise(across(everything(), .fns = function(x){
      y <- paste(unique(x), collapse = ",")
      ifelse(nchar(y) > 50, "n_char > 50", y)
    }), .groups = "drop")
  meta_data[, apply(meta_data, 2, function(x) all(x == "n_char > 50"))] <- NULL
  # ===================================
  # return
  return(list("meta_data"=meta_data,
              "meta_data_ext"=meta_data_ext))
}
# ===================================
# =================================== 
qcGex <- function(gex, min_gene=0, min_cell=0) {
  dims_i <- dim(gex)
  while (any(Matrix::colSums(gex != 0) < min_gene) | any(Matrix::rowSums(gex != 0) < min_cell)) {
    # filtered out cells for which fewer than min_gene genes were detected,
    # and genes that were expressed in fewer than min_cell cells.
    gex <- gex[, Matrix::colSums(gex != 0) >= min_gene ]
    gex <- gex[Matrix::rowSums(gex != 0) >= min_cell , ]
  }
  dims_f <- dim(gex)
  # ===================================
  message("*** ", paste(dims_i[1] - dims_f[1], "gene(s) and", dims_i[2] - dims_f[2], "cell(s) removed."))
  # ===================================
  # returns
  return(list("gex" = gex))
}
# ===================================
# ===================================
# Function to Check Full Rank -> edgeR
is.full.rank <- function(mtx) {
  rank.matrix <- qr(mtx)$rank  # Compute the rank of the matrix
  min.dim <- min(dim(mtx))
  return(rank.matrix == min.dim)
}
# ===================================
# =================================== 