# ===================================
# ===================================
date <- format(Sys.time(), "%Y-%m-%d")
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
  gex <- GetAssayData(sobj, assay="RNA", slot="counts")
  adj.mtx <- attr(sobj, which="graphs")[["RNA_nn"]]
  # ===================================
  sargent_anot <- sargentAnnotation(gex=gex, 
                                    gene.sets=gsets, 
                                    gene.sets.neg=gsets.ng,
                                    adjacent.mtx=adj.mtx)
  # print(sargent_anot)
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
    # print(sargent_anot)
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
    RNA@counts@Dimnames[[1]] <- newnames
    RNA@data@Dimnames[[1]] <- newnames
    # RNA@scale.data@Dimnames[[1]] <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
# ===================================
# =================================== 