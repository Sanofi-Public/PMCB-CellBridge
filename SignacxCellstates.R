# ===================================
# Run Louvian clustering and annotation 
# based on the signacX package
# ===================================
signacxCellstates <- function(sobj, package.path, opt) {
  # generate cell type labels
  if (opt$species == "hs") {
    set.seed(42)
    cr <- Signac(E=sobj, set.seed=TRUE, seed="42", graph.used="nn")
    set.seed(42)
    crx <- GenerateLabels(cr, E=sobj, graph.used="nn")
  } 
  if (opt$species %in% c("mm", "mf")) {
    # if mouse, convert gene names to human
    if (opt$species == "mm") {
      xxTOhs <- read.table(file.path(package.path, "Mus_musculus.txt"))      
    } 
    # if macaca-fascicularis, convert gene names to human
    if (opt$species == "mf") {
      xxTOhs <- read.table(file.path(package.path, "Macaca_fascicularis.txt"))      
    } 
    sobj_tmp <- sobj
    conv_gene <- sapply(rownames(sobj_tmp), function(x){
      converta2(gene=x, conv=xxTOhs)
    })
    sobj_tmp <- RenameGenesSeurat(obj=sobj_tmp, newnames=as.vector(conv_gene))
    # generate cell type labels
    set.seed(42)
    cr <- Signac(E=sobj_tmp, set.seed=TRUE, seed="42", graph.used="nn")
    set.seed(42)
    crx <- GenerateLabels(cr, E=sobj_tmp, graph.used="nn")
  } 
  # print(table(crx$CellStates))
  # ===================================
  # sort from largest to smallest cluster
  cls <- crx$clusters
  new_rank <- table(cls)
  new_rank <- rank(dplyr::desc(new_rank), ties.method = "first") - 1
  # print(new_rank)
  # print(sort(table(cls)))
  crx$clusters <- as.character(as.vector(new_rank[match(cls, names(new_rank))]))
  # print(sort(table(crx$clusters)))
  # ===================================
  mtds <- c("clusters", "CellTypes", "CellStates")
  names(mtds) <- c("signacx_clusters", "signacx_celltypes", "signacx_cellstates")
  for (mtd in as.vector(mtds)) {
    x <- as.character(crx[[mtd]])
    names(x) <- Cells(sobj)
    sobj <- AddMetaData(
      object = sobj,
      metadata = x[Cells(sobj)],
      col.name = names(mtds)[which(mtds == mtd)]
    )
  }
  # ===================================
  # returns
  signacx.res <- new("signacxRes",
                      sobj=sobj)
  return(signacx.res)
  # ===================================
}

# https://github.com/mathewchamberlain/SignacX/blob/master/R/helper_functions.R
# data.dir = '/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_PMCB/NOURI.Nima/work/repos/RP/'
# edges = CID.LoadEdges(data.dir = data.dir) # we do not need; prodices a data frame
# head(edges)
# class(edges)


# ===================================
# edges <- sobj@misc$edges
# edges <- data.frame(do.call(rbind, edges))
# colnames(edges) <- c("V1", "V2")
# if (min(edges) == 0) { edges <- edges + 1 }
# ===================================
# # identify clusters
# set.seed(42)
# louvains <- CID.Louvain(edges=edges)
# ===================================
# # sort from largest to smallest cluster
# new_rank <- table(louvains)
# new_rank <- rank(dplyr::desc(new_rank), ties.method = "first") - 1
# # print(new_rank)
# # print(sort(table(louvains)))
# louvains <- as.character(as.vector(new_rank[match(louvains, names(new_rank))]))
# # print(sort(table(louvains)))
# ===================================
# names(louvains) <- Cells(sobj)
# sobj <- AddMetaData(
#   object = sobj,
#   metadata = louvains[Cells(sobj)],
#   col.name = "signacx_clusters"
# )
# ===================================
