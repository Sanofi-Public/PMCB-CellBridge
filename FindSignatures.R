# ===================================
# Find markers
# ===================================

# find markers
findSignatures <- function(sobj, opt) {
  # ===================================
  if (opt$mrk_top_n == 0) {
    sobj@misc$signatures <- "none"
    # ===================================
    # returns
    signatures.res <- new("signaturesRes",
                          sobj=sobj)
    return(signatures.res)  
  }
  # ===================================
  clusters <- c("seurat_clusters", "signacx_clusters")
  names(clusters) <- c("Seurat clusters", "Signacx clusters")
  all_markers <- data.frame()
  for (cls in as.vector(clusters)) {
    Idents(sobj) <- sobj@meta.data[[cls]]
    lvls <- as.character(unique(Idents(sobj)))
    lvls <- lvls[mixedorder(lvls)]
    Idents(sobj) <- factor(Idents(sobj), levels=lvls, ordered=TRUE)
    DefaultAssay(sobj) <- "RNA"
    message(paste("find signatures of", length(unique(sobj@meta.data[[cls]])), 
                  names(clusters[as.vector(clusters) == cls])), sep = "\n")
    markers <- FindAllMarkers(sobj, 
                              slot = "data",
                              assay = "RNA",
                              # Log fold-change of the average expression between the two groups. 
                              # Positive values indicate that the feature is more highly expressed in the first group.
                              # A threshold on base 2 log fold-change is applied
                              logfc.threshold = opt$mrk_logfc, 
                              test.use = opt$mrk_test,
                              min.pct = opt$mrk_min_pct,
                              only.pos = opt$mrk_only_pos,
                              verbose = TRUE)
    markers$type <- cls
    rownames(markers) <- NULL
    all_markers <- all_markers %>%
      dplyr::bind_rows(markers)
  }
  # ===================================
  sobj@misc$signatures <- all_markers
  # ===================================
  # returns
  signatures.res <- new("signaturesRes",
                        sobj=sobj)
  return(signatures.res)
  # ===================================
}