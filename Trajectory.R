# ===================================
# Run sargent annotation 
# ===================================
Trajectory <- function(sobj, opt) {
  # ===================================
  sds_ls <- list()
  if (opt$trajectory != "none"){
    if (opt$trajectory == "umap" | opt$trajectory == "both") {
      # slingshot
      cat(paste("*** running slingshot on umap..."), sep = "\n")
      sds <- slingshot(data = Embeddings(sobj, "umap"), 
                       reducedDim = "umap",
                       clusterLabels = sobj$seurat_clusters)    
      sds_ls[[length(sds_ls) + 1]] <- sds
      names(sds_ls)[length(sds_ls)] <- paste("seurat_clusters", "umap", sep="-")
      # addining inferred lineages to the seurat object
      pt <- as.data.frame(slingPseudotime(sds))
      metas <- colnames(pt)
      for (meta in metas) {
        y <- setNames(pt[[meta]], as.character(rownames(pt)))
        # print(head(y))
        sobj <- AddMetaData(
          object = sobj,
          metadata = y,
          col.name = paste("seurat_clusters", "umap", meta, sep="-")
        )
        # print(dim(sobj@meta.data))
      } 
    }
    if ((opt$trajectory == "tsne" | opt$trajectory == "both") & opt$tsne) {
      # slingshot
      cat(paste("*** running slingshot on tsne ..."), sep = "\n")
      sds <- slingshot(data = Embeddings(sobj, "tsne"), 
                       reducedDim = "tsne",
                       clusterLabels = sobj$seurat_clusters)    
      sds_ls[[length(sds_ls) + 1]] <- sds
      names(sds_ls)[length(sds_ls)] <- paste("seurat_clusters", "tsne", sep="-")
      # addining inferred lineages to the seurat object
      pt <- as.data.frame(slingPseudotime(sds))
      metas <- colnames(pt)
      for (meta in metas) {
        y <- setNames(pt[[meta]], as.character(rownames(pt)))
        # print(head(y))
        sobj <- AddMetaData(
          object = sobj,
          metadata = y,
          col.name = paste("seurat_clusters", "tsne", meta, sep="-")
        )
        # print(dim(sobj@meta.data))
      }
    }
  } else {
    sds_ls <- NULL
  }
  # ===================================
  if (!is.null(sds_ls)) {
    # Extract top variable genes
    gns <- VariableFeatures(FindVariableFeatures(sobj, 
                                                 selection.method="vst", 
                                                 nfeatures=opt$traj_var_gene))
    # Extract count matrix
    cnts_mtx <- GetAssayData(sobj, assay="RNA", layer="counts")
    # Extract cluster information
    # Idents(sobj) <- sobj$seurat_clusters
    clst_info <- Idents(sobj)
    # start edgeR
    dge_df <- data.frame()
    for (i in seq_along(sds_ls)) {
      name_sds <- names(sds_ls)[i]
      sds <- sds_ls[[i]]
      pseudotime <- slingPseudotime(sds, na=FALSE)
      lineages <- slingLineages(sds)
      for (j in seq_along(lineages)) {
        # get lineage details
        lnge_clst <- lineages[[j]]
        lnge_nm <- names(lineages)[j]
        cat(paste("*** running edgeR on", 
                  lnge_nm, paste0("(",j,"/", length(lineages),")"), "from",
                  name_sds, paste0("(",i,"/", length(sds_ls),")"), "..."), sep = "\n")
        # Subset cells based on lineage
        cells_in_lnge <- names(clst_info)[clst_info %in% lnge_clst]
        lnge_cnts <- cnts_mtx[rownames(cnts_mtx) %in% gns, cells_in_lnge]
        # Quality Check
        lnge_cnts <- qcGex(gex=lnge_cnts, min_gene=5, min_cell=5)$gex 
        cells_in_lnge <- intersect(cells_in_lnge, colnames(lnge_cnts))
        # Prepare data for edgeR
        dge <- DGEList(lnge_cnts)
        dge <- calcNormFactors(dge)
        # Design matrix and model fitting
        design <- model.matrix(~ pseudotime[cells_in_lnge, lnge_nm])
        dge <- estimateDisp(dge, design)
        fit <- glmQLFit(dge, design)
        # Conduct differential expression
        qlf <- glmQLFTest(fit, coef=2)
        # Extract top tags (you can adjust the number of tags based on your needs)
        top_genes <- topTags(qlf, n=Inf, p.value = 0.05)$table %>%  # Use n=Inf to include all genes
          tibble::rownames_to_column(var = "Gene") %>%
          dplyr::mutate(Lineage = lnge_nm, 
                        type = name_sds) 
        # Combine results
        dge_df <- dge_df %>%
          dplyr::bind_rows(top_genes) 
      } 
    }
  }
  # ===================================
  sobj@misc$traj_genes <- dge_df
  # ===================================
  # returns
  traj.res <- new("trajectoryRes",
                  sobj=sobj,
                  sds_ls=sds_ls)
  return(traj.res)
  # ===================================
}