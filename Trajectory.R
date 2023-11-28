# ===================================
# Run sargent annotation 
# ===================================
Trajectory <- function(sobj, opt) {
  # ===================================
  sds_ls <- NULL
  if (opt$trajectory != "none"){
    sds_ls <- list()
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
  }
  # ===================================
  if (!is.null(sds_ls)) {
    # Extract count matrix
    cnts_mtx <- GetAssayData(sobj, assay="RNA", layer="counts")
    # Extract cluster information
    # Idents(sobj) <- sobj$seurat_clusters
    clst_info <- Idents(sobj)
    # start edgeR
    traj_dges <- data.frame()
    for (i in seq_along(sds_ls)) {
      name_sds <- names(sds_ls)[i]
      sds <- sds_ls[[i]]
      pseudotime <- slingPseudotime(sds, na=FALSE)
      lineages <- slingLineages(sds)
      for (j in seq_along(lineages)) {
        # get lineage details
        lnge_clst <- lineages[[j]]
        lnge_nm <- names(lineages)[j]
        msg <- paste("*** running edgeR on", 
                     lnge_nm, paste0("(",j,"/", length(lineages),")"), "from",
                     name_sds, paste0("(",i,"/", length(sds_ls),")"), "...")
        message(msg)
        # Subset cells based on lineage
        cells_in_lnge <- names(clst_info)[clst_info %in% lnge_clst]
        lnge_cnts <- cnts_mtx[, cells_in_lnge]
        # Extract top variable genes
        gns <- VariableFeatures(FindVariableFeatures(CreateSeuratObject(lnge_cnts), 
                                                     selection.method="vst", 
                                                     nfeatures=opt$traj_var_gene))
        lnge_cnts <- lnge_cnts[rownames(lnge_cnts) %in% gns, ]
        # Quality Check
        lnge_cnts <- qcGex(gex=lnge_cnts, min_gene=5, min_cell=5)$gex 
        cells_in_lnge <- intersect(cells_in_lnge, colnames(lnge_cnts))
        # Design matrix and model fitting
        design <- model.matrix(~ pseudotime[cells_in_lnge, lnge_nm])
        # Check if the design matrix is of full rank
        if (!is.full.rank(design)) {
          message(paste("*** WARNING MESSAGE: Design matrix for lineage", lnge_nm, "is not of full rank. Skipping..."))
          next  # Skip to the next iteration
        }
        # Prepare data for edgeR
        dge <- DGEList(lnge_cnts)
        dge <- calcNormFactors(dge)
        dge <- estimateDisp(dge, design)
        fit <- glmQLFit(dge, design)
        # Conduct differential expression
        # > edgeR is primarily used for identifying genes whose expression levels 
        # significantly differ between two or more conditions or groups. In the 
        # context of trajectory analysis, these "conditions" can be thought of as 
        # different points or segments along a trajectory.
        # > In the results from edgeR, the log fold change (logFC) indicates how gene 
        # expression changes with pseudotime. A positive logFC indicates that a gene’s 
        # expression increases with advancing pseudotime, while a negative logFC suggests 
        # that expression decreases.
        qlf <- glmQLFTest(fit, coef=2)
        # Extract top tags (you can adjust the number of tags based on your needs)
        top_genes <- topTags(qlf, n=Inf, p.value=0.05)$table %>%  # Use n=Inf to include all genes
          tibble::rownames_to_column(var = "Gene") %>%
          dplyr::mutate(Lineage = lnge_nm, 
                        type = name_sds) 
        # Combine results
        traj_dges <- traj_dges %>%
          dplyr::bind_rows(top_genes) 
      } 
    }
    # ===================================
    sobj@misc$traj_genes <- traj_dges  
  }
  # ===================================
  # returns
  traj.res <- new("trajectoryRes",
                  sobj=sobj,
                  sds_ls=sds_ls)
  return(traj.res)
  # ===================================
}