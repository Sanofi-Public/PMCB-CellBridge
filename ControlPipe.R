# ===================================
# script to call pipeline functions
# ===================================
controlPipe <- function(package.path, project.path, opt, pipe_version) {
  # ===================================
  message("** pipeline started")
  # ===================================
  ui <- paste("cellbridge", 
              paste0("v", pipe_version), 
              create_unique_ids(1, char_len = 15), 
              sep = "_")
  message(paste0("** created unique id '", ui, "'"))
  # ===================================
  dir.create(file.path(project.path, "outputs"), showWarnings=FALSE)
  outputs.path <- file.path(project.path, "outputs")
  # ===================================
  message("** running readIn function")
  readin.res <- readIn(project.path=project.path, 
                       package.path=package.path,
                       opt=opt)
  # ===================================
  if (opt$scr_th != 0) { message("** running scrubletR function") }
  scrublet.res <- scrubletR(obj_ls=readin.res@obj_ls, 
                            opt=opt)
  # ===================================
  message("** running qcFilter function")
  qcfilter.res <- qcFilter(obj_ls=readin.res@obj_ls, 
                           scrub_ls=scrublet.res@obj_scrub_ls,
                           opt=opt)
  # ===================================
  message("** running mergeSamples function")
  merge.res <- mergeSamples(sobj_ls=qcfilter.res@sobj_flt_ls, 
                            meta_data=readin.res@meta_data,
                            meta_data_ext=readin.res@meta_data_ext, 
                            opt=opt)
  # ===================================
  message("** running seuratPipe function")
  seurat.res <- seuratPipe(sobj=merge.res@sobj,
                           opt=opt)
  # ===================================
  message("** running forceAtlas function")
  fatlas.res <- forceAtlas(sobj=seurat.res@sobj, 
                           opt=opt)
  # ===================================
  message("** running signacxCellstates function")
  signacx.res <- signacxCellstates(sobj=fatlas.res@sobj,
                                   package.path=package.path,
                                   opt=opt)
  # ===================================
  message("** running sargentCellstates function")
  sargent.res <- sargentCellstates(sobj=signacx.res@sobj, 
                                   opt=opt)
  # ===================================
  message("** running findSignatures function")
  signatures.res <- findSignatures(sobj=sargent.res@sobj, 
                                   opt=opt)
  # ===================================
  message("** running Report function")
  rmarkdown::render(input=file.path(package.path, "Report.Rmd"),
                    output_file=paste(opt$project, ui, "summary.html", sep = "_"),
                    output_dir=outputs.path, 
                    params=list(obj_ls=readin.res@obj_ls,
                                sobj_ls=qcfilter.res@sobj_flt_ls, 
                                scrub_ls=scrublet.res@obj_scrub_ls,
                                fsobj=signatures.res@sobj,
                                meta_data=readin.res@meta_data,
                                pre.qc.summ=readin.res@obj_summ,
                                post.qc.summ=qcfilter.res@flt_summ,
                                opt=opt,
                                ui=ui),
                    quiet=TRUE)
  # ===================================
  signatures.res@sobj@misc$args <- opt
  signatures.res@sobj@misc$identifier <- ui
  signatures.res@sobj@misc$date <- format(Sys.time(), '%d %B, %Y')
  # ===================================
  message("** saving middle files")
  res <- list(readin.res, qcfilter.res, scrublet.res, opt)
  names(res) <- c("readin.res", "qcfilter.res", "scrublet.res", "opt")
  saveRDS(res, 
          file=file.path(outputs.path, paste(opt$project, ui, "middle-object.rds", sep="_")))
  # ===================================
  message("** saving final object")
  saveRDS(signatures.res@sobj, 
          file=file.path(outputs.path, paste(opt$project, ui, "final-object.rds", sep="_")))
  # ===================================
  message("** pipeline completed")
  # ===================================  
}
# ===================================