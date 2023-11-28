# ===================================
# script to call QC pipeline functions
# ===================================
controlQC <- function(package.path, project.path, opt) {
  # ===================================
  message("** QC pipeline started")
  # ===================================
  dir.create(file.path(project.path, "outputs"), showWarnings=FALSE)
  outputs.path <- file.path(project.path, "outputs")
  # ===================================
  message("** running readIn function")
  readin.res <- readIn(project.path=project.path, 
                       package.path=package.path,
                       opt=opt)
  # ===================================
  if (!is.null(opt$scr_th)) { message("** running scrubletR function") }
  scrublet.res <- scrubletR(obj_ls=readin.res@obj_ls, 
                            opt=opt)
  # ===================================
  message("** running Report function")
  rmarkdown::render(input=file.path(package.path, "ReportQC.Rmd"),
                    output_file=paste(opt$project, "cellbridge-QC.html", sep = "_"),
                    output_dir=outputs.path, 
                    params=list(obj_ls=readin.res@obj_ls,
                                scrub_ls=scrublet.res@obj_scrub_ls,
                                meta_data=readin.res@meta_data,
                                pre.qc.summ.rna=readin.res@rna_summ,
                                opt=opt),
                    quiet=TRUE)
  # ===================================
  message("** QC pipeline completed")
  # ===================================  
}
# ===================================