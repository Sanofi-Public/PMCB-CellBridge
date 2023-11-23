############ let us get the bioMart conversion
shh <- suppressPackageStartupMessages
shh(library('biomaRt'))
shh(library("Orthology.eg.db"))
shh(library("org.Hs.eg.db"))
shh(library("org.Mm.eg.db"))
# ===================================
message("** Building mouse genome convertor")
# ===================================
`%nin%` = Negate(`%in%`)
# ===================================
message("*** Fetching ensembls...")
human <- useMart(host = "https://jul2023.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mouse <- useMart(host = "https://jul2023.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
message("*** Ensembls fetched.")
# ===================================
# this will only convert gene names; it will NOT convert ENSEMBL symbols
# if there is no mouse name, i doubt there is a human name.
orthologs_mouse <- getLDS(attributes = c("mgi_symbol"), mart = mouse, 
                          attributesL = c("hgnc_symbol"), martL = human)

colnames(orthologs_mouse) <- c("MouseSymbol", "HumanSymbol")

orthologs_mouse[orthologs_mouse==""]<-NA
orthologs_mouse <- na.omit(orthologs_mouse)

human_nonunique <- names(table(orthologs_mouse$HumanSymbol)[table(orthologs_mouse$HumanSymbol)> 1])
mouse_nonunique <- names(table(orthologs_mouse$MouseSymbol)[table(orthologs_mouse$MouseSymbol)> 1])

orthologs_biomart <- subset(orthologs_mouse, HumanSymbol %nin% human_nonunique & MouseSymbol %nin% mouse_nonunique)
# ===================================
# let us get the NCBI conversion
musKeys <- keys(Orthology.eg.db, "Mus.musculus")
mouse2human <- select(Orthology.eg.db, musKeys, "Homo.sapiens", "Mus.musculus")

mouse2human <- na.omit(mouse2human)

mouse_genes <- mapIds(org.Mm.eg.db, as.character(mouse2human$Mus.musculus), "SYMBOL", "ENTREZID")
human_genes <- mapIds(org.Hs.eg.db, as.character(mouse2human$Homo.sapiens), "SYMBOL", "ENTREZID")

final_df <- data.frame(MouseSymbol=as.vector(mouse_genes), HumanSymbol=as.vector(human_genes))
# ===================================
# let us remove genes that are not 1:1 according to biomart!
orthologs_ncbi <- subset(final_df, HumanSymbol %nin% human_nonunique & MouseSymbol %nin% mouse_nonunique)
additional_mouse <- setdiff(orthologs_ncbi$MouseSymbol, orthologs_biomart$MouseSymbol)
additional_human <- setdiff(orthologs_ncbi$HumanSymbol, orthologs_biomart$HumanSymbol)
additional_ortho <- subset(final_df, HumanSymbol %in% additional_human & MouseSymbol %in% additional_mouse)
# ===================================
# combine the two
mmTOhs <- rbind(orthologs_biomart, additional_ortho)
row.names(mmTOhs) <- mmTOhs$MouseSymbol
mmTOhs$MouseSymbol <- NULL
# ===================================
message(paste("** saving Mus_musculus.txt ->", getwd())) # -> /src/cellbridge
write.table(mmTOhs, file.path(getwd(), "Mus_musculus.txt"), quote=F)
# ===================================
