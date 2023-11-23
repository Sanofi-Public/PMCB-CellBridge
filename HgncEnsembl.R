############ let us get the bioMart conversion
shh <- suppressPackageStartupMessages
shh(library('biomaRt'))
# ===================================
message("*** Fetching ensembls...")
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"), verbose = TRUE)
human <- useMart(host = "https://jul2023.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# ===================================
ensTOgn <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), mart = human)
# ===================================
dups <- ensTOgn$ensembl_gene_id[which(duplicated(ensTOgn$ensembl_gene_id))]
ensTOgn <- ensTOgn[!ensTOgn$ensembl_gene_id %in% dups, ]
# ===================================
ensTOgn[ensTOgn == ""] <- NA
ensTOgn <- na.omit(ensTOgn)
# ===================================
# dups <- ensTOgn$hgnc_symbol[which(duplicated(ensTOgn$hgnc_symbol))]
# ensTOgn <- ensTOgn[!ensTOgn$hgnc_symbol %in% dups, ]
# ===================================
row.names(ensTOgn) <- ensTOgn$ensembl_gene_id
ensTOgn$ensembl_gene_id <- NULL
colnames(ensTOgn)[colnames(ensTOgn) == "hgnc_symbol"] <- "HumanSymbol"
# ===================================
message(paste("** saving Hs_ensembl.txt ->", getwd())) # -> /src/cellbridge
write.table(ensTOgn, file.path(getwd(), "Hs_ensembl.txt"), quote=F)
# ===================================