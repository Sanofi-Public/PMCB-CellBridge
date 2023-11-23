############ let us get the bioMart conversion
shh <- suppressPackageStartupMessages
shh(library('biomaRt'))
# ===================================
message("** Building cyno monkey genome convertor")
# ===================================
`%nin%` = Negate(`%in%`)
# ===================================
message("*** Fetching ensembls...")
human <- useEnsembl(host = "https://jul2023.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
monkey <- useEnsembl(host = "https://jul2023.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mfascicularis_gene_ensembl")
message("*** Ensembls fetched.")
# ===================================
# this will only convert gene names; it will NOT convert ENSEMBL symbols
# if there is no monkey name, i doubt there is a human name.
orthologs_monkey <- getLDS(attributes = c("external_gene_name"),
                           mart = monkey, attributesL = c("hgnc_symbol"), martL = human)

colnames(orthologs_monkey) <- c("MonkeySymbol", "HumanSymbol")

orthologs_monkey[orthologs_monkey==""]<-NA
orthologs_monkey <- na.omit(orthologs_monkey)

human_nonunique <- names(table(orthologs_monkey$HumanSymbol)[table(orthologs_monkey$HumanSymbol)> 1])
monkey_nonunique <- names(table(orthologs_monkey$MonkeySymbol)[table(orthologs_monkey$MonkeySymbol)> 1])

mfTOhs <- subset(orthologs_monkey, HumanSymbol %nin% human_nonunique & MonkeySymbol %nin% monkey_nonunique)
rownames(mfTOhs) <- mfTOhs$MonkeySymbol
mfTOhs$MonkeySymbol <- NULL
# ===================================
message(paste("** saving Macaca_fascicularis.txt ->", getwd())) # -> /src/cellbridge
write.table(mfTOhs, file.path(getwd(), "Macaca_fascicularis.txt"), quote=F)
# ===================================