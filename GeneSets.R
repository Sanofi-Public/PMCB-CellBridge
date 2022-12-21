# ===================================
# GENE SETS
# ===================================
GENESETS <- list()
GENESETS.neg <- list()
LABEL.trans <- list() 
# ===================================
# PBMC GENE SETS
# ===================================
gpos <- list("TNK" = c("CD4", "CD3D", "CD3E", "CD3G", "IL7R", "CD8A", "CD8B", 
                       "FOXP3", "TIGIT", "CD27", "NCAM1", "KLRF1", "GNLY", 
                       "NKG7", "TNFRSF18"),
             # ===================================
             "TNK.T" = c("CD4", "CD3D", "CD3E", "CD3G", "IL7R", "CD8A", 
                         "CD8B", "FOXP3", "TIGIT", "CD27", "NKG7", "GNLY"),
             "TNK.NK" = c("KLRF1", "KLRC1", "KLRC2", "NCAM1", "FGFBP2", 
                          "FCGR3A", "CX3CR1", "NCR1", "FCER1G", "KDELC1", 
                          "NKG7", "GNLY", "KLRD1", "CD7"),
             # ===================================
             "TNK.T.CD4" = c("CD3D", "IL7R", "CD4"),
             "TNK.T.CD8" = c("CD3D", "IL7R", "CD8A", "CD8B"),
             # ===================================
             "TNK.T.CD4.NAIVE" = c("CCR7", "LEF1", "SELL"),
             "TNK.T.CD4.MEM" = c("S100A4", "PRDM1"),
             "TNK.T.CD4.Treg" = c("FOXP3", "TIGIT", "IL2RA", "IKZF2", "TNFRSF9"),
             # ===================================
             "TNK.T.CD8.NAIVE" = c("CCR7", "LEF1", "SELL"),
             "TNK.T.CD8.MEM" = c("S100A4", "PRDM1"),
             # ===================================
             "BPC" = c("CD79A", "CD79B", "MS4A1", "CD19", "CD24", "CD27", "CD37", 
                       "CD72", "CR2", 
                       "FCRL1", "FCRL2", "FCRL5", "FCRLA", 
                       "JCHAIN", "MZB1", "BCL11A", "SPIB", "IGHA1","IGHG1","TNFRSF17"),
             "BPC.NAIVE" = c("CD72", "CD69", "CD24", "IGHM", "IGHD"),
             "BPC.MEM" = c("CD72", "CD69", "CD24", "IGHG1", "IGHG2", "IGHG3", 
                           "IGHA1", "IGHA2","CD27", "IGHE", "MME"),
             "BPC.Plasma" = c("PRDM1", "XBP1", "JCHAIN", "MZB1"),
             # ===================================
             "MPh" = c("CD4", "CD14","LYZ","FCGR3A","MS4A7","VCAN","FCN1",   
                       "ITGAM", "MARCO", "ITGAX", "CD1C", "CST3", "FCER1A", 
                       "CLEC10A", "TCF4", "IRF7", "IL3RA", "LILRA4",
                       "JCHAIN", "BCL11A", "SPIB", "MZB1"),
             # ===================================
             "MPh.Mon" = c("CD14", "LYZ", "FCGR3A", "MS4A7", "CD163"),
             "MPh.cDC" = c("CD1C", "CST3", "FCER1A", "IRF8", "IRF7", 
                           "CLEC9A", "CLEC10A", "LILRA4"),
             "MPh.pDC" = c("CST3", "FCER1A", "TCF4", "IL3RA", "LILRA4",
                           "JCHAIN", "BCL11A", "SPIB", "MZB1"),
             # "MPh.Mac" = c("ITGAM", "ITGAX", "MARCO", "CD68", "CD86", "CD163"),
             # ===================================
             "MPh.Mon.CD14" = c("CD14", "LYZ"),
             "MPh.Mon.CD16" = c("FCGR3A", "MS4A7"),
             # ===================================
             "Megakaryocyte" = c("PPBP","PF4", "ITGA2B"),
             # ===================================
             "Erythrocyte" = 	c("HBB", "HBA1", "HBA2"),
             # ===================================
             "HSC" = c("CD34", "PROM1", "SPINK2", "HOPX", "CLEC9A", "SOX4")#,
             # ===================================
             # "Granulocyte" = c("S100A8", "S100A9", "IFITM2", "FCGR3B", 
             #                   "MS4A2", "CPA3", "TPSAB1", "FCER1A", "MPO")
)
# ===================================
gneg <- list("TNK.T.CD4" = c("CD8A", "CD8B"),
             "TNK.T.CD8" = c("CD4"),
             # ===================================
             "TNK.T.CD4.MEM" = c("CCR7", "LEF1", "SELL"),
             "TNK.T.CD4.NAIVE" = c("S100A4", "PRDM1"),
             # ===================================
             "TNK.T.CD8.MEM" = c("CCR7", "LEF1", "SELL"),
             "TNK.T.CD8.NAIVE" = c("S100A4", "PRDM1"),
             # ===================================
             "BPC.Naive" = c("MZB1", "PRDM1", "XBP1", "IGHG1", "IGHG2", 
                             "IGHG3", "IGHA1", "CD27", "MME"), # , "GPR183", "AIM2", "TNFRSF13B","SPIB", "IL10RA", "IGHE"
             "BPC.Mem" = c("MZB1", "PRDM1", "XBP1", "IGHM", "IGHD"),    
             "BPC.Plasma" = c("CD72", "CD69", "CD24"),
             # ===================================
             "MPh.Mon" = c("MARCO", 
                           "JCHAIN", "BCL11A", "SPIB", "MZB1"),
             "MPh.cDC" = c("MARCO", "CD14", "FCGR3A",
                           "JCHAIN", "BCL11A", "SPIB", "MZB1"),
             "MPh.pDC" = c("MARCO", "CD14", "FCGR3A",
                           "CD1C", "IRF8", "IRF7", "CLEC9A", "CLEC10A"),
             # ===================================
             "MPh.Mon.CD14" = c("FCGR3A", "MS4A7"),
             "MPh.Mon.CD16" = c("CD14", "LYZ")
             # ===================================
)
# ===================================
newlbls <- c("BPC.MEM" = "B.memory", 
             "BPC.NAIVE" = "B.naive", 
             "BPC.Plasma"= "PC", 
             "Erythrocyte" = "Erythrocyte", 
             "HSC" = "HSC",
             "Megakaryocyte" = "Megakaryocyte",
             "MPh.cDC" = "cDC",
             # "MPh.Mac" = "Macrophages",
             "MPh.Mon.CD14" = "Mon.Classical",
             "MPh.Mon.CD16" = "Mon.NonClassical", 
             "MPh.pDC" = "pDC",
             "TNK.NK" = "NK",
             "TNK.T.CD4.MEM" = "T.CD4.memory",
             "TNK.T.CD4.NAIVE" = "T.CD4.naive",
             "TNK.T.CD4.Treg" = "T.regs",
             "TNK.T.CD8.MEM" = "T.CD8.memory",   
             "TNK.T.CD8.NAIVE" = "T.CD8.naive",
             "TBD" = "Unclassified")
# ===================================
GENESETS[[length(GENESETS) + 1]] <- gpos
names(GENESETS)[length(GENESETS)] <- "pbmc"
GENESETS.neg[[length(GENESETS.neg) + 1]] <- gneg
names(GENESETS.neg)[length(GENESETS.neg)] <- "pbmc"
LABEL.trans[[length(LABEL.trans) + 1]] <- newlbls
names(LABEL.trans)[length(LABEL.trans)] <- "pbmc"
# ===================================
# ===================================
# CNS GENE SETS
# ===================================
# ===================================
gpos <- list("Astrocyte" = c("GFAP", "AQP4", "GJA1", "ALDH1L1", "HPSE2", "RFX4", 
                             "EAAT1", "SOX9"),
             "Endothelial" = c("FLT1", "CLDN5", "ITM2A", "EGFL7", "VWF", "ARL15", 
                               "MFSD2A", "SLC7A5", "PECAM1", "CDH5"),
             "Ependymal" = c("DNAH11", "CFAP299", "DYNLRB2", "ENKUR", "DNAH12", 
                             "CFAP43", "CCDC153", "TEKT1", "TMEM212"),
             "Fibroblast" = c("CEMIP", "ABCA10", "FBLN1", "KCNMA1", "SLC4A4", 
                              "COL1A1", "COL6A1", "LUM"),
             "Oligodendrocyte" = c("BCAS1", "MOG", "MAG", "OPALIN", "KLK6", 
                                   "CDKN1C", "HAPLN2", "MBP", "GPR37"),
             "Microglia" = c("CTSS", "CSF1R", "C1QA", "TGFBR1", "FCGR3A", "CD53", 
                             "TMEM119", "P2RY12", "CX3CR1", "TREM2", "SPI1", 
                             "AIF1", "MYO1F"),
             "Pericyte" = c("PDGFRB", "MUSTN1", "TAGLN", "ACTA2", "MYH11", "MYOCD", 
                            "HIGD1B", "NDUFA4L2", "NOTCH3", "BGN", "GRM8"),
             "Neuron" = c("SNAP25", "RBFOX3", "CHAT", "TUBB3", "DCX", "NRGN", 
                          "NRG1", "MYT1L", "SLC17A6", "SLC32A1", "GAD1", "GAD2", 
                          "PARM1", "GRIA1", "TCERG1L", "NEFL", "NEFH", "SYN1", 
                          "RBFOX1"),
             "OPC" = c("PDGFRA", "CSPG4", "OLIG1", "LHFPL3", "GPR17", "PCDH15", 
                       "VCAN"),
             "Immune" = c("CD8A", "THEMIS", "SKAP1", "STAT4", "SLAMF1", 
                          "IL7R", "NKG7", "PTPRC", "CD53")
)
x <- GENESETS[["pbmc"]]
names(x) <- paste("Immune", names(x), sep = ".")
gpos <- c(gpos, x)
# ===================================
x <- GENESETS.neg[["pbmc"]]
names(x) <- paste("Immune", names(x), sep = ".")
gneg <- x
# ===================================
newlbls <- c("Astrocyte" = "Astrocyte", 
             "Endothelial" = "Endothelial", 
             "Ependymal"= "Ependymal", 
             "Fibroblast" = "Fibroblast", 
             "Oligodendrocyte" = "Oligodendrocyte",
             "Microglia" = "Microglia",
             "Pericyte" = "Pericyte",
             "Neuron" = "Neuron",
             "OPC" = "OPC",
             "Immune" = "Immune",
             "Immune.BPC.MEM" = "B.memory", 
             "Immune.BPC.NAIVE" = "B.naive", 
             "Immune.BPC.Plasma"= "PC", 
             "Immune.Erythrocyte" = "Erythrocyte", 
             "Immune.HSC" = "HSC",
             "Immune.Megakaryocyte" = "Megakaryocyte",
             "Immune.MPh.cDC" = "cDC",
             # "Immune.MPh.Mac" = "Macrophages",
             "Immune.MPh.Mon.CD14" = "Mon.Classical",
             "Immune.MPh.Mon.CD16" = "Mon.NonClassical", 
             "Immune.MPh.pDC" = "pDC",
             "Immune.TNK.NK" = "NK",
             "Immune.TNK.T.CD4.MEM" = "T.CD4.memory",
             "Immune.TNK.T.CD4.NAIVE" = "T.CD4.naive",
             "Immune.TNK.T.CD4.Treg" = "T.regs",
             "Immune.TNK.T.CD8.MEM" = "T.CD8.memory",   
             "Immune.TNK.T.CD8.NAIVE" = "T.CD8.naive",
             "TBD" = "Unclassified")
# ===================================
GENESETS[[length(GENESETS) + 1]] <- gpos
names(GENESETS)[length(GENESETS)] <- "cns"
GENESETS.neg[[length(GENESETS.neg) + 1]] <- gneg
names(GENESETS.neg)[length(GENESETS.neg)] <- "cns"
LABEL.trans[[length(LABEL.trans) + 1]] <- newlbls
names(LABEL.trans)[length(LABEL.trans)] <- "cns"
# ===================================
# ===================================
# NASAL GENE SETS
# ===================================
# ===================================
gpos <- list("Endothelial" = c("AQP1","GNG11","ACKR1","HLA-DPB1"),
             "Pericytes" = c("NOTCH3","ACTA2","MYL9","RERGL","PDGFRB"),
             "Fibroblasts" = c("TSPAN8","GNG11","FBLN1","DCN","C1R",
                               "PDGFRB","CST3","IL32"),
             "SmoothMuscle" = c("PALLD","AQP1","C1R","DES","CNN1","ACTA2",
                                 "MYL9"),
             "Epithelial" = c("KRT5","KRT19","MUC5AC","SPDEF","SCGB1A1",
                              "SCGB3A1","BPIFB1","ANPEP","EPCAM","CAPS",
                              "TP63"),
             "Epithelial.Basal" = c("KRT5","TP63","DLK2","KRT19","CYP24A1","SERPINB4",
                         "FABP5","TFCP2L1","S100A9","MUC1","S100P","PSCA"),
             "Epithelial.Secretory" = c("KRT19","SCGB1A1","TSPAN8","MUC5AC","S100P",
                             "MUC5B","FABP5","TFCP2L1","S100A9","MUC1","PSCA",
                             "C15orf48"),
             "Epithelial.Deuterosomal" = c("SCGB1A1","MUC1","CDC20Bshort","HES6","DEUP1",
                                           "FOXJ1","PIFO"),
             "Epithelial.Multiciliated" = c("SCGB1A1","FOXJ1","PIFO","RYR3","MUC5AC",
                                            "MUC1","S100P","PSCA","TRPXL","PALLD","CYP24A1"),
             
             "Epithelial.Serous" = c("SCGB1A1","LTF","LYZ","PIP","AZGP1","SSR4"),
             "Epithelial.Mucous" = c("SCGB1A1","TSPAN8","S100P","LYZ","PIP","AZGP1",
                                     "MUC5B","BPIFB2","CST3","SSR4"),
             "Epithelial.Ionocytes" = c("FOXI1","CFTR"),
             "Immune" = c("CD8A", "THEMIS", "SKAP1", "STAT4", "SLAMF1", "IL7R", "CCR7", "S100A4", "MS4A1",
                          "NKG7", "PTPRC", "CD53")
)
x <- GENESETS[["pbmc"]]
names(x) <- paste("Immune", names(x), sep = ".")
gpos <- c(gpos, x)
# ===================================
x <- GENESETS.neg[["pbmc"]]
names(x) <- paste("Immune", names(x), sep = ".")
gneg <- x
# ===================================
newlbls <- c("Endothelial" = "Endothelial",
             "Pericytes" = "Pericytes",
             "Fibroblasts" = "Fibroblasts",
             "SmoothMuscle" = "SmoothMuscle",
             "Epithelial" = "Epithelial",
             "Epithelial.Basal" = "Basal",
             "Epithelial.Secretory" = "Secretory",
             "Epithelial.Deuterosomal" = "Deuterosomal",
             "Epithelial.Multiciliated" = "Multiciliated",
             "Epithelial.Serous" = "Serous",
             "Epithelial.Mucous" = "Mucous",
             "Epithelial.Ionocytes" = "Ionocytes",
             "Immune" = "Immune",
             "Immune.BPC.MEM" = "B.memory", 
             "Immune.BPC.NAIVE" = "B.naive", 
             "Immune.BPC.Plasma"= "PC", 
             "Immune.Erythrocyte" = "Erythrocyte", 
             "Immune.HSC" = "HSC",
             "Immune.Megakaryocyte" = "Megakaryocyte",
             "Immune.MPh.cDC" = "cDC",
             # "Immune.MPh.Mac" = "Macrophages",
             "Immune.MPh.Mon.CD14" = "Mon.Classical",
             "Immune.MPh.Mon.CD16" = "Mon.NonClassical", 
             "Immune.MPh.pDC" = "pDC",
             "Immune.TNK.NK" = "NK",
             "Immune.TNK.T.CD4.MEM" = "T.CD4.memory",
             "Immune.TNK.T.CD4.NAIVE" = "T.CD4.naive",
             "Immune.TNK.T.CD4.Treg" = "T.regs",
             "Immune.TNK.T.CD8.MEM" = "T.CD8.memory",   
             "Immune.TNK.T.CD8.NAIVE" =  "T.CD8.naive",
             "TBD" = "Unclassified")
# ===================================
GENESETS[[length(GENESETS) + 1]] <- gpos
names(GENESETS)[length(GENESETS)] <- "nasal"
GENESETS.neg[[length(GENESETS.neg) + 1]] <- gneg
names(GENESETS.neg)[length(GENESETS.neg)] <- "nasal"
LABEL.trans[[length(LABEL.trans) + 1]] <- newlbls
names(LABEL.trans)[length(LABEL.trans)] <- "nasal"
# ===================================