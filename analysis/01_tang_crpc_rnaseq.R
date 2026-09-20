# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  tang_compute_scores = TRUE,
  tang_raw_counts_confirmed = FALSE,
  tang_run_upstream = FALSE,
  inputs = list(
    tang_counts = "data/tang/GSE181374_totalDf_rlog.csv",
    tang_gene_map = "data/tang/gene_map.csv",
    tang_meta_rds = "data/tang/metaData.rds",
    tang_metadata = "data/tang/totalDf_meta.csv",
    tang_rlog = "data/tang/rlog.rds"
  )
)

# Supporting functions --------------------------------------------------------
input_file <- function(key) {
  if (!key %in% names(config$inputs)) stop("Unconfigured input: ", key)
  path <- config$inputs[[key]]
  if (!file.exists(path)) stop("Missing external input '", key, "': ", path,
                              "\nEdit the input paths at the top of this script. Data are not bundled.")
  path
}

output_dir <- function(topic) {
  path <- file.path(config$output_root, topic)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
}

require_columns <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
}

assert_ids <- function(ids, label) {
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop(label, " must contain unique, nonmissing identifiers")
}

record_session <- function(out) capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

export_prism <- function(value, group, out, stem) {
  keep <- !is.na(group)
  values <- split(value[keep], as.character(group[keep]))
  if (!length(values)) stop("No groups for GraphPad export")
  max_n <- max(lengths(values))
  wide <- as.data.frame(lapply(values, function(v) c(v, rep(NA_real_, max_n - length(v)))), check.names = FALSE)
  write.csv(wide, file.path(out, paste0(stem, "_GraphPad.csv")), row.names = FALSE, na = "")
  invisible(wide)
}

# Analysis and visualization ---------------------------------------------------
# Topic: Tang CRPC RNA-seq; Figure 1B.
# Heatmap columns follow the AR/WNT/NE/SCL model order.
# Includes Figure 1C and S1B analysis tables for GraphPad. Does not draw those panels.
require_packages(c("pheatmap"))
if (isTRUE(config$tang_run_upstream)) require_packages("DESeq2")
out <- output_dir("01_tang_crpc_rnaseq")

if (isTRUE(config$tang_run_upstream)) {
  # 1. Load and align counts and model annotations.
  counts <- as.matrix(read.csv(input_file("tang_counts"), row.names = 1, check.names = FALSE))
  meta <- read.csv(input_file("tang_metadata"), row.names = 1, check.names = FALSE)
  require_columns(meta, "group", "Tang metadata")
  assert_ids(rownames(meta), "Model IDs")
  assert_ids(colnames(counts), "Count matrix sample IDs")
  if (!all(rownames(meta) %in% colnames(counts))) stop("Metadata IDs missing from count matrix")
  counts <- counts[, rownames(meta), drop = FALSE]
  if (!isTRUE(config$tang_raw_counts_confirmed)) {
    stop("Confirm raw-count input in the settings above: input filename contains 'rlog'. ",
         "Do not pass transformed expression to DESeq2.")
  }
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0)) stop("Invalid count matrix")

  # 2. DESeq2 filtering and variance-stabilizing transformation.
  dds <- DESeq2::DESeqDataSetFromMatrix(round(counts), meta, design = ~group)
  dds <- dds[rowSums(DESeq2::counts(dds) > 0) > 10, ]
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::DESeq(dds)
  expr <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
  write.csv(data.frame(size_factor = DESeq2::sizeFactors(dds)), file.path(out, "size_factors.csv"))

} else {
  expr <- as.matrix(readRDS(input_file("tang_rlog")))
  meta <- as.data.frame(readRDS(input_file("tang_meta_rds")))
}
require_columns(meta, "group", "Tang metadata")
# 3. Map versioned Ensembl identifiers to gene symbols.
gene_map <- read.csv(input_file("tang_gene_map"), check.names = FALSE)
require_columns(gene_map, c("ensembl_gene_id_version", "external_gene_name", "chromosome_name"), "Gene map")
gene_map <- gene_map[!startsWith(gene_map$chromosome_name, "HG"), ]
groups <- list(
  AR = c("AR", "KLK2", "NKX3-1"),
  NE = c("CHGA", "SYP", "ENO2"),
  `Canonical WNT` = c("AXIN2", "TCF7L2", "CTNNB1"),
  `Non-canonical WNT` = c("WNT5A", "FZD6", "FOSL1"),
  `Stem cell` = c("CD44", "TACSTD2", "ATXN1"),
  Basal = c("TP63", "ITGB4", "KRT5", "CDH3")
)
genes <- unname(unlist(groups))
selected <- gene_map[gene_map$external_gene_name %in% genes, ]
if (anyDuplicated(selected$external_gene_name)) stop("Resolve duplicate gene mappings explicitly")
selected <- selected[match(genes, selected$external_gene_name), ]
if (anyNA(selected$ensembl_gene_id_version) || !all(selected$ensembl_gene_id_version %in% rownames(expr)))
  stop("One or more figure genes are absent from the transformed counts/gene map")

# Display order and labels for the 40 CRPC models.
sample_order <- c("MSKPCa2", "22Rv1", "C4.2", "LNCaP", "VCaP", "PDX_X0009aS1p3", "PDX_X0010aS2p2", "PDX_X0016aS2p1", "PDX_X0016bS2p1", "BM111", "ST278", "PDX_X0034aS1p4", "MSKPCa1", "MSKPCa16", "PM1078", "PM1262", "MSKEF1", "PARCB1", "PARCB3", "PARCB6", "PARCB8", "MSKPCa10", "MSKPCa14", "PM154", "H660", "ST280", "PRNEX01aS1", "MSKPCa11", "MSKPCa12", "MSKPCa13", "MSKPCa15", "MSKPCa17", "MSKPCa3", "MSKPCa8", "MSKPCa9", "PM155", "DU145", "PC3", "BM110", "ST262")
labels <- c("MSKPCa2", "22Rv1", "C4.2", "LNCaP", "VCaP", "PDX09a", "PDX10a", "PDX16a", "PDX16b", "MSKPCa19", "MSKPCa22", "PDX34a", "MSKPCa1", "MSKPCa16", "WCM1078", "WCM1262", "MSKEF1", "PARCB1", "PARCB3", "PARCB6", "PARCB8", "MSKPCa10", "MSKPCa14", "WCM154", "H660", "MSKPCa24", "PDX01a", "MSKPCa11", "MSKPCa12", "MSKPCa13", "MSKPCa15", "MSKPCa17", "MSKPCa3", "MSKPCa8", "MSKPCa9", "WCM155", "DU145", "PC3", "MSKPCa18", "MSKPCa20")
if (!all(sample_order %in% colnames(expr))) stop("Missing model in publication order")
# Standardize each gene across the 40 models in the heatmap.
if (nrow(meta) < 40L || !setequal(rownames(meta)[1:40], sample_order))
  stop("First 40 metadata models differ from the specified heatmap population")
mat <- expr[selected$ensembl_gene_id_version, rownames(meta)[1:40], drop = FALSE]
rownames(mat) <- genes
z <- t(scale(t(mat)))
if (any(!is.finite(z))) stop("Constant or non-finite expression among selected genes")
z <- z[, sample_order, drop = FALSE]
group <- as.character(meta[sample_order, "group"])
group[group == "AR_dependent"] <- "ARPC"

# Leave CDH3 outside the pathway annotation blocks.
pathway <- factor(c(rep(names(groups)[1:5], each = 3), rep("Basal", 3), NA), levels = names(groups))
row_anno <- data.frame(Pathway = pathway, row.names = genes)
col_anno <- data.frame(Group = group, row.names = sample_order)
anno_colors <- list(
  Group = c(ARPC = "#C84058", NEPC = "#517B5F", WNT = "#97A8BC", SCL = "#531A93"),
  Pathway = c(AR = "#b54d4d", NE = "#4db565", `Canonical WNT` = "#4d86b5",
              `Non-canonical WNT` = "#913e7d", `Stem cell` = "#854db5", Basal = "#f3a415")
)
# Match expression values, annotations and labels to the display order.
p <- pheatmap::pheatmap(
  z, labels_col = labels, labels_row = genes,
  annotation_col = col_anno, annotation_row = row_anno, annotation_colors = anno_colors,
  fontsize_col = 8, fontsize_row = 10, gaps_col = c(12, 16, 27), gaps_row = c(3, 6, 9),
  color = colorRampPalette(c("#053061", "#5ba2cb", "white", "#ea9173", "#6a011f"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, silent = TRUE
)
pdf(file.path(out, "Fig1B_heatmap.pdf"), width = 10, height = 4)
print(p)
dev.off()
write.csv(mat[, sample_order], file.path(out, "Fig1B_vst_expression.csv"))
write.csv(z, file.path(out, "Fig1B_row_zscores.csv"))

# 4. Calculate pathway signature scores and export GraphPad tables.
AR_sigName <- c("PTGER4","FKBP5","KLK2","CENPN","MAF","ACSL3","HERC3","ZBTB10","EAF2","ABCC4","C1orf116","PMEPA1","MED28","MPHOSPH9","TMPRSS2","KLK3","NKX3-1","NNMT","ADAM7","ELL2")
NEPC_sigName <- c("ASXL3","CAND2","ETV5","GPX2","JAKMIP2","KIAA0408","SOGA3","TRIM9","BRINP1","C7orf76","GNAO1","KCNB2","KCND2","LRRC16B","MAP10","NRSN1","PCSK1","PROX1","RGS7","SCG3","SEC11C","SEZ6","ST8SIA3","SVOP","SYT11","AURKA","DNMT1","EZH2","MYCN")
WNT_sigName <- c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E","CTNNB1","CUL1","DKK1","DKK4","DLL1","DVL2","FRAT1","FZD1","FZD8","GNAI1","HDAC11","HDAC2","HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1","MAML1","MYC","NCOR2","NCSTN","NKD1","NOTCH1","NOTCH4","NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53","WNT1","WNT5B","WNT6")
SCL_sigName <- c("ABI3BP","ACER2","ACTA2","ACTG2","ACVR2A","ADAMTS1","ADAMTS2","ADARB1","AEBP1","AGPAT4","AHI1","AKAP2","AKT3","ALDH1L2","AMOTL1","ANGPTL2","ANKRD1","ANTXR1","APOE","AQP9","ARC","ARHGAP20","ARHGAP24","ARHGAP25","ARHGEF25","ARNTL","ARSI","ASPHD2","ATP2A2","AXL","BACH1","BAG3","BCAR1","BCL2L11","BCOR","BMP1","BMP7","BNC1","BTBD11","BVES","C10orf54","C12orf24","C13orf33","C14orf149","C14orf37","C17orf81","C19orf12","C1QTNF4","C21orf7","C2orf40","C3orf54","C3orf58","C3orf64","C4orf49","C5orf25","C6orf145","C6orf228","C8orf84","C9orf3","CACNB4","CADM1","CALD1","CALML3","CALU","CARD10","CAV1","CBLB","CCDC3","CCDC85B","CCND2","CD36","CD70","CDC42EP2","CDH13","CDH3","CDH4","CDKN1A","CHST3","CHST7","CHSY3","CLIP3","CLMP","CNN1","CNP","CNRIP1","COL12A1","COL14A1","COL16A1","COL17A1","COL18A1","COL23A1","COL4A1","COL4A2","COL5A1","COL5A2","COL7A1","COL9A2","CPNE8","CPXM1","CRISPLD1","CRLF1","CRYAB","CSDC2","CSPG4","CSRP2","CTGF","CTNNAL1","CXCL14","CYGB","D4S234E","DCBLD2","DCHS1","DCUN1D3","DKK3","DLK2","DLL1","DMWD","DOCK10","DPYSL3","DST","DUSP6","DUSP7","DZIP1L","EBF3","EDARADD","EDNRB","EEPD1","EFCAB1","EFNB1","EGFR","EGR2","EGR3","EID3","ELK3","ELOVL4","ENC1","ENPP2","EPAS1","EPDR1","EPHB1","ERF","ETS1","ETV5","EVC","EXT1","FABP5","FAM101B","FAM132A","FAM176A","FAM184A","FAM70B","FAS","FBLN1","FBLN7","FBXO30","FERMT2","FEZ1","FGFRL1","FGL2","FHL1","FHOD3","FJX1","FLNC","FLRT2","FMOD","FOXP1","FST","FXYD1","FZD8","GEM","GJA1","GJC1","GNAI1","GNB4","GNG11","GOLIM4","GPC3","GPR124","GPR176","GPR3","GPR87","GPSM1","GRASP","GSN","GYLTL1B","GYPC","HAS2","HDAC4","HEG1","HGFAC","HRAS","HS3ST3A1","HSPB2","HSPG2","HTRA1","ICAM1","ID4","IGFBP2","IGFBP3","IGFBP4","IGFBP6","IL17B","IL17RD","IL1B","IL24","IL6","IL6ST","IRX4","ISM1","ITGA1","ITGA6","ITGA9","ITGB1","ITGB4","ITM2A","JAG1","JAM2","JAM3","KANK4","KCNIP3","KCNMA1","KCNMB1","KDELC1","KIAA0889","KLHDC5","KLHL21","KLHL29","KRT14","KRT16","KRT5","KRT75","LAG3","LAMA1","LAMA3","LAMB1","LAMB3","LAMC1","LBH","LCA5","LCAT","LEP","LEPRE1","LEPREL1","LGALS1","LGALS7","LGR6","LHFP","LIFR","LIMA1","LIMS2","LMOD1","LPHN1","LRCH2","LRP1","LRP4","LRRC8C","LRRN1","LTBP4","LUZP1","MALT1","MAMDC2","MAOB","MATN2","MBNL1","MCAM","MEF2C","MEG3","MEST","MFNG","MIA","MICAL2","MME","MMP2","MPDZ","MRGPRF","MRVI1","MSRB3","MSX1","MTSS1","MXRA7","MYC","MYH11","MYL9","MYLK","MYOCD","NBL1","NDN","NETO2","NGF","NGFR","NLGN2","NNAT","NNMT","NPTX2","NRCAM","NRG1","NRP1","NRP2","NT5E","NTF3","NTRK2","NUDT10","NUDT11","NXN","ODZ3","OSBPL6","OSR1","OXTR","PAMR1","PARD6G","PCBP4","PCDH18","PCDH19","PCDH7","PCDHGC3","PCOLCE","PDGFA","PDLIM4","PDLIM7","PDPN","PEG3","PELO","PGF","PHLDA3","PHLDB1","PKD1","PKD2","PKNOX2","PKP1","PLA2G7","PLCH2","PLEKHA4","PLS3","PLXNA2","PODN","POPDC2","POSTN","POU3F1","PPAP2A","PPAP2B","PPP1R14A","PPP1R16B","PPP1R18","PPP1R3C","PPP2R2B","PRDM1","PRICKLE1","PRICKLE2","PRNP","PROS1","PRRX1","PRX","PSD2","PTGS2","PTPLA","PTPRE","PTPRT","PVRL3","PXN","QKI","QRICH2","RAB34","RAPGEF1","RARB","RARRES2","RASIP1","RASL12","RBPMS","RCN3","RCSD1","RECK","RELN","RFX2","RGNEF","RHOJ","RND3","RNF165","RUSC2","SCARF2","SCHIP1","SCML2","SCN4B","SDK2","SDPR","SEC24D","SEMA3C","SEMA5A","SERPINF1","SERPING1","SERPINH1","SGCB","SGIP1","SH2D5","SH3TC1","SHE","SIAH2","SKI","SLC12A4","SLC1A3","SLC1A5","SLC25A4","SLC27A3","SLC27A6","SLC2A3","SLC38A5","SLC4A3","SLC6A8","SLCO3A1","SLIT2","SLIT3","SMTN","SNAI2","SNCA","SNTB2","SOBP","SORBS1","SORCS1","SOX11","SPARC","SPHK1","SPRED1","SRGN","SRPX","SSBP2","SSH1","STAC","STAC2","STARD8","STXBP4","SULF1","SVEP1","SYDE1","SYNM","TACC1","TAGLN","TBX2","TCF4","TCF7L1","TCOF1","TES","TGFB1I1","TGFBR3","THBS1","THSD1","THY1","TIE1","TIMP3","TINAGL1","TM7SF3","TMEM121","TMEM178","TMEM201","TMEM204","TMEM47","TMEM64","TNS1","TNS4","TOX","TP63","TPM2","TPST1","TRIM29","TRIM9","TRO","TRPC1","TSHZ2","TSHZ3","TSKU","TSPY26P","TSPYL2","TTL","TTYH2","TUBB6","TWIST2","UCN2","UNC45A","UPP1","VCAN","VGLL3","VIM","VIT","VSNL1","WIF1","WIPF1","WTIP","YAF2","ZC3H12B","ZNF219","ZNF423")
Basal_sigName <- c("TP63", "TRIM29", "ITGB4", "KRT5", "KRT14", "CDH3")
Luminal_sigName <- c("EPCAM", "KRT8", "KRT18", "HPN", "DPP4")
C_WNT_sigName <- c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E","CTNNB1","DKK1","DKK2","FZD1","FZD2","FZD3","FZD4","FZD5","GSK3B","LEF1","LRP5","LRP6","MYC","TCF7","WNT3A")
N_WNT_sigName <- c("CELSR1","DAAM1","DVL1","DVL2","FZD6","FZD7","FZD8","JNK1","JNK2","NFATC1","NFATC2","NRH2","PRICKLE1","PRICKLE2","RHOA","ROCK1","ROCK2","ROR1","ROR2","WNT5A")
YAP_sigName <- c("AJUBA","AMOTL2","AREG","CTGF","CYR61","DCHS1","FAT1","ITGB2","LATS1","LATS2","MST1","MST2","NEGR1","SAV1","TAZ","TEAD1","TEAD2","WWTR1","YAP1","ZEB1")
if (isTRUE(config$tang_compute_scores)) {
  require_packages("GSVA")
  signature_names <- c("AR", "WNT", "NEPC", "SCL", "Basal", "Luminal", "C_WNT", "N_WNT", "YAP")
  sets <- lapply(signature_names, function(s) {
    symbols <- get(paste0(s, "_sigName"))
    unique(gene_map$ensembl_gene_id_version[gene_map$external_gene_name %in% symbols])
  })
  names(sets) <- paste0(signature_names, "_score")
  # Preserve mx.diff=FALSE and gene-set sizes 5..500 across GSVA API versions.
  if (exists("gsvaParam", envir = asNamespace("GSVA"), inherits = FALSE)) {
    param <- GSVA::gsvaParam(expr, sets, minSize = 5, maxSize = 500,
                            maxDiff = FALSE, kcdf = "Gaussian")
    gs <- GSVA::gsva(param, verbose = FALSE)
  } else {
    gs <- GSVA::gsva(expr, sets, mx.diff = FALSE, min.sz = 5, max.sz = 500,
                     method = "gsva", kcdf = "Gaussian")
  }
  if (!all(names(sets) %in% rownames(gs))) stop("Some gene sets were dropped; check gene map coverage")
  meta[, colnames(t(gs))] <- t(gs)[rownames(meta), , drop = FALSE]
  write.csv(meta, file.path(out, "pathway_score.csv"))
  for (score in c("Basal_score", "C_WNT_score", "N_WNT_score")) {
    export_prism(meta[[score]], meta$group, out, paste0("Fig1C_S1B_", score))
  }
}
cdh3_table <- data.frame(sample = sample_order, group = group,
                         CDH3_vst = as.numeric(mat["CDH3", sample_order]),
                         CDH3_zscore = as.numeric(z["CDH3", ]))
write.csv(cdh3_table, file.path(out, "Fig1C_CDH3_source.csv"), row.names = FALSE)
export_prism(cdh3_table$CDH3_zscore, cdh3_table$group, out, "Fig1C_CDH3")
# Import exported tables into GraphPad Prism for the corresponding bar plots.
record_session(out)
