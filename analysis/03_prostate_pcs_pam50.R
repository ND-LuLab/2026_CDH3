# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  prostate_cohorts = c("GSE46691", "GSE62667", "GSE72291", "GSE79958"),
  prostate_input_dir = "data/prostate_subtypes"
)

# Supporting functions --------------------------------------------------------
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

save_figure <- function(plot, out, stem, width = 5, height = 4) {
  ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggplot2::ggsave(file.path(out, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, bg = "white")
}

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
# Topic: prostate cohort PCS/PAM50 concordance and CDH3 expression.
# Fig1D Sankey + box-plot tables; Fig1E tables; FigS1C Sankeys + tables.
# PCS-associated signatures use 37-gene singscore sets followed by
# cohort-wise score standardization. It is not the full published PCS classifier.
require_packages(c("singscore", "ggplot2", "ggalluvial", "dplyr", "svglite"))
library(ggplot2)
library(ggalluvial)
library(dplyr)
out <- output_dir("03_prostate_pcs_pam50")

# Optional upstream processing: CEL RMA and probe-to-symbol aggregation.
# The main workflow starts from normalized gene-expression matrices.
# Supply sample exclusions explicitly through exclude_sample_ids.
prepare_microarray <- function(cel_dir, probe_csv, output_csv, exclude_sample_ids = character()) {
  require_packages(c("oligo", "Biobase"))
  files <- list.files(cel_dir, pattern = "[.]CEL([.]gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No CEL files in ", cel_dir)
  raw <- oligo::read.celfiles(files)
  expr <- Biobase::exprs(oligo::rma(raw, background = TRUE, normalize = TRUE))
  if (!all(exclude_sample_ids %in% colnames(expr))) stop("Unknown excluded sample ID")
  expr <- expr[, !colnames(expr) %in% exclude_sample_ids, drop = FALSE]
  probe <- read.csv(probe_csv)
  require_columns(probe, c("ID", "gene"), "Probe dictionary")
  gene <- gsub(" ", "", probe$gene[match(rownames(expr), probe$ID)])
  keep <- !is.na(gene)
  # Average probes that map to the same gene.
  result <- aggregate(expr[keep, , drop = FALSE], by = list(gene = gene[keep]), FUN = mean)
  write.csv(result, output_csv, row.names = FALSE)
  invisible(result)
}

# 1. PCS-associated gene signatures.
sets <- list(
  PCS1 = c("STMN1", "MCM4", "CCNB1", "CDC6", "CDKN3", "EZH2", "TPX2", "FOXM1", "KIF11", "HMMR", "MKI67", "KNTC1"),
  PCS2 = c("RAB3B", "SLC4A4", "ANK3", "GJB1", "SLC12A2"),
  PCS3 = c("CFD", "COL6A1", "PTGDS", "LTBP4", "SOCS3", "SPEG", "GABRP", "PENK", "SMARCD3", "CLIP3", "ACTC1", "ASPA", "COL4A6", "CYP4B1", "ROR2", "SGCA", "SLC2A5", "PAGE4", "ACOX2", "C16orf45")
)
pam50_colors <- c(Basal = "#e7970c", LumA = "#62767f", LumB = "#97a8bc")

# 2. Per-cohort scoring, subtype comparison, and exports.
for (cohort in config$prostate_cohorts) {
  expression_path <- file.path(config$prostate_input_dir, paste0("Normalized clean ", cohort, ".csv"))
  scores_path <- file.path(config$prostate_input_dir, paste0(cohort, "_pam50scores.txt"))
  if (!all(file.exists(c(expression_path, scores_path)))) stop("Missing expression/PAM50 scores for ", cohort)
  dt <- read.csv(expression_path, check.names = FALSE)
  require_columns(dt, "gene", "Normalized expression")
  assert_ids(dt$gene, "Gene symbols")
  mat <- as.matrix(dt[, setdiff(names(dt), "gene"), drop = FALSE])
  rownames(mat) <- dt$gene
  assert_ids(colnames(mat), "Expression sample IDs")
  if (!"CDH3" %in% rownames(mat)) stop("CDH3 missing for ", cohort)
  if (!is.numeric(mat) || any(!is.finite(mat))) stop("Expression must be finite numeric values")
  ranked <- singscore::rankGenes(mat)
  scores <- sapply(sets, function(genes) singscore::simpleScore(ranked, upSet = genes)$TotalScore)
  rownames(scores) <- colnames(mat)
  standardized <- scale(scores)
  if (any(!is.finite(standardized))) stop("Constant or nonfinite PCS score")
  pcs <- colnames(scores)[apply(standardized, 1, which.max)]

  # Assign PAM50 subtypes using Basal, LumA and LumB scores.
  pam <- read.delim(scores_path, row.names = 1, check.names = FALSE)
  require_columns(pam, c("Basal", "LumA", "LumB"), "PAM50 scores")
  assert_ids(rownames(pam), "PAM50 score sample IDs")
  if (!all(colnames(mat) %in% rownames(pam))) stop("PAM50 sample IDs do not match expression")
  pam <- pam[colnames(mat), c("Basal", "LumA", "LumB"), drop = FALSE]
  if (any(!is.finite(as.matrix(pam)))) stop("Invalid PAM50 scores")
  calls <- colnames(pam)[apply(pam, 1, which.max)]
  # Align subtype assignments and expression by sample ID.
  df <- data.frame(sample = colnames(mat), CDH3 = as.numeric(mat["CDH3", ]),
                   PAM50 = calls, PCS = pcs, CDH3.zscore = as.numeric(scale(mat["CDH3", ])))
  write.csv(scores, file.path(out, paste0(cohort, "_PCS_scores.csv")))
  write.csv(standardized, file.path(out, paste0(cohort, "_PCS_scores_standardized.csv")))
  write.csv(df, file.path(out, paste0(cohort, "_CDH3_PCS_PAM50_subtype.csv")), row.names = FALSE)
  export_prism(df$CDH3.zscore, df$PAM50, out, paste0(cohort, "_CDH3_PAM50"))
  export_prism(df$CDH3.zscore, df$PCS, out, paste0(cohort, "_CDH3_PCS"))

  # 3. Draw Sankey plots and export expression tables for GraphPad box plots.
  df$PCS <- factor(df$PCS, levels = c("PCS3", "PCS2", "PCS1"))
  df_count <- df %>% count(PCS, PAM50, name = "Freq")
  write.csv(df_count, file.path(out, paste0(cohort, "_Sankey_counts.csv")), row.names = FALSE)
  p <- ggplot(df_count, aes(axis1 = PCS, axis2 = PAM50, y = Freq)) +
    geom_alluvium(aes(fill = PAM50), width = 1/12) +
    geom_stratum(width = 1/12, fill = "gray80", color = "black") +
    scale_fill_manual(values = pam50_colors) +
    scale_x_discrete(limits = c("Prostate Subtyping Classifier\n(PCS)", "Breast Cancer Model\n(PAM50)"), expand = c(.1, .1)) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4, hjust = -0.5,
              position = position_nudge(x = -0.35)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line = element_blank(), axis.text.x = element_text(size = 11),
          plot.title = element_blank(), legend.position = "right", aspect.ratio = 1) + ylab(NULL)
  save_figure(p, out, paste0(cohort, "_Sankey"), 6, 4)
  ggsave(file.path(out, paste0(cohort, "_Sankey.svg")), p, width = 5, height = 5)
}
record_session(out)
