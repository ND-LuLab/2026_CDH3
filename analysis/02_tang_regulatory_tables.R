# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  inputs = list(
    tang_gene_map = "data/tang/gene_map.csv",
    tang_kd_counts = "data/tang_kd/counts.tsv",
    tang_kd_samples = "data/tang_kd/sample_sheet.csv",
    tang_rlog = "data/tang/rlog.rds",
    tang_tf_list = "data/tang/dorothea_tf_list.csv"
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
# Topic: CDH3 transcription-factor correlations and YAP/TAZ knockdown; Figures 3A/I.
# Run one stage: Rscript analysis/02_tang_regulatory_tables.R correlation
#           or: Rscript analysis/02_tang_regulatory_tables.R knockdown
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[1] %in% c("correlation", "knockdown"))
  stop("Choose a stage: correlation or knockdown")
out <- output_dir(paste0("02_tang_regulatory_tables/", args[1]))

if (args[1] == "correlation") {
  # External transformed expression, gene map, and TF list are required.
  # Provide a transcription-factor list with a column named tf.
  expr <- as.matrix(readRDS(input_file("tang_rlog")))
  mapping <- read.csv(input_file("tang_gene_map"), check.names = FALSE)
  require_columns(mapping, c("ensembl_gene_id_version", "external_gene_name"), "Gene map")
  symbols <- mapping$external_gene_name[match(rownames(expr), mapping$ensembl_gene_id_version)]
  if (sum(symbols == "CDH3", na.rm = TRUE) != 1L) stop("Expected one mapped CDH3 row")
  cdh3 <- expr[which(symbols == "CDH3"), ]
  result <- t(apply(expr, 1, function(x) {
    keep <- is.finite(x) & is.finite(cdh3)
    if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(cdh3[keep]) == 0)
      return(c(correlation = NA_real_, p_value = NA_real_))
    test <- cor.test(x[keep], cdh3[keep], method = "pearson")
    c(correlation = unname(test$estimate), p_value = test$p.value)
  }))
  result <- data.frame(ensembl_id = rownames(expr), gene = symbols, result)
  tf <- read.csv(input_file("tang_tf_list"))
  require_columns(tf, "tf", "Transcription-factor list")
  tf_result <- result[result$gene %in% tf$tf, ]
  tf_result <- tf_result[order(tf_result$correlation, decreasing = TRUE, na.last = TRUE), ]
  write.csv(result, file.path(out, "CDH3_all_gene_correlations.csv"), row.names = FALSE)
  write.csv(tf_result, file.path(out, "CDH3_TF_correlations.csv"), row.names = FALSE)
  write.csv(head(tf_result[is.finite(tf_result$correlation), ], 10),
            file.path(out, "Fig3A_top10_TF_GraphPad.csv"), row.names = FALSE)
}

if (args[1] == "knockdown") {
  require_packages(c("edgeR", "limma"))
  # Assign experimental groups using the sample sheet.
  # Input count matrix: first column gene ID; remaining columns sample IDs.
  counts <- as.matrix(read.delim(input_file("tang_kd_counts"), row.names = 1, check.names = FALSE))
  sample_sheet <- read.csv(input_file("tang_kd_samples"), check.names = FALSE)
  require_columns(sample_sheet, c("sample", "group"), "Knockdown sample sheet")
  assert_ids(sample_sheet$sample, "Knockdown sample IDs")
  if (!all(sample_sheet$sample %in% colnames(counts))) stop("Unmatched knockdown sample IDs")
  counts <- counts[, sample_sheet$sample, drop = FALSE]
  group <- factor(sample_sheet$group, levels = c("scr", "TAZ_KD", "YAP_KD", "YAP_TAZ_KD"))
  if (anyNA(group) || any(table(group) == 0)) stop("Expected scr, TAZ_KD, YAP_KD and YAP_TAZ_KD groups")
  if (any(!is.finite(counts)) || any(counts < 0)) stop("Counts must be nonnegative and finite")
  x <- edgeR::DGEList(counts = counts, group = group)
  keep <- edgeR::filterByExpr(x, group = group)
  x <- x[keep, , keep.lib.sizes = FALSE]
  x <- edgeR::calcNormFactors(x, method = "TMM")
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  contrasts <- limma::makeContrasts(SCvsTAZ = TAZ_KD - scr, SCvsYAP = YAP_KD - scr,
                                   SCvsTAZ_YAP = YAP_TAZ_KD - scr, levels = design)
  v <- limma::voom(x, design, plot = FALSE)
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design), contrasts))
  logcpm <- edgeR::cpm(x, log = TRUE)
  cpm <- edgeR::cpm(x, log = FALSE)
  write.csv(cpm, file.path(out, "TMM_CPM.csv"))
  write.csv(logcpm, file.path(out, "TMM_logCPM.csv"))
  for (contrast in colnames(contrasts))
    write.csv(limma::topTable(fit, coef = contrast, number = Inf),
              file.path(out, paste0(contrast, "_limma_results.csv")))
  mapping <- read.csv(input_file("tang_gene_map"), check.names = FALSE)
  require_columns(mapping, c("ensembl_gene_id_version", "external_gene_name"), "Gene map")
  symbols <- mapping$external_gene_name[match(rownames(cpm), mapping$ensembl_gene_id_version)]
  selected <- which(symbols %in% c("CDH3", "ROR2", "WNT5A"))
  for (i in selected) {
    df <- data.frame(sample = colnames(cpm), group = group, CPM = as.numeric(cpm[i, ]),
                     logCPM = as.numeric(logcpm[i, ]))
    write.csv(df, file.path(out, paste0(symbols[i], "_sample_values.csv")), row.names = FALSE)
    export_prism(df$CPM, df$group, out, paste0(symbols[i], "_TMM_CPM"))
  }
  # Exported values are TMM CPM/logCPM; relative-to-control normalization is not applied.
  # Apply any required control normalization in GraphPad before plotting relative expression.
}
record_session(out)
