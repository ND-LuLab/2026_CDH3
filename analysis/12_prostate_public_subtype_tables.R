# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  inputs = list(
    fhcrc_subtype_samples = "data/fhcrc_subtypes/data_clinical_sample.txt",
    fhcrc_subtype_zscores = "data/fhcrc_subtypes/data_mRNA_median_all_sample_Zscores.txt",
    tcga_prad_zscores = "data/tcga_prad/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
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
# Topic: TCGA PRAD / metastatic FHCRC PAM50-CDH3 associations; FigS1D/E.
# Stop after table export. Make the box plots in GraphPad Prism.
require_packages("genefu")
out <- output_dir("12_prostate_public_subtype_tables")
for (cohort in c("TCGA_PRAD", "FHCRC_metastatic")) {
  key <- if (cohort == "TCGA_PRAD") "tcga_prad_zscores" else "fhcrc_subtype_zscores"
  expr <- read.delim(input_file(key), check.names = FALSE)
  require_columns(expr, c("Hugo_Symbol", "Entrez_Gene_Id"), "Expression")
  if (cohort == "FHCRC_metastatic") {
    path <- input_file("fhcrc_subtype_samples")
    skip <- which(!startsWith(readLines(path), "#"))[1] - 1L
    sm <- read.delim(path, skip = skip, check.names = FALSE)
    require_columns(sm, c("SAMPLE_ID", "SAMPLE_TYPE"), "FHCRC samples")
    primary <- sm$SAMPLE_ID[sm$SAMPLE_TYPE == "Primary"]
    expr <- expr[, !names(expr) %in% primary, drop = FALSE]
  }
  # Disambiguate gene symbols and use expression z-scores for classification.
  expr <- expr[!is.na(expr$Hugo_Symbol) & expr$Hugo_Symbol != "", ]
  Gene <- make.names(expr$Hugo_Symbol, unique = TRUE)
  annot <- data.frame(Gene = Gene)
  rownames(expr) <- Gene
  mat <- t(as.matrix(expr[, setdiff(names(expr), c("Hugo_Symbol", "Entrez_Gene_Id")), drop = FALSE]))
  pred <- genefu::molecular.subtyping(sbt.model = "pam50", data = mat, annot = annot, do.mapping = FALSE)
  prob <- as.data.frame(pred$subtype.proba)
  # Restrict subtype assignment to Basal, LumA and LumB.
  require_columns(prob, c("Basal", "LumA", "LumB"), "PAM50 probabilities")
  prob <- prob[, c("Basal", "LumA", "LumB"), drop = FALSE]
  subtype <- names(prob)[apply(prob, 1, which.max)]
  if (sum(expr$Hugo_Symbol == "CDH3") != 1L) stop("Expected exactly one CDH3 expression row")
  df <- data.frame(sample = rownames(mat), CDH3 = as.numeric(mat[, "CDH3"]), PAM50 = subtype)
  write.csv(df, file.path(out, paste0(cohort, "_CDH3_source.csv")), row.names = FALSE)
  write.csv(prob, file.path(out, paste0(cohort, "_PAM50_probabilities.csv")))
  export_prism(df$CDH3, df$PAM50, out, paste0(cohort, "_CDH3_PAM50"))
}
record_session(out)
