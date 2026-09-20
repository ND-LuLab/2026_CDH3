# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  inputs = list(
    aurora_expr = "data/aurora/GSE193103_salmon_gene_normalized.matrix_RAP101_plus_Normals24.txt",
    fhcrc_clinical = "data/fhcrc/data_clinical_sample.txt",
    fhcrc_expr = "data/fhcrc/standardized_expression.csv"
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

record_session <- function(out) capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

# Analysis and visualization ---------------------------------------------------
# Topic: Primary/metastatic comparisons; Figure S1H (FHCRC), S1I (AURORA).
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyplots)

pal <- c(Normal = "#7FC97F", Primary = "#377EB8", Metastasis = "#E41A1C")
theme_xrot <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# RAP GSE193103
rap_expr <- input_file("aurora_expr")
out_rap <- output_dir("11_primary_metastatic_comparisons/aurora")
dir.create(out_rap, recursive = TRUE, showWarnings = FALSE)

dt <- fread(rap_expr, header = TRUE, sep = "\t", check.names = FALSE)
setnames(dt, colnames(dt)[1], "gene")
mat <- as.matrix(dt[, -1])
rownames(mat) <- dt$gene
storage.mode(mat) <- "numeric"
mat <- log2(pmax(mat, 0) + 1)

ids <- colnames(mat)
is_normal  <- grepl("(BrainN|LungN|LivN|-NL-|-Normal|^Norm)", ids, ignore.case = TRUE) | grepl("^9830-", ids)
is_primary <- grepl("PrimT|PriT|PT[Cc]ore|PrimaryT", ids) & !is_normal
is_metast  <- !is_primary & !is_normal
tt <- ifelse(is_normal, "Normal", ifelse(is_primary, "Primary", ifelse(is_metast, "Metastasis", NA_character_)))

rap <- data.table(
  sample = ids,
  patient = sub("^([A-Za-z0-9]+)-.*$", "\\1", ids),
  tissue_type = tt,
  expr = mat["CDH3", ids]
)
rap <- rap[!is.na(tissue_type)]
rap[, tissue_type := factor(tissue_type, levels = c("Normal", "Primary", "Metastasis"))]
fwrite(rap, file.path(out_rap, "r2q2_aurora_rap_3col.csv"))

n <- rap[, .N, by = tissue_type][order(match(tissue_type, levels(rap$tissue_type)))]
lab <- setNames(paste0(n$tissue_type, "\n(n=", n$N, ")"), as.character(n$tissue_type))
du <- copy(rap)
du[, x := factor(tissue_type, levels = names(lab), labels = lab)]
yr <- range(du$expr, na.rm = TRUE)
p_np <- suppressWarnings(wilcox.test(du[tissue_type == "Normal", expr], du[tissue_type == "Primary", expr], exact = FALSE)$p.value)
p_nm <- suppressWarnings(wilcox.test(du[tissue_type == "Normal", expr], du[tissue_type == "Metastasis", expr], exact = FALSE)$p.value)
p_pm <- suppressWarnings(wilcox.test(du[tissue_type == "Primary", expr], du[tissue_type == "Metastasis", expr], exact = FALSE)$p.value)
annot <- ifelse(c(p_np, p_nm, p_pm) < 0.001, "p<0.001", paste0("p=", formatC(c(p_np, p_nm, p_pm), digits = 3, format = "g")))
step <- 0.13 * diff(yr)
ypos <- yr[2] + step * 1:3
p_un <- du %>%
  tidyplot(x = x, y = expr, color = tissue_type) %>%
  add_boxplot(alpha = 0, linewidth = 0.4) %>%
  adjust_colors(pal) %>%
  adjust_x_axis_title("") %>%
  adjust_y_axis_title("log2 normalized CDH3 (salmon)") %>%
  remove_legend() %>%
  add(theme_xrot) %>%
  add_data_points_jitter(jitter_width = 0.18, alpha = 0.85, size = 0.9) %>%
  add(geom_signif(
    comparisons = list(c(lab[["Normal"]], lab[["Primary"]]),
                       c(lab[["Normal"]], lab[["Metastasis"]]),
                       c(lab[["Primary"]], lab[["Metastasis"]])),
    annotations = annot, y_position = ypos, tip_length = 0.005,
    textsize = 2.8, vjust = -0.3, color = "black", size = 0.3, fontface = "plain"
  )) %>%
  add(expand_limits(y = max(ypos) + step * 1.5)) %>%
  add_caption("Wilcoxon rank-sum (unpaired)") %>%
  adjust_size(width = 50, height = 55)
save_plot(p_un, file.path(out_rap, "r2q2_aurora_rap_3col_overall_unpaired.pdf"), units = "mm", view_plot = FALSE, padding = 0.05, device = cairo_pdf)
save_plot(p_un, file.path(out_rap, "r2q2_aurora_rap_3col_overall_unpaired.png"), units = "mm", view_plot = FALSE, padding = 0.05, dpi = 600)

multi <- rap[, uniqueN(tissue_type), by = patient][V1 >= 2L, patient]
dp <- rap[patient %in% multi]
n <- dp[, .N, by = tissue_type][order(match(tissue_type, levels(dp$tissue_type)))]
lab <- setNames(paste0(n$tissue_type, "\n(n=", n$N, ")"), as.character(n$tissue_type))
dp[, x := factor(tissue_type, levels = names(lab), labels = lab)]
yr <- range(dp$expr, na.rm = TRUE)
nrm <- dp[tissue_type == "Normal", .(y = mean(expr, na.rm = TRUE)), by = patient]
pri <- dp[tissue_type == "Primary", .(y = mean(expr, na.rm = TRUE)), by = patient]
met <- dp[tissue_type == "Metastasis", .(y = mean(expr, na.rm = TRUE)), by = patient]
pid_np <- intersect(nrm$patient, pri$patient)
pid_nm <- intersect(nrm$patient, met$patient)
pid_pm <- intersect(pri$patient, met$patient)
comps <- list(); annot <- character()
if (length(pid_np) >= 3L) {
  p_np <- suppressWarnings(wilcox.test(nrm[match(pid_np, patient)]$y, pri[match(pid_np, patient)]$y, paired = TRUE, exact = FALSE)$p.value)
  comps[[length(comps) + 1L]] <- c(lab[["Normal"]], lab[["Primary"]])
  annot <- c(annot, if (p_np < 0.001) "p<0.001" else paste0("p=", formatC(p_np, digits = 3, format = "g")))
}
if (length(pid_nm) >= 3L) {
  p_nm <- suppressWarnings(wilcox.test(nrm[match(pid_nm, patient)]$y, met[match(pid_nm, patient)]$y, paired = TRUE, exact = FALSE)$p.value)
  comps[[length(comps) + 1L]] <- c(lab[["Normal"]], lab[["Metastasis"]])
  annot <- c(annot, if (p_nm < 0.001) "p<0.001" else paste0("p=", formatC(p_nm, digits = 3, format = "g")))
}
if (length(pid_pm) >= 3L) {
  p_pm <- suppressWarnings(wilcox.test(pri[match(pid_pm, patient)]$y, met[match(pid_pm, patient)]$y, paired = TRUE, exact = FALSE)$p.value)
  comps[[length(comps) + 1L]] <- c(lab[["Primary"]], lab[["Metastasis"]])
  annot <- c(annot, if (p_pm < 0.001) "p<0.001" else paste0("p=", formatC(p_pm, digits = 3, format = "g")))
}
step <- 0.13 * diff(yr)
ypos <- yr[2] + step * seq_along(annot)
p_pr <- dp %>%
  tidyplot(x = x, y = expr, color = tissue_type) %>%
  add_boxplot(alpha = 0, linewidth = 0.4) %>%
  adjust_colors(pal) %>%
  adjust_x_axis_title("") %>%
  adjust_y_axis_title("log2 normalized CDH3 (salmon)") %>%
  remove_legend() %>%
  add(theme_xrot) %>%
  add(geom_line(aes(group = patient), color = "grey75", alpha = 0.5, linewidth = 0.3, inherit.aes = TRUE)) %>%
  add_data_points(alpha = 0.9, size = 0.9) %>%
  add(geom_signif(
    comparisons = comps, annotations = annot, y_position = ypos, tip_length = 0.005,
    textsize = 2.8, vjust = -0.3, color = "black", size = 0.3, fontface = "plain"
  )) %>%
  add(expand_limits(y = max(ypos) + step * 1.5)) %>%
  add_caption("Wilcoxon signed-rank (paired, on patient-mean)") %>%
  adjust_size(width = 50, height = 55)
save_plot(p_pr, file.path(out_rap, "r2q2_aurora_rap_3col_overall_paired.pdf"), units = "mm", view_plot = FALSE, padding = 0.05, device = cairo_pdf)
save_plot(p_pr, file.path(out_rap, "r2q2_aurora_rap_3col_overall_paired.png"), units = "mm", view_plot = FALSE, padding = 0.05, dpi = 600)

# FHCRC Kumar 2016
fhc_expr <- input_file("fhcrc_expr")
fhc_clin <- input_file("fhcrc_clinical")
out_fhc <- output_dir("11_primary_metastatic_comparisons/fhcrc")
dir.create(out_fhc, recursive = TRUE, showWarnings = FALSE)

expr <- fread(fhc_expr, header = TRUE, check.names = FALSE)
gcol <- intersect(c("GeneSymbol", "Gene_symbol", "gene_name", "Hugo_Symbol"), colnames(expr))[1]
scols <- setdiff(colnames(expr), c(gcol, "Entrez", "Entrez_Gene_Id"))
cdh3 <- as.numeric(expr[which(expr[[gcol]] == "CDH3")[1], scols, with = FALSE])
names(cdh3) <- scols

skip <- 0L
for (ln in readLines(fhc_clin)) {
  if (!startsWith(ln, "#")) break
  skip <- skip + 1L
}
cl <- fread(fhc_clin, sep = "\t", skip = skip, header = TRUE, na.strings = c("", "NA", "[Not Available]"))
cl <- cl[, .(SAMPLE_ID = as.character(SAMPLE_ID), PATIENT_ID = as.character(PATIENT_ID),
             SAMPLE_TYPE = as.character(SAMPLE_TYPE))]

fhc <- data.table(sample = names(cdh3), expr = unname(cdh3))
fhc <- merge(fhc, cl, by.x = "sample", by.y = "SAMPLE_ID")
fhc <- fhc[SAMPLE_TYPE %in% c("Primary", "Metastasis")]
fhc[, tissue_type := factor(SAMPLE_TYPE, levels = c("Primary", "Metastasis"))]
setnames(fhc, "PATIENT_ID", "patient")
fhc <- fhc[, .(sample, patient, tissue_type, expr)]
fwrite(fhc, file.path(out_fhc, "r2q2_fhcrc_2col.csv"))

n <- fhc[, .N, by = tissue_type][order(match(tissue_type, levels(fhc$tissue_type)))]
lab <- setNames(paste0(n$tissue_type, "\n(n=", n$N, ")"), as.character(n$tissue_type))
du <- copy(fhc)
du[, x := factor(tissue_type, levels = names(lab), labels = lab)]
yr <- range(du$expr, na.rm = TRUE)
p_pm <- suppressWarnings(wilcox.test(du[tissue_type == "Primary", expr], du[tissue_type == "Metastasis", expr], exact = FALSE)$p.value)
annot <- if (p_pm < 0.001) "p<0.001" else paste0("p=", formatC(p_pm, digits = 3, format = "g"))
step <- 0.13 * diff(yr)
ypos <- yr[2] + step
p_un <- du %>%
  tidyplot(x = x, y = expr, color = tissue_type) %>%
  add_boxplot(alpha = 0, linewidth = 0.4) %>%
  adjust_colors(pal) %>%
  adjust_x_axis_title("") %>%
  adjust_y_axis_title("log2 normalized CDH3 (microarray)") %>%
  remove_legend() %>%
  add(theme_xrot) %>%
  add_data_points_jitter(jitter_width = 0.18, alpha = 0.85, size = 0.9) %>%
  add(geom_signif(
    comparisons = list(c(lab[["Primary"]], lab[["Metastasis"]])),
    annotations = annot, y_position = ypos, tip_length = 0.005,
    textsize = 2.8, vjust = -0.3, color = "black", size = 0.3, fontface = "plain"
  )) %>%
  add(expand_limits(y = max(ypos) + step * 1.5)) %>%
  add_caption("Wilcoxon rank-sum (unpaired)") %>%
  adjust_size(width = 50, height = 55)
save_plot(p_un, file.path(out_fhc, "r2q2_fhcrc_2col_overall_unpaired.pdf"), units = "mm", view_plot = FALSE, padding = 0.05, device = cairo_pdf)
save_plot(p_un, file.path(out_fhc, "r2q2_fhcrc_2col_overall_unpaired.png"), units = "mm", view_plot = FALSE, padding = 0.05, dpi = 600)

multi <- fhc[, uniqueN(tissue_type), by = patient][V1 >= 2L, patient]
dp <- fhc[patient %in% multi]
n <- dp[, .N, by = tissue_type][order(match(tissue_type, levels(dp$tissue_type)))]
lab <- setNames(paste0(n$tissue_type, "\n(n=", n$N, ")"), as.character(n$tissue_type))
dp[, x := factor(tissue_type, levels = names(lab), labels = lab)]
yr <- range(dp$expr, na.rm = TRUE)
a <- dp[tissue_type == "Primary", .(y = mean(expr, na.rm = TRUE)), by = patient]
b <- dp[tissue_type == "Metastasis", .(y = mean(expr, na.rm = TRUE)), by = patient]
pid <- intersect(a$patient, b$patient)
p_pm <- suppressWarnings(wilcox.test(a[match(pid, patient)]$y, b[match(pid, patient)]$y, paired = TRUE, exact = FALSE)$p.value)
annot <- if (p_pm < 0.001) "p<0.001" else paste0("p=", formatC(p_pm, digits = 3, format = "g"))
step <- 0.13 * diff(yr)
ypos <- yr[2] + step
p_pr <- dp %>%
  tidyplot(x = x, y = expr, color = tissue_type) %>%
  add_boxplot(alpha = 0, linewidth = 0.4) %>%
  adjust_colors(pal) %>%
  adjust_x_axis_title("") %>%
  adjust_y_axis_title("log2 normalized CDH3 (microarray)") %>%
  remove_legend() %>%
  add(theme_xrot) %>%
  add(geom_line(aes(group = patient), color = "grey75", alpha = 0.5, linewidth = 0.3, inherit.aes = TRUE)) %>%
  add_data_points(alpha = 0.9, size = 0.9) %>%
  add(geom_signif(
    comparisons = list(c(lab[["Primary"]], lab[["Metastasis"]])),
    annotations = annot, y_position = ypos, tip_length = 0.005,
    textsize = 2.8, vjust = -0.3, color = "black", size = 0.3, fontface = "plain"
  )) %>%
  add(expand_limits(y = max(ypos) + step * 1.5)) %>%
  add_caption("Wilcoxon signed-rank (paired, on patient-mean)") %>%
  adjust_size(width = 50, height = 55)
save_plot(p_pr, file.path(out_fhc, "r2q2_fhcrc_2col_overall_paired.pdf"), units = "mm", view_plot = FALSE, padding = 0.05, device = cairo_pdf)
save_plot(p_pr, file.path(out_fhc, "r2q2_fhcrc_2col_overall_paired.png"), units = "mm", view_plot = FALSE, padding = 0.05, dpi = 600)

record_session(output_dir("11_primary_metastatic_comparisons"))
