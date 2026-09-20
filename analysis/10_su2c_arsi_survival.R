# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  inputs = list(
    su2c_clinical = "data/su2c/Abida2019_TableS1.xlsx",
    su2c_expr = "data/su2c/fpkm_polya_standardized.txt",
    su2c_sample = "data/su2c/data_clinical_sample.txt"
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
# Topic: SU2C time on first-line ARSI; Figure S1G.
library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(patchwork)
library(readxl)

out <- output_dir("10_su2c_arsi_survival")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pal <- c(`CDH3-low` = "#377EB8", `CDH3-high` = "#E41A1C")

expr_path <- input_file("su2c_expr")
samp_path <- input_file("su2c_sample")
abida_path <- input_file("su2c_clinical")

dt <- fread(expr_path, sep = "\t", header = TRUE, check.names = FALSE)
gcol <- intersect(c("Gene_symbol", "Hugo_Symbol", "gene_name", "GeneSymbol", "gene_symbol"), colnames(dt))[1]
scols <- setdiff(colnames(dt), c(gcol, "Entrez_Gene_Id", "Entrez_ID"))
cdh3 <- as.numeric(dt[which(dt[[gcol]] == "CDH3")[1], scols, with = FALSE])
names(cdh3) <- scols
cdh3 <- log2(pmax(cdh3, 0) + 1)

skip <- which(!startsWith(readLines(samp_path), "#"))[1] - 1L
smap <- fread(samp_path, sep = "\t", skip = skip, header = TRUE,
              na.strings = c("", "NA", "[Not Available]", "[Not Applicable]", "[Unknown]"))
smap[, SAMPLE_ID := as.character(SAMPLE_ID)]
smap[, PATIENT_ID := as.character(PATIENT_ID)]

abida <- as.data.table(read_excel(abida_path, sheet = 2))
setnames(abida, c("Subject ID", "Time on first-line ARSI", "Off ARSI?"),
         c("PATIENT_ID", "time_on_arsi", "off_arsi"), skip_absent = TRUE)
abida[, PATIENT_ID := as.character(PATIENT_ID)]

df <- data.table(SAMPLE_ID = names(cdh3), cdh3 = unname(cdh3))
df <- merge(df, smap[, .(SAMPLE_ID, PATIENT_ID)], by = "SAMPLE_ID")
df <- merge(df, abida[, .(PATIENT_ID, time = as.numeric(time_on_arsi),
                          event = as.integer(off_arsi))], by = "PATIENT_ID")
df <- unique(df, by = "PATIENT_ID")[!is.na(cdh3) & !is.na(time) & !is.na(event) & time > 0]

# Calculate lower and upper tertile thresholds.
q <- quantile(df$cdh3, c(1 / 3, 2 / 3), na.rm = TRUE)
df[, median := ifelse(cdh3 > median(cdh3, na.rm = TRUE), "high",
               ifelse(cdh3 < median(cdh3, na.rm = TRUE), "low", NA_character_))]
df[, tertile := ifelse(cdh3 >= q[[2L]], "high",
                ifelse(cdh3 <= q[[1L]], "low", NA_character_))]
fwrite(df[, .(sample_id = SAMPLE_ID, patient_id = PATIENT_ID, cdh3, time, event, median, tertile)],
       file.path(out, "SU2C_mCRPC_2019_TimeOnARSI_source.csv"))
saveRDS(df, file.path(out, "SU2C_mCRPC_2019_TimeOnARSI_source.rds"))

km_arsi <- function(df, grp, fname) {
  d <- df[!is.na(df[[grp]])]
  d2 <- data.frame(time = d$time, event = d$event, grp = factor(d[[grp]], levels = c("low", "high")))
  cx <- coxph(Surv(time, event) ~ grp, data = d2)
  s <- summary(cx)
  ph <- cox.zph(cx)
  sd <- survdiff(Surv(time, event) ~ grp, data = d2)
  lr_p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  fit <- survfit(Surv(time, event) ~ grp, data = d2)
  n_high <- sum(d2$grp == "high")
  n_low <- sum(d2$grp == "low")
  hr <- unname(s$conf.int[1, "exp(coef)"])
  ci_lo <- unname(s$conf.int[1, "lower .95"])
  ci_hi <- unname(s$conf.int[1, "upper .95"])
  cox_p <- unname(s$coefficients[1, "Pr(>|z|)"])
  raw_strata <- sub("^grp=", "", names(fit$strata))
  lab_strata <- ifelse(raw_strata == "low", "CDH3-low", "CDH3-high")
  ndf <- data.frame(time = fit$time, surv = fit$surv, strata = rep(lab_strata, fit$strata))
  origin <- data.frame(time = 0, surv = 1, strata = unique(ndf$strata))
  ndf <- rbind(origin, ndf)
  ndf$strata <- factor(ndf$strata, levels = c("CDH3-low", "CDH3-high"))
  ann <- sprintf(
    paste0("CDH3-low n = %d; CDH3-high n = %d\n",
           "HR = %.2f (95%% CI %.2f-%.2f)\n",
           "Cox p = %s; Log-rank p = %s"),
    n_low, n_high, hr, ci_lo, ci_hi,
    formatC(cox_p, digits = 3, format = "g"),
    formatC(lr_p, digits = 3, format = "g")
  )
  tmax <- max(d2$time, na.rm = TRUE)
  x_limits <- c(0, ceiling(tmax / 5) * 5)
  brks <- pretty(x_limits, n = 6)
  brks <- brks[brks >= x_limits[1] & brks <= x_limits[2]]
  n_at_low <- colSums(outer(d2$time[d2$grp == "low"], brks, ">="))
  n_at_high <- colSums(outer(d2$time[d2$grp == "high"], brks, ">="))
  rt <- rbind(
    data.frame(strata = factor("CDH3-low", levels = c("CDH3-high", "CDH3-low")), time = brks, n_at = n_at_low),
    data.frame(strata = factor("CDH3-high", levels = c("CDH3-high", "CDH3-low")), time = brks, n_at = n_at_high)
  )
  axis_color <- pal[as.character(levels(rt$strata))]
  km <- ggplot(ndf, aes(x = time, y = surv, color = strata)) +
    geom_step(linewidth = 0.6) +
    scale_color_manual(values = pal, name = NULL) +
    coord_cartesian(ylim = c(0, 1.05)) +
    scale_x_continuous(limits = x_limits, breaks = brks, expand = expansion(mult = c(0.06, 0.02))) +
    labs(x = "Time on first-line ARSI (months)", y = "Probability remaining on ARSI", subtitle = ann) +
    theme_classic(13) +
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 13, color = "black"),
          legend.position = "inside", legend.position.inside = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = "white", color = NA),
          legend.title = element_blank(),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 11),
          plot.subtitle = element_text(size = 11, color = "black", hjust = 1,
                                       lineheight = 1.0, margin = margin(b = 2)),
          plot.title.position = "plot",
          plot.margin = margin(2, 3, 1, 2))
  risk_tbl <- ggplot(rt, aes(x = time, y = strata)) +
    geom_text(aes(label = n_at), size = 4.25, color = "black", vjust = 0.5, family = "sans") +
    scale_x_continuous(limits = x_limits, breaks = brks, expand = expansion(mult = c(0.06, 0.02))) +
    labs(x = NULL, y = NULL, title = "Number at risk") +
    theme_classic(12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 11, color = axis_color, hjust = 1, face = "plain",
                                     margin = margin(r = 4)),
          plot.title = element_text(size = 12, color = "black", hjust = 0, margin = margin(b = 2, l = 18)),
          plot.title.position = "plot",
          plot.margin = margin(1, 2, 1, 2))
  p <- free(km, type = "label", side = "l") / risk_tbl + plot_layout(heights = c(3.6, 1))
  p[[1]] <- p[[1]] +
    theme(axis.title.y.left = element_text(hjust = 0.5, vjust = 0.5, margin = margin(r = 6)))
  p_titled <- p + plot_annotation(
    title = "SU2C mCRPC 2019",
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )
  ggsave(paste0(fname, "_clean.pdf"), p, width = 4.5, height = 4.0, device = cairo_pdf)
  ggsave(paste0(fname, "_clean.png"), p, width = 4.5, height = 4.0, dpi = 600, bg = "white")
  ggsave(paste0(fname, "_titled.pdf"), p_titled, width = 4.32, height = 3.84, device = cairo_pdf)
  ggsave(paste0(fname, "_titled.png"), p_titled, width = 4.32, height = 3.84, dpi = 600, bg = "white")
  file.copy(paste0(fname, "_clean.pdf"), paste0(fname, ".pdf"), overwrite = TRUE)
  file.copy(paste0(fname, "_clean.png"), paste0(fname, ".png"), overwrite = TRUE)
  saveRDS(p, paste0(fname, "_plot.rds"))
  data.table(
    cohort = "SU2C_mCRPC_2019", endpoint = "TimeOnARSI", stratification = grp,
    n_low = n_low, n_high = n_high,
    events_low = sum(d2$event[d2$grp == "low"]),
    events_high = sum(d2$event[d2$grp == "high"]),
    HR_high_vs_low = hr, CI_low = ci_lo, CI_high = ci_hi,
    cox_wald_p = cox_p, cox_ph_assumption_p = unname(ph$table["grp", "p"]),
    log_rank_p = lr_p, source_data = "SU2C_mCRPC_2019_TimeOnARSI_source.csv"
  )
}

stats_med <- km_arsi(df, "median", file.path(out, "SU2C_mCRPC_2019_TimeOnARSI_median_KM"))
stats_ter <- km_arsi(df, "tertile", file.path(out, "SU2C_mCRPC_2019_TimeOnARSI_tertile_KM"))
fwrite(rbind(stats_med, stats_ter), file.path(out, "SU2C_mCRPC_2019_TimeOnARSI_stats.csv"))

record_session(out)
