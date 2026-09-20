# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  inputs = list(
    brca_calling = "data/brca/TCGA_PAM50_calling.csv",
    brca_expr = "data/brca/TCGA_rsem_zscore.txt"
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
# Topic: TCGA BRCA biomarker/PAM50 concordance and CDH3 expression; Figures 1F/S1F.
require_packages(c("ggplot2", "ggalluvial", "dplyr", "svglite"))
library(ggplot2)
layout.dir <- paste0(output_dir("04_brca_subtype_sankey"), "/")
calling <- read.csv(file = input_file("brca_calling"),
                      header = T, check.names = F,
                      stringsAsFactors = F)

call <- calling$FINAL_CALL

ER <- rep('',nrow(calling))
ER[which(calling$er_status_by_ihc == 'Positive')] <- 'ER+'

PR <- rep('',nrow(calling))
PR[which(calling$pr_status_by_ihc == 'Positive')] <- 'PR+'

HER2 <- rep('',nrow(calling))
HER2[which(calling$HER2.newly.derived == 'Positive')] <- 'HER2+'

Subtype <- paste0(ER,PR,HER2)

Subtype[which(calling$`Triple Negative Status` == 'Yes')] <- 'TNBC'
Subtype[which(Subtype == '')] <- 'Unknown'

table(Subtype)

calling$Subtype <- Subtype
calling.tumor <- calling[which(calling$`Tumor or Normal`=='Tumor' & calling$`PAM50 and Claudin-low (CLOW) Molecular Subtype`!='CLOW'),]
calling.tumor.final <- calling.tumor[which(calling.tumor$Subtype != 'Unknown'),]
table(calling$Subtype)

## make Sankey plot
library(ggalluvial)
library(dplyr)
# compare PAM50 and histology
df <- calling.tumor.final

# Select PAM50 calls from the molecular-subtype field.
require_columns(df, "PAM50 and Claudin-low (CLOW) Molecular Subtype", "BRCA calls")
df$PAM50 <- df$`PAM50 and Claudin-low (CLOW) Molecular Subtype`

# change HER2E into HER2
df$PAM50[which(df$PAM50=='HER2E')] <- 'HER2'
df$PAM50 <- factor(df$PAM50, levels = c("Basal", "HER2", "LumA", "LumB","normal-like"))

df$Subtype <- factor(df$Subtype, levels = c("TNBC","HER2+","ER+HER2+",
                                            "PR+HER2+","ER+PR+HER2+","ER+","PR+","ER+PR+"))

# Count transitions between biomarker-defined groups and PAM50 subtypes.
df_count <- df %>%
  count(Subtype, PAM50) %>%
  rename(Freq = n)

# Define PAM50 colors from the boxplot
pam50_colors <- c("Basal" = "#e7970c",  # golden yellow
                  "LumA" = "#62767f",   # grayish
                  "LumB" = "#97a8bc",
                  "HER2" = "#517b5f",
                  "normal-like" = "#555b53")   # bluish gray

p <- ggplot(df_count,
            aes(axis1 = Subtype, axis2 = PAM50, y = Freq)) +
  geom_alluvium(aes(fill = PAM50), width = 1/12) +
  geom_stratum(width = 1/12, fill = "gray80", color = "black") +
  scale_fill_manual(values = pam50_colors) +
  scale_x_discrete(limits = c("Histology Classification", "PAM50 Subtyping"),
                   expand = c(.1, .1)) +

  # Default stratum text (no manual hjust or position_nudge)
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4, hjust = 0) +
  # Clean theme
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 11),
    plot.title = element_blank(),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  ylab(NULL)

jpeg(filename = paste0(layout.dir,'TNBC_Subtype_PAM50_Sankey_reordered_new.jpg'), height = 5, width = 7.5, res = 300, units = 'in')
print(p)
dev.off()
ggsave(paste0(layout.dir,'TNBC_Histology_PAM50_Sankey_reordered_new.svg'), plot=p, width=5, height=5)

write.csv(df_count, file = paste0(layout.dir,'count_histology_pam50_new.csv'))

# CDH3 table export for Fig1F box plots and FigS1F (GraphPad Prism).
rna.rsem <- read.table(input_file("brca_expr"), sep = "\t", header = TRUE, check.names = FALSE)
# extract the CDH3 from rnaseq
cdh3 <- t(rna.rsem[which(rna.rsem$Hugo_Symbol == 'CDH3'),-c(1,2)])
patient <- rownames(cdh3)
patient <- sub("^(([^-]*-){2}[^-]*).*", "\\1", patient)
df$CLID.clean <- sub("-[^-]*$", "", df$CLID)

ind.match <- match(df$CLID.clean, patient)

df.cdh3 <- data.frame(CDH3 = cdh3[ind.match], subtype = df$Subtype,
                      PAM50 = df$PAM50,
                     patient = df$CLID.clean)
colnames(df.cdh3) <- c('CDH3','Subtype','PAM50','Patient')
df_duplicates <- df.cdh3[df.cdh3$Patient %in% df.cdh3$Patient[duplicated(df.cdh3$Patient)], ]
write.csv(df.cdh3, file = paste0(layout.dir,'subtype_CDH3_new.csv'))
# Use the first expression sample per matched patient; export duplicate/unmatched records.
write.csv(df_duplicates, file.path(layout.dir, "duplicate_patient_rows.csv"), row.names = FALSE)
write.csv(df.cdh3[is.na(df.cdh3$CDH3), ], file.path(layout.dir, "unmatched_expression_rows.csv"), row.names = FALSE)
export_prism(df.cdh3$CDH3, df.cdh3$PAM50, layout.dir, "Fig1F_CDH3_PAM50")
export_prism(df.cdh3$CDH3, df.cdh3$Subtype, layout.dir, "FigS1F_CDH3_biomarkers")
# Finish box plots/statistical annotations in GraphPad Prism using the figure legend.
record_session(layout.dir)
