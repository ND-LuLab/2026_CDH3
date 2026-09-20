# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  inputs = list(
    wong_all = "data/single_cell/wong.rds"
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

save_figure <- function(plot, out, stem, width = 5, height = 4) {
  ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggplot2::ggsave(file.path(out, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, bg = "white")
}

# Analysis and visualization ---------------------------------------------------
# Topic: Wong single-cell analysis and visualization. FigS2B-E.
# Start from a processed Seurat object. Dataset-specific cluster IDs and UMAP
# coordinates are essential: these manual gates must not be applied to reclustered data.
require_packages(c("Seurat", "ggplot2", "RColorBrewer", "patchwork"))
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(grid)
dataset <- "wong"
out.dir <- paste0(output_dir("06_wong_prostate_single_cell"), "/")
ds <- readRDS(input_file("wong_all"))
require_columns(ds@meta.data, c("seurat_clusters", "predicted_label"), "Seurat metadata")
if (!"umap" %in% names(ds@reductions)) stop("Stored UMAP coordinates are required")
DefaultAssay(ds) <- "RNA"
# 1. Annotate cell populations and select the epithelial compartment.
ds@meta.data$celltype <- as.character(ds@meta.data$predicted_label)
Idents(ds) <- ds@meta.data$seurat_clusters
ds <- subset(ds, idents = setdiff(levels(Idents(ds)), c(26,25,21,19,23,20)))
ds@meta.data$celltype[which(ds@meta.data$celltype == "B cell")] <- 'Plasma cell'
ds@meta.data$celltype[which(ds@meta.data$celltype == "Cytotoxic T cell")] <- 'T cell'
ds@meta.data$celltype[which(ds@meta.data$celltype == "Monocyte")] <- 'Neutrophil'
ds@meta.data$celltype[which(ds@meta.data$celltype == "Macrophage")] <- 'Macrophage/Monocyte'

p_colors <- colorRampPalette(brewer.pal(12, "Paired"))(10)
ds@meta.data$celltype <- factor(ds@meta.data$celltype, levels = c("Basal epithelial cell", "Luminal epithelial cell",
                                                                 "T cell","Plasma cell","Macrophage/Monocyte","Neutrophil","Natural killer cell",
                                                                 "Fibroblast","Endothelial cell","Mast cell"))

saveRDS(ds,paste0(out.dir,dataset,"_reannotated.rds"))

my_colors <- p_colors
my_colors[7] <- p_colors[1]
my_colors[c(1,2)] <- c('orange2','skyblue3')

p <- DimPlot(ds, group.by = "celltype", cols = my_colors) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = "") & NoAxes()

pdf(paste0(out.dir,dataset,"_annotation_label.pdf"), width = 6, height = 4)
print(LabelClusters(p, id = "celltype", color = 'black', size = 3, repel = T,  box.padding = 1))
dev.off()

pdf(paste0(out.dir,dataset,"_annotation.pdf"), width = 5.5, height = 3.5)
print(DimPlot(ds, group.by = "celltype", label = F, cols = my_colors) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = "") & NoAxes())
dev.off()

umap_coords <- Embeddings(ds, "umap")
ds@meta.data$UMAP_1 <- umap_coords[, 1]
ds@meta.data$UMAP_2 <- umap_coords[, 2]

ds.epi <- subset(ds,subset = seurat_clusters %in% c(1,6,9,10,13,15,16) &
                  UMAP_1 > -5 & UMAP_2 > 0 & UMAP_2 < 12 &
                  UMAP_2 > 4 - UMAP_1)

celltype <- rep('Luminal',nrow(ds.epi@meta.data))
celltype[which(ds.epi$seurat_clusters %in% c(6,13,16))] <- "Basal"
ds.epi@meta.data$celltype <- celltype

p <- DimPlot(ds.epi,group.by = "celltype", cols = c('orange2','skyblue3')) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = "") & NoAxes() +
  theme(legend.direction = "horizontal", legend.position = 'bottom', legend.justification = "center",
        legend.text = element_text(size = 20), legend.key.size = unit(2, 'lines'))
pdf(paste0(out.dir,dataset,"_epi.pdf"), width = 4, height = 4)
print(p)
dev.off()

saveRDS(ds.epi,paste0(out.dir,dataset,"_epi.rds"))

gene <- "CDH3"
p <- FeaturePlot(ds.epi, features = gene, order = T, pt.size = 0.5) +
  scale_color_gradientn(colors = c("grey90", "orange", "red", "darkred"),limits=c(0,3.5)) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = paste0(gene,"\n")) & NoAxes() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box = "horizontal", legend.key.width = unit(0.5, 'cm'))
pdf(paste0(out.dir,dataset,"_epi_",gene,".pdf"), width = 4, height = 4)
print(p)
dev.off()

dp <- DotPlot(ds.epi,features = c("CDH3"), group.by = "celltype") +
  scale_color_gradient2(low = "#5ba2cb",   # Color for low values
                        mid = "white",   # Color for mid values
                        high = "red3",    # Color for high values
                        midpoint = 0     # Midpoint of the color scale
  ) +
  scale_size_continuous(range = c(2, 10)) +
  theme_classic() +
  theme(axis.line = element_blank(),        # Remove axis lines
        axis.ticks = element_blank(),       # Remove axis ticks
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title = element_blank()) +
  guides(
    color = guide_colorbar(title = "Avg Exp", title.theme = element_text(size = 14)),  # Customize color legend
    size = guide_legend(title = "CDH3\n\nPercent Exp", title.theme = element_text(size = 14))  # Customize size legend
  )
pdf(paste0(out.dir,dataset,"_epi_CDH3_dot.pdf"), width = 3, height = 4)
print(dp)
dev.off()

# 2. Export annotated cell identities and plotting values.
write.csv(ds@meta.data, file.path(out.dir, "all_cell_annotations.csv"))
write.csv(ds.epi@meta.data, file.path(out.dir, "epithelial_annotations.csv"))

# Visualize WNT5A and ROR2 in the epithelial compartment.
for (gene in c("WNT5A", "ROR2")) {
  if (!gene %in% rownames(ds.epi)) stop("Missing requested gene: ", gene)
  p <- FeaturePlot(ds.epi, features = gene, order = TRUE, pt.size = 0.5) +
    scale_color_gradientn(colors = c("grey90", "orange", "red", "darkred"), limits = c(0, 3.5)) +
    labs(x = "UMAP1", y = "UMAP2", title = "", color = paste0(gene, "\n")) & NoAxes() +
    theme(legend.position = "bottom", legend.justification = "center",
          legend.box = "horizontal", legend.key.width = unit(0.5, "cm"))
  save_figure(p, out.dir, paste0(dataset, "_epi_", gene), 4, 4)
}

# 3. Plot canonical and noncanonical WNT expression.
require_packages(c("dplyr", "ggh4x"))
library(dplyr)
library(ggh4x)
  # Define gene sets
  gene_set_b <- c("WNT3A","FZD5","LRP5","LRP6","AXIN2","GSK3B","TCF7")
  gene_set_a <- c("WNT5A","FZD7","ROR2","DKK1","PRICKLE1","NFATC2","FOSL1")

  # Create gene-to-set mapping
  gene_map <- rbind(
    data.frame(gene = gene_set_a, set = "Noncanonical WNT"),
    data.frame(gene = gene_set_b, set = "Canonical WNT")
  )

  # Get DotPlot data
  p <- DotPlot(ds.epi, features = c(gene_set_a, gene_set_b), group.by = "celltype",
               scale = FALSE)

  # Merge and plot
  p <- p$data %>%
    left_join(gene_map, by = c("features.plot" = "gene")) %>%
    mutate(id = factor(id, levels = rev(levels(id))),
           set = factor(set, levels = c("Noncanonical WNT","Canonical WNT")),
           features.plot=factor(features.plot,levels=c(gene_set_a,gene_set_b))) %>%
    ggplot(aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    facet_wrap2(~set, scales = "free_x", nrow = 1,
                strip = strip_themed(
                  background_x = elem_list_rect(fill = c('orange2','skyblue3')),
                  text_x = elem_list_text(color = "white")
                )) +
    scale_color_gradient2(low = "white", high = "firebrick3",
                          name = "Avg Exp") +
    scale_size(name = "% Exp") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank())
  pdf(file = paste0(out.dir,dataset,"_WNT_dot_plot_bar.pdf"),width = 6, height = 2)
  print(p)
  dev.off()

# Export dot-plot values and basal/luminal expression comparisons.
dot_data <- DotPlot(ds.epi, features = c(gene_set_a, gene_set_b), group.by = "celltype", scale = FALSE)$data
write.csv(dot_data, file.path(out.dir, "WNT_dotplot_values.csv"), row.names = FALSE)
Idents(ds.epi) <- "celltype"
compare <- FindMarkers(ds.epi, ident.1 = "Basal", ident.2 = "Luminal")
genelist <- c("CDH3", gene_set_a, gene_set_b)
avg <- AverageExpression(ds.epi, features = genelist, group.by = "celltype")$RNA
ds.epi$all_cells <- "all"
avg_all <- AverageExpression(ds.epi, features = genelist, group.by = "all_cells")$RNA
values <- data.frame(avg, avg_all)
values <- cbind(values, compare[rownames(values), c("p_val", "p_val_adj")])
write.csv(values[genelist, ], file.path(out.dir, "selected_basal_luminal_expression.csv"))

record_session(out.dir)
