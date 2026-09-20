# Local settings ---------------------------------------------------------------
# Run from your repository root. Edit these input paths for your environment.
# Inputs and supporting functions are defined in this file.
config <- list(
  output_root = "results",
  inputs = list(
    xu_all = "data/single_cell/Xu.rds"
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

# Analysis and visualization ---------------------------------------------------
# Topic: Xu single-cell analysis and visualization. FigS2I-J.
# Start from a processed Seurat object. Dataset-specific cluster IDs and UMAP
# coordinates are essential: these manual gates must not be applied to reclustered data.
require_packages(c("Seurat", "ggplot2", "RColorBrewer", "patchwork"))
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(grid)
dataset <- "Xu"
out.dir <- paste0(output_dir("09_xu_breast_single_cell"), "/")
ds <- readRDS(input_file("xu_all"))
require_columns(ds@meta.data, c("seurat_clusters", "predicted_label"), "Seurat metadata")
if (!"umap" %in% names(ds@reductions)) stop("Stored UMAP coordinates are required")
DefaultAssay(ds) <- "RNA"
# 1. Annotate cell populations and select the epithelial compartment.
ds@meta.data$celltype <- as.character(ds@meta.data$predicted_label)
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters %in% c(16,2,11,10))] <- 'Basal epithelial cells'
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters %in% c(4,6))] <- 'Luminal epithelial cells'
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters == 12)] <- 'B cells'
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters == 5)] <- 'NK cells'
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters == 13)] <- 'pDCs'
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters == 15)] <- 'cDCs'
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters %in% c(17,1))] <- 'T cells'
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters == 7)] <- 'Neutrophils'
ds@meta.data$celltype[which(ds@meta.data$seurat_clusters == 9)] <- 'Macrophages/Monocytes'
Idents(ds) <- ds@meta.data$seurat_clusters

p_colors <- colorRampPalette(brewer.pal(12, "Paired"))(15)
ds@meta.data$celltype <- factor(ds@meta.data$celltype, levels = c("Basal epithelial cells", "Luminal epithelial cells",
                                                                 "T cells","NK cells","B cells","Macrophages/Monocytes","Neutrophils","pDCs","cDCs",
                                                                 "Fibroblasts","Endothelial cells","Mast cells","Pericytes","Proliferating cells"))
my_colors <- p_colors
my_colors[10] <- p_colors[1]
my_colors[c(1,2)] <- c('orange2','skyblue3')

p <- DimPlot(ds, group.by = "celltype", cols = my_colors) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = "") & NoAxes()

pdf(paste0(out.dir,dataset,"_annotation_label.pdf"), width = 6, height = 4)
print(LabelClusters(p, id = "celltype", color = 'black', size = 3, repel = T,  box.padding = 1))
dev.off()

pdf(paste0(out.dir,dataset,"_annotation.pdf"), width = 6.5, height = 3.5)
print(DimPlot(ds, group.by = "celltype", label = F, cols = my_colors) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = "") & NoAxes())
dev.off()

umap_coords <- Embeddings(ds, "umap")
ds@meta.data$UMAP_1 <- umap_coords[, 1]
ds@meta.data$UMAP_2 <- umap_coords[, 2]

ds.epi <- subset(ds, subset = seurat_clusters %in% c(2,10,11,16,4,6) & UMAP_1 < 0 & UMAP_2 > -7)

celltype <- rep('Luminal',nrow(ds.epi@meta.data))
celltype[which(ds.epi$seurat_clusters %in% c(16,2,11,10))] <- "Basal"
ds.epi@meta.data$celltype <- celltype

p <- DimPlot(ds.epi,group.by = "celltype", cols = c('orange2','skyblue3')) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = "") & NoAxes() +
  theme(legend.direction = "horizontal", legend.position = 'bottom', legend.justification = "center",
        legend.text = element_text(size = 20), legend.key.size = unit(2, 'lines'))
pdf(paste0(out.dir,dataset,"_epi.pdf"), width = 4, height = 4)
print(p)
dev.off()

saveRDS(ds.epi,paste0(out.dir,dataset,"_epi.rds"))

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
    size = guide_legend(title = "Percent Exp", title.theme = element_text(size = 14))  # Customize size legend
  )
pdf(paste0(out.dir,dataset,"_epi_CDH3_dot.pdf"), width = 3, height = 4)
print(dp)
dev.off()

gene <- "CDH3"
p <- FeaturePlot(ds, features = gene, order = T, pt.size = 0.5) +
  scale_color_gradientn(colors = c("grey90", "orange", "red", "darkred"),limits=c(0,3.5)) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = paste0(gene,"\n")) & NoAxes() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box = "horizontal", legend.key.width = unit(0.5, 'cm'))
pdf(paste0(out.dir,dataset,"_",gene,".pdf"), width = 4, height = 4)
print(p)
dev.off()

gene <- "CDH3"
p <- FeaturePlot(ds.epi, features = gene, order = T, pt.size = 0.5) +
  scale_color_gradientn(colors = c("grey90", "orange", "red", "darkred"),limits=c(0,3.5)) +
  labs(x = "UMAP1",y = "UMAP2",title = "",color = paste0(gene,"\n")) & NoAxes() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box = "horizontal", legend.key.width = unit(0.5, 'cm'))
pdf(paste0(out.dir,dataset,"_epi_",gene,".pdf"), width = 4, height = 4)
print(p)
dev.off()

# 2. Export annotated cell identities and plotting values.
write.csv(ds@meta.data, file.path(out.dir, "all_cell_annotations.csv"))
write.csv(ds.epi@meta.data, file.path(out.dir, "epithelial_annotations.csv"))

record_session(out.dir)
