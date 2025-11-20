.libPaths("/home/heshidian/R/x86_64-pc-linux-gnu-library/4.3")
library(Seurat)
library(ggplot2)
library(ggforce)
library(scales)
library(cowplot)
library(grid)
library(gridExtra)
library(ggh4x)
library(egg)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(dplyr)
library(SeuratDisk)
library(spacexr)

seurat_out <- readRDS("./pythonProject/seurat_out_annotated_MAGIC.rds")

Celltype_color <- c("#A4C9DD", "#D5231E", "#399938", "#EF7C1C",
                    "#F19695", "#C6B0D2", "#F5BB6F", "#483D8B",
                    "#ADD487", "#2572A9", "#1E90FF")
names(Celltype_color) <- levels(seurat_out$Celltype)

HSD_spatial_dimplot_highlight_2 <- function(seurat_obj = seurat_out,
                                            group_color = Celltype_color,
                                            group_feature = "Celltype",
                                            highlight_member = c("Proliferative_progenitor cells",
                                                                 "Mural cells",
                                                                 "Endothelial cells"),
                                            highlight_order = NULL,
                                            panel_width = 13) {
  xlab <- "Spatial 1"
  ylab <- "Spatial 2"
  embeddings <- Embeddings(seurat_obj, reduction = "spatial")
  embeddings <- as.data.frame(embeddings)
  x_length <- range(embeddings[,1])[2] - range(embeddings[,1])[1]
  y_length <- range(embeddings[,2])[2] - range(embeddings[,2])[1]
  y_x_ratio <- y_length / x_length
  features <- group_feature
  if (features %in% colnames(seurat_obj@meta.data)) {
    embeddings$features <- seurat_obj@meta.data[,features]
  } else {
    stop()
  }
  if (!all(highlight_member %in% embeddings$features)) {
    stop("The highlighted group members are not completely exsited.")
  }
  embeddings$features <- as.character(embeddings$features)
  embeddings$features <- ifelse(embeddings$features %in% highlight_member,
                                embeddings$features, "Others")
  if (is.null(highlight_order)) {
    highlight_order <- levels(seurat_obj@meta.data[,features])
    highlight_order <- highlight_order[highlight_order %in% highlight_member]
  }
  group_color <- group_color[highlight_order]
  group_color <- c(group_color, "gray")
  names(group_color)[length(group_color)] <- "Others"
  highlight_order <- c(highlight_order, "Others")
  embeddings$features <- factor(as.character(embeddings$features),
                                levels = highlight_order)
  n <- nrow(embeddings)
  p <- ggplot(data = embeddings,
              aes(x = embeddings[,1],
                  y = embeddings[,2],
                  color = features)) +
    geom_point(size = 0.000001) +
    scale_color_manual(values = group_color) +
    theme_bw() +
    ggtitle("") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"), # 设置绘图区域背景色
      axis.text = element_blank(),
      axis.title = element_text(size = 10, colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 12, colour = "black", face = "bold"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 10, colour = "black")
    ) + 
    labs(x = xlab, y = ylab, color = features) +
    guides(color = guide_legend(override.aes = list(size = 1.5)))
  p <- set_panel_size(p,
                      width  = unit(panel_width, "cm"),
                      height = unit(panel_width * y_x_ratio, "cm"))
  return(p)
}

p <- HSD_spatial_dimplot_highlight_2()
pdf("./pythonProject/15.Spatial_dimplot_highlight_selected_celltypes.pdf", height = 10, width = 10)
grid.draw(p)
dev.off()


p <- HSD_spatial_dimplot_highlight_2(seurat_obj = seurat_out,
                                     group_color = Celltype_color,
                                     group_feature = "Celltype",
                                     highlight_member = c("AnMCs",
                                                          "Progenitor cells",
                                                          "Proliferative_progenitor cells",
                                                          "Chondrocytes",
                                                          "Hypertrophic chondrocytes")
)
pdf("./highlight_selected_celltypes.pdf", height = 10, width = 10)
grid.draw(p)
dev.off()



HSD_spatial_featureplot_gene <- function(seurat_obj,
                                         features,
                                         panel_width = 13, title = NULL,
                                         assay = "SCT",
                                         max_cut = 0.9,
                                         pal_gradiant = NULL) {
  if (!assay %in% names(seurat_obj@assays)) {
    stop("The specified assay is not found in the Seurat object.")
  }
  xlab <- "Spatial 1"
  ylab <- "Spatial 2"
  embeddings <- Embeddings(seurat_obj, reduction = "spatial")
  embeddings <- as.data.frame(embeddings)
  x_length <- range(embeddings[,1])[2] - range(embeddings[,1])[1]
  y_length <- range(embeddings[,2])[2] - range(embeddings[,2])[1]
  y_x_ratio <- y_length / x_length
  # features <- features
  # if (is.null(pt.size)) {
  #   pt.size <- (4.5 / x_length) * 10
  # }
  if (features %in% colnames(seurat_obj@meta.data)) {
    embeddings$features <- seurat_obj@meta.data[,features]
  } else if (features %in% rownames(seurat_obj@assays[[assay]]@data)) {
    embeddings$features <- seurat_obj@assays[[assay]]@data[features,]
  } else {
    stop("The specified feature is not found in either meta.data or in the rownames of the object.")
  }
  n <- nrow(embeddings)
  if (is.null(pal_gradiant)) {
    colours <- c("#10357A", "#1072A7", "#50B777",
                 "#EDAC37", "#E2712E", "#FF0000")
  } else {
    colours <- pal_gradiant
  }
  max_value <- max(embeddings$features)
  embeddings$features[embeddings$features > max_value*max_cut] <- max_value*max_cut
  p <- ggplot() +
    geom_point(data = embeddings, aes(x = embeddings[,1],
                                      y = embeddings[,2],
                                      color = embeddings[,"features"]),
               size = 0.001) +
    scale_color_gradientn(colours = colours) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "white"), # 设置绘图区域背景色
      # axis.text = element_text(size = 10, colour = "black"),
      axis.text = element_blank(),
      axis.title = element_text(size = 10, colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 12, colour = "black", face = "bold"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 10, colour = "black")
      # aspect.ratio = 230/130
      # plot.background = element_rect(fill = "#10357A", color = NA)    # 设置整个图的背景色
    ) + 
    labs(x = xlab, y = ylab, color = features)
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      labs(color = "Value") + 
      theme(plot.title = element_text(hjust = 0.5, size = 12, colour = "black", face = "bold"))
  }
  p <- set_panel_size(p,
                      width  = unit(panel_width, "cm"),
                      height = unit(panel_width * y_x_ratio, "cm"))
  return(p)
}

seurat_out_selected <- subset(seurat_out,
                              subset = Celltype %in% c("AnMCs",
                                                       "Progenitor cells",
                                                       "Proliferative_progenitor cells",
                                                       "Chondrocytes",
                                                       "Hypertrophic chondrocytes"))
p <- HSD_spatial_featureplot_gene(seurat_obj = seurat_out_selected,
                                  features = "HIF1A",
                                  panel_width = 13,
                                  title = "HIF1A",
                                  max_cut = 0.8,
                                  assay = "MAGIC_Spatial")
pdf("HIF1A.pdf", width = 10, height = 12)
grid.newpage()
grid.draw(p)
dev.off()


p <- HSD_spatial_featureplot_gene(seurat_obj = seurat_out_selected,
                                  features = "VEGFA",
                                  panel_width = 13,
                                  title = "VEGFA",
                                  max_cut = 0.8,
                                  assay = "MAGIC_Spatial")
dev.off()
pdf("VEGFA.pdf", width = 10, height = 12)
grid.newpage()
grid.draw(p)
dev.off()



p <- HSD_spatial_featureplot_gene(seurat_obj = seurat_out_selected,
                                  features = "ANGPT1",
                                  panel_width = 13,
                                  title = "ANGPT1",
                                  max_cut = 0.8,
                                  assay = "MAGIC_Spatial")
dev.off()
pdf("ANGPT1.pdf", width = 10, height = 12)
grid.newpage()
grid.draw(p)
dev.off()


p <- HSD_spatial_featureplot_gene(seurat_obj = seurat_out_selected,
                                  features = "ANGPT2",
                                  panel_width = 13,
                                  title = "ANGPT2",
                                  max_cut = 0.5,
                                  assay = "MAGIC_Spatial")
dev.off()
pdf("ANGPT2.pdf", width = 10, height = 12)
grid.newpage()
grid.draw(p)
dev.off()



p <- HSD_spatial_featureplot_gene(seurat_obj = seurat_out_selected,
                                  features = "PDGFB",
                                  panel_width = 13,
                                  title = "PDGFB",
                                  max_cut = 0.3,
                                  assay = "MAGIC_Spatial")
dev.off()
pdf("PDGFB.pdf", width = 10, height = 12)
grid.newpage()
grid.draw(p)
dev.off()


seurat_out_selected <- subset(seurat_out,
                              subset = Celltype %in% c("AnMCs",
                                                       "Proliferative_progenitor cells",
                                                       "Progenitor cells"))

pdf("PTH1R.pdf", width = 10, height = 12)
p <- HSD_spatial_featureplot_gene(seurat_obj = seurat_out_selected,
                                  features = "PTH1R",
                                  panel_width = 13,
                                  title = "PTH1R",
                                  max_cut = 0.9,
                                  assay = "MAGIC_Spatial")
grid.newpage()
grid.draw(p)
dev.off()

