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

HSD_spatial_featureplot_gene <- function(seurat_obj,
                                         features,
                                         panel_width = 13, title = NULL,
                                         assay = "SCT",
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
  p <- ggplot() +
    geom_point(data = embeddings, aes(x = embeddings[,1],
                                      y = embeddings[,2],
                                      color = embeddings[,"features"]),
               size = 0.001) +
    scale_color_gradientn(colours = colours) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "#10357A", color = "black"), # 设置绘图区域背景色
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
HSD_spatial_dimplot <- function(seurat_obj,
                                features, group_color = NULL,
                                group_order = NULL,
                                panel_width = 13) {
  xlab <- "Spatial 1"
  ylab <- "Spatial 2"
  embeddings <- Embeddings(seurat_obj, reduction = "spatial")
  embeddings <- as.data.frame(embeddings)
  x_length <- range(embeddings[,1])[2] - range(embeddings[,1])[1]
  y_length <- range(embeddings[,2])[2] - range(embeddings[,2])[1]
  y_x_ratio <- y_length / x_length
  features <- features
  if (features %in% colnames(seurat_obj@meta.data)) {
    embeddings$features <- seurat_obj@meta.data[,features]
  } else {
    stop()
  }
  group_size <- length(unique(as.character(embeddings$features)))
  if (is.null(group_color)) {
    group_color <- scales::hue_pal()(group_size)
  }
  if (!is.null(group_order)) {
    embeddings$features <- factor(as.character(embeddings$features),
                                  levels = group_order)
  }
  n <- nrow(embeddings)
  p <- ggplot(data = embeddings, aes(x = embeddings[,1], y = embeddings[,2], color = features)) +
    geom_point(size = 0.000001) +
    scale_color_manual(values = group_color) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "#D3D3D3", color = "black"), # 设置绘图区域背景色
      # axis.text = element_text(size = 10, colour = "black"),
      axis.text = element_blank(),
      axis.title = element_text(size = 10, colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 12, colour = "black"),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 10, colour = "black")
      # aspect.ratio = 230/130
      # plot.background = element_rect(fill = "#10357A", color = NA)    # 设置整个图的背景色
    ) + 
    labs(x = xlab, y = ylab, color = features) +
    guides(color = guide_legend(override.aes = list(size = 1.5)))
  p <- set_panel_size(p,
                      width  = unit(panel_width, "cm"),
                      height = unit(panel_width * y_x_ratio, "cm"))
  return(p)
}
HSD_spatial_dimplot_highlight <- function(seurat_obj,
                                          features,
                                          group_order = NULL,
                                          panel_width = 13) {
  xlab <- "Spatial 1"
  ylab <- "Spatial 2"
  embeddings <- Embeddings(seurat_obj, reduction = "spatial")
  embeddings <- as.data.frame(embeddings)
  x_length <- range(embeddings[,1])[2] - range(embeddings[,1])[1]
  y_length <- range(embeddings[,2])[2] - range(embeddings[,2])[1]
  y_x_ratio <- y_length / x_length
  features <- features
  if (features %in% colnames(seurat_obj@meta.data)) {
    embeddings$features <- seurat_obj@meta.data[,features]
  } else {
    stop()
  }
  group_size <- length(unique(as.character(embeddings$features)))
  if (!is.null(group_order)) {
    embeddings$features <- factor(as.character(embeddings$features),
                                  levels = group_order)
  }
  n <- nrow(embeddings)
  p <- c()
  for (x in as.character(levels(embeddings$features))) {
    temp <- embeddings
    temp$features <- as.character(temp$features)
    temp$features <- ifelse(temp$features == x, as.character(x), "Others")
    temp$features <- factor(temp$features,
                            levels = c(x, "Others"))
    p1 <- ggplot(data = temp,
                 aes(x = temp[,1],
                     y = temp[,2],
                     color = features)) +
      geom_point(size = 0.000001) +
      scale_color_manual(values = c("red", "gray")) +
      theme_bw() +
      ggtitle(paste0(features, ": ", x)) +
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
    p1 <- set_panel_size(p1,
                         width  = unit(panel_width, "cm"),
                         height = unit(panel_width * y_x_ratio, "cm"))
    p <- c(p, list(p1))
  }
  
  return(p)
}
HSD_ggplot2_theme_1 <- function(p, title = "") {
  p <- p + theme_bw() +
    ggtitle(title) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 12, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10, colour = "black"))
  p
}


pdf("14.Celltype_marker_featureplot(AnMCs)_spatial.pdf",
    width = 10, height = 10)
for (i in c("SFRP2", "IGF1", "PRRX1", "THY1", "TWIST2")) {
  p <- HSD_spatial_featureplot_gene(seurat_obj = seurat_out,
                                    features = i,
                                    panel_width = 13,
                                    title = i,
                                    assay = "SCT",
                                    pal_gradiant = c("gray", "blue", "red"))
  grid.newpage()
  grid.draw(p)
}
dev.off()


