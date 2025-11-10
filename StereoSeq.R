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

HSD_spatial_featureplot <- function(seurat_obj,
                                features, limits,
                                panel_width = 13, title = NULL) {
  xlab <- "Spatial 1"
  ylab <- "Spatial 2"
  embeddings <- Embeddings(seurat_obj, reduction = "spatial")
  embeddings <- as.data.frame(embeddings)
  x_length <- range(embeddings[,1])[2] - range(embeddings[,1])[1]
  y_length <- range(embeddings[,2])[2] - range(embeddings[,2])[1]
  y_x_ratio <- y_length / x_length
  features <- features
  # if (is.null(pt.size)) {
  #   pt.size <- (4.5 / x_length) * 10
  # }
  if (features %in% colnames(seurat_obj@meta.data)) {
    embeddings$features <- seurat_obj@meta.data[,features]
  } else if (features %in% rownames(seurat_obj)) {
    embeddings$features <- seurat_obj@assays$Spatial@data[features,]
  } else {
    stop()
  }
  n <- nrow(embeddings)
  p <- ggplot() +
    geom_point(data = embeddings, aes(x = embeddings[,1],
                                      y = embeddings[,2],
                                      color = features),
               size = 0.001) +
    scale_color_gradientn(colours = c("#10357A", "#1072A7", "#50B777",
                                     "#EDAC37", "#E2712E", "#FF0000"),
                         limits = c(limits[1], limits[2]),
                         oob = squish,
                         labels = label_number(scale_cut = cut_short_scale()) # 将数值格式化为 K（千）、M（百万）等
    ) +
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
HSD_spatial_embed_reflect <- function(seurat_obj) {
  spatial_embeddings <- seurat_obj@reductions[["spatial"]]@cell.embeddings
  mid <- median(spatial_embeddings[,2])
  spatial_embeddings <- split.data.frame(spatial_embeddings,
                                         f = list(spatial_embeddings[,1]))
  spatial_embeddings <- lapply(spatial_embeddings,
                               function(x){
                                 x[,2] <- 2*mid - x[,2]
                                 as.data.frame(x)
                               })
  spatial_embeddings <- data.table::rbindlist(spatial_embeddings)
  spatial_embeddings <- as.data.frame(spatial_embeddings)
  rownames(spatial_embeddings) <- rownames(seurat_obj@reductions[["spatial"]]@cell.embeddings)
  spatial_embeddings <- as.matrix(spatial_embeddings)
  seurat_obj@reductions[["spatial"]]@cell.embeddings <- as.matrix(spatial_embeddings)
  return(seurat_obj)
}

setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/Stereo-seq/pythonProject")
seurat_out <- readRDS("seurat_out_Remove_fur.rds")
image <- seurat_out@images
seurat_out@images <- list()
seurat_out <- HSD_spatial_embed_reflect(seurat_out)
seurat_out <- SCTransform(seurat_out, assay = "Spatial",
                          verbose = TRUE)
seurat_out <- RunPCA(seurat_out, assay = "SCT", verbose = TRUE)
seurat_out <- FindNeighbors(seurat_out, reduction = "pca", dims = 1:30)
seurat_out <- FindClusters(seurat_out, resolution = 0.935)
seurat_out <- RunUMAP(seurat_out, reduction = "pca", dims = 1:30,
                      min.dist = 0.01)

p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                             limits = c(0, 45000),
                             features = "nCount_Spatial")
dev.off()
pdf("1.nCount_Spatial.pdf",
    width = 10, height = 10)
grid.newpage()
grid.draw(p)
dev.off()

set.seed(125)
cluster_colors <- randomcoloR::distinctColorPalette(40)

p <- DimPlot(seurat_out, group.by = "seurat_clusters",
        cols = cluster_colors, label = T)
p <- HSD_ggplot2_theme_1(p)
pdf("2.clusters_umap.pdf",
    width = 5.5, height = 4.7)
print(p)
dev.off()

p <- HSD_spatial_dimplot(seurat_out,
                         features = "seurat_clusters",
                         group_color = cluster_colors)
dev.off()
pdf("2.clusters_spatial.pdf",
    width = 10, height = 10)
grid.newpage()
grid.draw(p)
dev.off()

table(seurat_out$seurat_clusters)

p <- HSD_spatial_dimplot_highlight(seurat_obj = seurat_out,
                                   features = "seurat_clusters",
                                   panel_width = 13)
pdf("3.clusters_spatial_splited.pdf",
    height = 10)
for (i in 1:length(p)) {
  grid.newpage()
  grid.draw(p[[i]])
}
dev.off()

p <- DimPlot(seurat_out, group.by = "seurat_clusters",
             cols = cluster_colors, label = T,
             split.by = "seurat_clusters", ncol = 4)
p <- HSD_ggplot2_theme_1(p)
pdf("3.clusters_umap_splited.pdf",
    width = 16, height = 20)
print(p)
dev.off()


# SPOTlight, Spatial
{
  # 提取表达矩阵
  counts <- GetAssayData(seurat_out, slot = "counts", assay = "Spatial")
  # 提取空间坐标
  spatial_coords <- seurat_out@reductions$spatial
  # 提取元数据
  metadata <- seurat_out@meta.data
  spe_obj <- SpatialExperiment(
    assays = list(counts = counts),
    colData = metadata,
    spatialCoords = as.matrix(spatial_coords@cell.embeddings),
    imgData = NULL
  )
  
  sce <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/RNA_1.RNA_singlet.rds")
  sce <- as.SingleCellExperiment(sce)
  sce <- logNormCounts(sce)
  genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(sce))
  dec <- modelGeneVar(sce, subset.row = genes)
  plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
  curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
  hvg <- getTopHVGs(dec, n = 3000)
  colLabels(sce) <- colData(sce)$Celltype_rename
  # Compute marker genes
  mgs <- scoreMarkers(sce, subset.row = genes)
  mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.7
    x <- x[x$mean.AUC > 0.7, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)
  table(mgs_df$cluster)
  
  # split cell indices by identity
  idx <- split(seq(ncol(sce)), sce$Celltype_rename)
  n_cells <- 100
  cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
      n_cells <- n
    sample(i, n_cells)
  })
  sce <- sce[, unlist(cs_keep)]
  table(sce$Celltype_rename)
  res <- SPOTlight(
    x = sce,
    y = seurat_out,
    groups = as.character(sce$Celltype_rename),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")
  mod <- res$NMF
  plotTopicProfiles(
    x = mod,
    y = sce$Celltype_rename,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
  plotTopicProfiles(
    x = mod,
    y = sce$Celltype_rename,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6)
  library(NMF)
  sign <- basis(mod)
  colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
  
  mat <- res[["mat"]]
  mat[is.nan(mat)] <- 0
  plotInteractions(mat,
                   which = "heatmap",
                   metric = "jaccard")
  plotInteractions(mat, which = "network")
  saveRDS(mat, "SPOTlight_Spatial_result.rds")
  ct <- colnames(mat)
  # mat[mat < 0.1] <- 0
  paletteMartin <- c(
    "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
    "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
    "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
  
  pal <- colorRampPalette(paletteMartin)(length(ct))
  names(pal) <- ct
  
  seurat_out$SPOTlight_spatial_ss <- res[[2]][colnames(spe_obj)]
  quantile(seurat_out$SPOTlight_spatial_ss, na.rm = TRUE)
  
  p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                               limits = c(0, 0.1),
                               features = "SPOTlight_spatial_ss")
  dev.off()
  pdf("SPOTlight_Spatial_residuals.pdf",
      width = 10, height = 10)
  grid.newpage()
  grid.draw(p)
  dev.off()
  
  # mat_2 <- mat[rowSums(mat) > 0,]
  # SPOTlight_MoransI <- RunMoransI(data = t(mat_2),
  #                                 pos = seurat_out@reductions[["spatial"]]@cell.embeddings[colnames(t(mat_2)),])
  SPOTlight_Spatial_result <- readRDS("SPOTlight_Spatial_result.rds")
  colnames(SPOTlight_Spatial_result) <- paste0("SPOTlight_Spatial: ",
                                               colnames(SPOTlight_Spatial_result))
  seurat_out@meta.data <- as.data.frame(cbind(seurat_out@meta.data,
                                              SPOTlight_Spatial_result[colnames(seurat_out),]))
  pdf("4.SPOTlight_spatial.pdf", height = 12)
  for (i in colnames(SPOTlight_Spatial_result)) {
    p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                                 limits = c(0, max(seurat_out@meta.data[,i])),
                                 features = i, title = i)
    grid.newpage()
    grid.draw(p)
  }
  dev.off()
}

# SPOTlight, SCT
{
  # 提取表达矩阵
  counts <- GetAssayData(seurat_out, slot = "counts", assay = "SCT")
  # 提取空间坐标
  spatial_coords <- seurat_out@reductions$spatial
  # 提取元数据
  metadata <- seurat_out@meta.data
  spe_obj <- SpatialExperiment(
    assays = list(counts = counts),
    colData = metadata,
    spatialCoords = as.matrix(spatial_coords@cell.embeddings),
    imgData = NULL
  )
  
  sce <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/RNA_1.RNA_singlet.rds")
  sce <- as.SingleCellExperiment(sce)
  sce <- logNormCounts(sce)
  genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(sce))
  dec <- modelGeneVar(sce, subset.row = genes)
  plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
  curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
  hvg <- getTopHVGs(dec, n = 3000)
  colLabels(sce) <- colData(sce)$Celltype_rename
  # Compute marker genes
  mgs <- scoreMarkers(sce, subset.row = genes)
  mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.7
    x <- x[x$mean.AUC > 0.7, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)
  table(mgs_df$cluster)
  
  # split cell indices by identity
  idx <- split(seq(ncol(sce)), sce$Celltype_rename)
  n_cells <- 100
  cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
      n_cells <- n
    sample(i, n_cells)
  })
  sce <- sce[, unlist(cs_keep)]
  table(sce$Celltype_rename)
  res <- SPOTlight(
    x = sce,
    y = seurat_out,
    groups = as.character(sce$Celltype_rename),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")
  mod <- res$NMF
  plotTopicProfiles(
    x = mod,
    y = sce$Celltype_rename,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
  plotTopicProfiles(
    x = mod,
    y = sce$Celltype_rename,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6)
  library(NMF)
  sign <- basis(mod)
  colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
  
  mat <- res[["mat"]]
  mat[is.nan(mat)] <- 0
  plotInteractions(mat,
                   which = "heatmap",
                   metric = "jaccard")
  plotInteractions(mat, which = "network")
  saveRDS(mat, "SPOTlight_SCT_result.rds")
  ct <- colnames(mat)
  # mat[mat < 0.1] <- 0
  paletteMartin <- c(
    "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
    "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
    "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
  
  pal <- colorRampPalette(paletteMartin)(length(ct))
  names(pal) <- ct
  
  seurat_out$SPOTlight_SCT_ss <- res[[2]][colnames(spe_obj)]
  quantile(seurat_out$SPOTlight_SCT_ss, na.rm = TRUE)
  
  p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                               limits = c(0, 0.1),
                               features = "SPOTlight_SCT_ss")
  dev.off()
  pdf("SPOTlight_SCT_residuals.pdf",
      width = 10, height = 10)
  grid.newpage()
  grid.draw(p)
  dev.off()
  
  # mat_2 <- mat[rowSums(mat) > 0,]
  # SPOTlight_MoransI <- RunMoransI(data = t(mat_2),
  #                                 pos = seurat_out@reductions[["spatial"]]@cell.embeddings[colnames(t(mat_2)),])
  SPOTlight_SCT_result <- readRDS("SPOTlight_SCT_result.rds")
  colnames(SPOTlight_SCT_result) <- paste0("SPOTlight_SCT: ",
                                               colnames(SPOTlight_SCT_result))
  seurat_out@meta.data <- as.data.frame(cbind(seurat_out@meta.data,
                                              SPOTlight_SCT_result[colnames(seurat_out),]))
  pdf("4.SPOTlight_SCT.pdf", height = 12)
  for (i in colnames(SPOTlight_SCT_result)) {
    p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                                 limits = c(0, max(seurat_out@meta.data[,i])),
                                 features = i, title = i)
    grid.newpage()
    grid.draw(p)
  }
  dev.off()
}

# Seurat
{
  sce <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/RNA_1.RNA_singlet.rds")
  sce.list <- SplitObject(sce, split.by = "sample")
  sce.list <- lapply(sce.list, function(x){
    x <- SCTransform(x, vst.flavor = "v2", verbose = TRUE #,
                     # vars.to.regress = "percent.mt"
                     )
    x <- RunPCA(x, npcs = 30, verbose = TRUE)
    x
  })
  features <- SelectIntegrationFeatures(object.list = sce.list,
                                        nfeatures = 3000)
  sce.list <- PrepSCTIntegration(object.list = sce.list,
                                 anchor.features = features)
  sce.list.anchors <- FindIntegrationAnchors(object.list = sce.list,
                                             normalization.method = "SCT",
                                             reduction = c("rpca"),
                                             anchor.features = features)
  combined.sce <- IntegrateData(anchorset = sce.list.anchors,
                                normalization.method = "SCT")
  DefaultAssay(combined.sce)
  combined.sce <- RunPCA(combined.sce, verbose = TRUE)
  combined.sce <- RunUMAP(combined.sce, reduction = "pca",
                          dims = 1:30, verbose = TRUE)
  DimPlot(combined.sce, group.by = "sample")
  DimPlot(combined.sce, label = T)
  
  
  anchors <- FindTransferAnchors(reference = combined.sce, query = seurat_out,
                                 normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, refdata = combined.sce$Celltype_rename,
                                    prediction.assay = TRUE,
                                    weight.reduction = seurat_out[["pca"]], dims = 1:30)
  seurat_out[["predictions"]] <- predictions.assay
  DefaultAssay(seurat_out) <- "predictions"
  # seurat_MoransI <- RunMoransI(data = seurat_out@assays[["predictions"]]@data,
  #                       pos = seurat_out@reductions[["spatial"]]@cell.embeddings)
  
  Seurat_result <- as.data.frame(t(seurat_out@assays[["predictions"]]@data))
  colnames(Seurat_result) <- paste0("Seurat: ",
                                    colnames(Seurat_result))
  seurat_out@meta.data <- as.data.frame(cbind(seurat_out@meta.data,
                                              Seurat_result[colnames(seurat_out),]))
  pdf("4.Seurat_spatial.pdf", height = 12)
  for (i in colnames(Seurat_result)) {
    p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                                 limits = c(0, max(seurat_out@meta.data[,i])),
                                 features = i, title = i)
    grid.newpage()
    grid.draw(p)
  }
  dev.off()
  
}

# RCTD, SCT
{
  # 创建 SpatialRNA 对象
  spatialRNA <- SpatialRNA(coords = as.data.frame(seurat_out@reductions[["spatial"]]@cell.embeddings),
                           counts = seurat_out@assays[["SCT"]]@counts)
  
  # 创建 Reference 对象
  reference <- Reference(counts = combined.sce@assays[["SCT"]]@counts,
                         cell_types = as.factor(combined.sce$Celltype_rename))
  # 创建 RCTD 对象，设置最大核心数以提高计算效率（根据您的计算资源）
  myRCTD <- create.RCTD(spatialRNA, reference, max_cores = 20)
  # 运行 RCTD 分析
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD, "myRCTD_SCTCount_full.rds")
  # 获取结果
  results <- myRCTD@results
  saveRDS(results, "RCTD_SCTCount_full_result.rds")
  results <- readRDS("RCTD_SCTCount_full_result.rds")
  # 提取权重矩阵（每个空间点的细胞类型比例）
  weights <- results$weights
  # RCTD_MoransI <- RunMoransI(data = t(weights),
  #                            pos = seurat_out@reductions[["spatial"]]@cell.embeddings)
  myRCTD_SCTCount_full <- readRDS("myRCTD_SCTCount_full.rds")
  myRCTD_SCTCount_full <- myRCTD_SCTCount_full@results
  myRCTD_SCTCount_full <- myRCTD_SCTCount_full[["weights"]]
  myRCTD_SCTCount_full <- as.data.frame(myRCTD_SCTCount_full)
  colnames(myRCTD_SCTCount_full) <- paste0("RCTD_SCTCount: ",
                                           colnames(myRCTD_SCTCount_full))
  seurat_out@meta.data <- as.data.frame(cbind(seurat_out@meta.data,
                                              myRCTD_SCTCount_full[colnames(seurat_out),]))
  pdf("4.RCTD_SCTCount_spatial.pdf", height = 12)
  for (i in colnames(myRCTD_SCTCount_full)) {
    p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                                 limits = c(0, max(seurat_out@meta.data[,i])),
                                 features = i, title = i)
    grid.newpage()
    grid.draw(p)
  }
  dev.off()
}

# RCTD, Spatial
{
  # 创建 SpatialRNA 对象
  spatialRNA <- SpatialRNA(coords = as.data.frame(seurat_out@reductions[["spatial"]]@cell.embeddings),
                           counts = seurat_out@assays[["Spatial"]]@counts)
  
  # 创建 Reference 对象
  reference <- Reference(counts = combined.sce@assays[["SCT"]]@counts,
                         cell_types = as.factor(combined.sce$Celltype_rename))
  # 创建 RCTD 对象，设置最大核心数以提高计算效率（根据您的计算资源）
  myRCTD <- create.RCTD(spatialRNA, reference, max_cores = 20)
  # 运行 RCTD 分析
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD, "myRCTD_SpatialCount_full.rds")
  # 获取结果
  results <- myRCTD@results
  saveRDS(results, "RCTD_SpatialCount_full_result.rds")
  results <- readRDS("RCTD_SpatialCount_full_result.rds")
  # 提取权重矩阵（每个空间点的细胞类型比例）
  weights <- results$weights
  # RCTD_MoransI <- RunMoransI(data = t(weights),
  #                            pos = seurat_out@reductions[["spatial"]]@cell.embeddings)
  myRCTD_SpatialCount_full <- readRDS("myRCTD_SpatialCount_full.rds")
  myRCTD_SpatialCount_full <- myRCTD_SpatialCount_full@results
  myRCTD_SpatialCount_full <- myRCTD_SpatialCount_full[["weights"]]
  myRCTD_SpatialCount_full <- as.data.frame(myRCTD_SpatialCount_full)
  colnames(myRCTD_SpatialCount_full) <- paste0("RCTD_SpatialCount: ",
                                           colnames(myRCTD_SpatialCount_full))
  myRCTD_SpatialCount_full <- myRCTD_SpatialCount_full[colnames(seurat_out),]
  seurat_out@meta.data <- as.data.frame(cbind(seurat_out@meta.data,
                                              myRCTD_SpatialCount_full[colnames(seurat_out),]))
  pdf("4.RCTD_SpatialCount_spatial.pdf", height = 12)
  for (i in colnames(myRCTD_SpatialCount_full)) {
    p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                                 limits = c(0, max(seurat_out@meta.data[,i])),
                                 features = i, title = i)
    grid.newpage()
    grid.draw(p)
  }
  dev.off()
}

# cell2location
{
  seurat_out_2 <- seurat_out
  DefaultAssay(seurat_out_2) <- "Spatial"
  seurat_out_2_C0_6 <- subset(seurat_out_2, subset = seurat_clusters %in% c("0", "1", "2", "3",
                                                                            "4", "5", "6"))
  seurat_out_2_C0_6@assays <- seurat_out_2_C0_6@assays["Spatial"]
  SaveH5Seurat(seurat_out_2_C0_6, filename = "seurat_out_2_C0_6.h5Seurat",
               overwrite = TRUE)
  Convert("seurat_out_2_C0_6.h5Seurat", dest = "h5ad", overwrite = TRUE)
  
  seurat_out_2_C7_16 <- subset(seurat_out_2, subset = seurat_clusters %in% c("7", "8", "9", "10", "11",
                                                                             "12", "13", "14", "15", "16"))
  seurat_out_2_C7_16@assays <- seurat_out_2_C7_16@assays["Spatial"]
  SaveH5Seurat(seurat_out_2_C7_16, filename = "seurat_out_2_C7_16.h5Seurat",
               overwrite = TRUE)
  Convert("seurat_out_2_C7_16.h5Seurat", dest = "h5ad", overwrite = TRUE)
  
  # 保存 Seurat 对象为 h5Seurat 格式
  sce <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/RNA_1.RNA_singlet.rds")
  sce <- CreateSeuratObject(counts = sce@assays$RNA@counts, meta.data = sce@meta.data)
  DefaultAssay(sce) <- "RNA"
  SaveH5Seurat(sce, filename = "scRNA-seq.h5Seurat",
               overwrite = TRUE)
  # 将 h5Seurat 文件转换为 h5ad 格式
  Convert("scRNA-seq.h5Seurat", dest = "h5ad",
          overwrite = TRUE)
  
  Celltype_map <- data.frame(number = 0:9,
                             label = levels(sce$Celltype_rename))
  rownames(Celltype_map) <- as.character(Celltype_map$number)
  st_cell2location_res_0_6 <- read.csv("st_cell2location_res_0-6.csv", row.names = 1)
  colnames(st_cell2location_res_0_6) <- paste0("Cell2location: ", levels(sce$Celltype_rename))
  st_cell2location_res_7_16 <- read.csv("st_cell2location_res_7-16.csv", row.names = 1)
  colnames(st_cell2location_res_7_16) <- paste0("Cell2location: ", levels(sce$Celltype_rename))
  
  st_cell2location_res <- as.data.frame(rbind(st_cell2location_res_0_6,
                                              st_cell2location_res_7_16))
  
  seurat_out@meta.data <- as.data.frame(cbind(seurat_out@meta.data,
                                              st_cell2location_res[colnames(seurat_out),]))
  pdf("4.Cell2location_spatial.pdf", height = 12)
  for (i in colnames(st_cell2location_res)) {
    p <- HSD_spatial_featureplot(seurat_obj = seurat_out,
                                 limits = c(0, max(seurat_out@meta.data[,i])),
                                 features = i, title = i)
    grid.newpage()
    grid.draw(p)
  }
  dev.off()
}


SPOTlight_Spatial_result <- readRDS("SPOTlight_Spatial_result.rds")
pdf("5.SPOTlight_Spatial_dotplot.pdf",
    height = 10, width = 6)
DotPlot(seurat_out, features = paste0("SPOTlight_Spatial: ", colnames(SPOTlight_Spatial_result))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("white", "#FFF0F5", "#F5DEB3", "#FF8C00", "red"))
dev.off()
SPOTlight_SCT_result <- readRDS("SPOTlight_SCT_result.rds")
pdf("5.SPOTlight_SCT_dotplot.pdf",
    height = 10, width = 6)
DotPlot(seurat_out, features = paste0("SPOTlight_SCT: ", colnames(SPOTlight_SCT_result))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("white", "#FFF0F5", "#F5DEB3", "#FF8C00", "red"))
dev.off()

pdf("5.RCTD_SpatialCount_dotplot.pdf",
    height = 10, width = 6)
DotPlot(seurat_out, features = colnames(myRCTD_SpatialCount_full)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("white", "#FFF0F5", "#F5DEB3", "#FF8C00", "red"))
dev.off()

pdf("5.RCTD_SCTCount_dotplot.pdf",
    height = 10, width = 6)
DotPlot(seurat_out, features = colnames(myRCTD_SCTCount_full)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("white", "#FFF0F5", "#F5DEB3", "#FF8C00", "red"))
dev.off()

pdf("5.Seurat_dotplot.pdf",
    height = 10, width = 6)
DotPlot(seurat_out, features = colnames(Seurat_result)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("white", "#FFF0F5", "#F5DEB3", "#FF8C00", "red"))
dev.off()

pdf("5.Cell2location_dotplot.pdf",
    height = 10, width = 6)
DotPlot(seurat_out, features = paste0("Cell2location: ", levels(sce$Celltype_rename))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c("white", "#FFF0F5", "#F5DEB3", "#FF8C00", "red"))
dev.off()


plist <- lapply(paste0("Cell2location: ", levels(sce$Celltype_rename)), function(x){
  FeaturePlot(seurat_out, features = x, order = T) +
    scale_color_gradientn(colours = c("#10357A", "#1072A7", "#50B777",
                                      "#EDAC37", "#E2712E", "#FF0000"))
})
p <- cowplot::plot_grid(plotlist = plist, nrow = 3)
pdf("6.Cell2location_umap_featurePlot.pdf",
    height = 14, width = 20)
p
dev.off()

plist <- lapply(colnames(Seurat_result), function(x){
  FeaturePlot(seurat_out, features = x, order = T) +
    scale_color_gradientn(colours = c("#10357A", "#1072A7", "#50B777",
                                      "#EDAC37", "#E2712E", "#FF0000"))
})
p <- cowplot::plot_grid(plotlist = plist, nrow = 3)
pdf("6.Seurat_umap_featurePlot.pdf",
    height = 14, width = 20)
p
dev.off()


plist <- lapply(colnames(myRCTD_SCTCount_full), function(x){
  FeaturePlot(seurat_out, features = x, order = T) +
    scale_color_gradientn(colours = c("#10357A", "#1072A7", "#50B777",
                                      "#EDAC37", "#E2712E", "#FF0000"))
})
p <- cowplot::plot_grid(plotlist = plist, nrow = 3)
pdf("6.RCTD_SCTCount_umap_featurePlot.pdf",
    height = 14, width = 20)
p
dev.off()


plist <- lapply(colnames(myRCTD_SpatialCount_full), function(x){
  FeaturePlot(seurat_out, features = x, order = T) +
    scale_color_gradientn(colours = c("#10357A", "#1072A7", "#50B777",
                                      "#EDAC37", "#E2712E", "#FF0000"))
})
p <- cowplot::plot_grid(plotlist = plist, nrow = 3)
pdf("6.RCTD_SpatialCount_umap_featurePlot.pdf",
    height = 14, width = 20)
p
dev.off()

SPOTlight_SCT_result
plist <- lapply(paste0("SPOTlight_SCT: ", colnames(SPOTlight_SCT_result)),
                function(x){
  FeaturePlot(seurat_out, features = x, order = T) +
    scale_color_gradientn(colours = c("#10357A", "#1072A7", "#50B777",
                                      "#EDAC37", "#E2712E", "#FF0000"))
})
p <- cowplot::plot_grid(plotlist = plist, nrow = 3)
pdf("6.SPOTlight_SCT_umap_featurePlot.pdf",
    height = 14, width = 20)
p
dev.off()

SPOTlight_Spatial_result
plist <- lapply(paste0("SPOTlight_Spatial: ", colnames(SPOTlight_Spatial_result)),
                function(x){
                  FeaturePlot(seurat_out, features = x, order = T) +
                    scale_color_gradientn(colours = c("#10357A", "#1072A7", "#50B777",
                                                      "#EDAC37", "#E2712E", "#FF0000"))
                })
p <- cowplot::plot_grid(plotlist = plist, nrow = 3)
pdf("6.SPOTlight_Spatial_umap_featurePlot.pdf",
    height = 14, width = 20)
p
dev.off()


# SPOTlight_Spatial
{
  SPOTlight_Spatial <- seurat_out@meta.data[,paste0("SPOTlight_Spatial: ",
                                                    levels(sce$Celltype_rename))]
  SPOTlight_Spatial <- colnames(SPOTlight_Spatial)[unlist(apply(SPOTlight_Spatial, 1, FUN = which.max))]
  seurat_out$SPOTlight_Spatial <- SPOTlight_Spatial
  Celltype_color <-  c("#004949", "#009292", "#ff6db6", "#ffb6db", 
                       "#006ddb", "#b66dff", "#6db6ff", "#924900",
                       "#db6d00", "#ffff6d")
  names(Celltype_color) <- paste0("SPOTlight_Spatial: ",
                                  levels(sce$Celltype_rename))
  p <- HSD_spatial_dimplot_highlight(seurat_obj = seurat_out,
                                     features = "SPOTlight_Spatial",
                                     panel_width = 13,
                                     group_order = names(Celltype_color))
  pdf("7.SPOTlight_Spatial_CelltypeAnno_splited.pdf",
      width = 9.5, height = 12)
  for (i in 1:length(p)) {
    grid.newpage()
    grid.draw(p[[i]])
  }
  dev.off()
  
  p <- HSD_spatial_dimplot(seurat_out,
                           features = "SPOTlight_Spatial",
                           group_color = Celltype_color)
  dev.off()
  pdf("8.SPOTlight_Spatial_CelltypeAnno_spatial.pdf",
      width = 10, height = 10)
  grid.newpage()
  grid.draw(p)
  dev.off()
  p <- DimPlot(seurat_out, group.by = "SPOTlight_Spatial",
               cols = Celltype_color, label = F)
  p <- HSD_ggplot2_theme_1(p)
  pdf("8.SPOTlight_Spatial_CelltypeAnno_umap.pdf",
      width = 8, height = 4.7)
  print(p)
  dev.off()
  
}

# SPOTlight_SCT
{
  SPOTlight_SCT <- seurat_out@meta.data[,paste0("SPOTlight_SCT: ",
                                                levels(sce$Celltype_rename))]
  SPOTlight_SCT <- colnames(SPOTlight_SCT)[unlist(apply(SPOTlight_SCT, 1, FUN = which.max))]
  seurat_out$SPOTlight_SCT <- SPOTlight_SCT
  Celltype_color <-  c("#004949", "#009292", "#ff6db6", "#ffb6db", 
                       "#006ddb", "#b66dff", "#6db6ff", "#924900",
                       "#db6d00", "#ffff6d")
  names(Celltype_color) <- paste0("SPOTlight_SCT: ",
                                  levels(sce$Celltype_rename))
  p <- HSD_spatial_dimplot_highlight(seurat_obj = seurat_out,
                                     features = "SPOTlight_SCT",
                                     panel_width = 13,
                                     group_order = names(Celltype_color))
  pdf("7.SPOTlight_SCT_CelltypeAnno_splited.pdf",
      width = 9.5, height = 12)
  for (i in 1:length(p)) {
    grid.newpage()
    grid.draw(p[[i]])
  }
  dev.off()
  
  p <- HSD_spatial_dimplot(seurat_out,
                           features = "SPOTlight_SCT",
                           group_color = Celltype_color)
  dev.off()
  pdf("8.SPOTlight_SCT_CelltypeAnno_spatial.pdf",
      width = 10, height = 10)
  grid.newpage()
  grid.draw(p)
  dev.off()
  p <- DimPlot(seurat_out, group.by = "SPOTlight_SCT",
               cols = Celltype_color, label = F)
  p <- HSD_ggplot2_theme_1(p)
  pdf("8.SPOTlight_SCT_CelltypeAnno_umap.pdf",
      width = 8, height = 4.7)
  print(p)
  dev.off()
  
}

# RCTD_Spatial
{
  RCTD_Spatial <- seurat_out@meta.data[,paste0("RCTD_SpatialCount: ",
                                                levels(sce$Celltype_rename))]
  RCTD_Spatial_2 <- c()
  for (i in 1:nrow(RCTD_Spatial)) {
    if (sum(is.na(RCTD_Spatial[i,])) > 0) {
      RCTD_Spatial_2 <- c(RCTD_Spatial_2, "Unidentified")
    } else {
      RCTD_Spatial_2 <- c(RCTD_Spatial_2,
                          colnames(RCTD_Spatial)[which.max(RCTD_Spatial[i,])])
    }
  }
  seurat_out$RCTD_Spatial <- RCTD_Spatial_2
  Celltype_color <-  c("#004949", "#009292", "#ff6db6", "#ffb6db", 
                       "#006ddb", "#b66dff", "#6db6ff", "#924900",
                       "#db6d00", "#ffff6d", "gray")
  names(Celltype_color) <- c(paste0("RCTD_SpatialCount: ",
                                    levels(sce$Celltype_rename)), "Unidentified")
  p <- HSD_spatial_dimplot_highlight(seurat_obj = seurat_out,
                                     features = "RCTD_Spatial",
                                     panel_width = 13,
                                     group_order = names(Celltype_color))
  dev.off()
  pdf("7.RCTD_Spatial_CelltypeAnno_splited.pdf",
      width = 9.5, height = 12)
  for (i in 1:length(p)) {
    grid.newpage()
    grid.draw(p[[i]])
  }
  dev.off()
  
  p <- HSD_spatial_dimplot(seurat_out,
                           features = "RCTD_Spatial",
                           group_color = Celltype_color)
  dev.off()
  pdf("8.RCTD_Spatial_CelltypeAnno_spatial.pdf",
      width = 10, height = 10)
  grid.newpage()
  grid.draw(p)
  dev.off()
  p <- DimPlot(seurat_out, group.by = "RCTD_Spatial",
               cols = Celltype_color, label = F)
  p <- HSD_ggplot2_theme_1(p)
  pdf("8.RCTD_Spatial_CelltypeAnno_umap.pdf",
      width = 8, height = 4.7)
  print(p)
  dev.off()
  
}

# RCTD_SCT
{
  RCTD_SCT <- seurat_out@meta.data[,paste0("RCTD_SCTCount: ",
                                           levels(sce$Celltype_rename))]
  RCTD_SCT_2 <- c()
  for (i in 1:nrow(RCTD_SCT)) {
    if (sum(is.na(RCTD_SCT[i,])) > 0) {
      RCTD_SCT_2 <- c(RCTD_SCT_2, "Unidentified")
    } else {
      RCTD_SCT_2 <- c(RCTD_SCT_2,
                      colnames(RCTD_SCT)[which.max(RCTD_SCT[i,])])
    }
  }
  seurat_out$RCTD_SCT <- RCTD_SCT_2
  Celltype_color <-  c("#004949", "#009292", "#ff6db6", "#ffb6db", 
                       "#006ddb", "#b66dff", "#6db6ff", "#924900",
                       "#db6d00", "#ffff6d", "gray")
  names(Celltype_color) <- c(paste0("RCTD_SCTCount: ",
                                    levels(sce$Celltype_rename)), "Unidentified")
  p <- HSD_spatial_dimplot_highlight(seurat_obj = seurat_out,
                                     features = "RCTD_SCT",
                                     panel_width = 13,
                                     group_order = names(Celltype_color))
  dev.off()
  pdf("7.RCTD_SCT_CelltypeAnno_splited.pdf",
      width = 9.5, height = 12)
  for (i in 1:length(p)) {
    grid.newpage()
    grid.draw(p[[i]])
  }
  dev.off()
  
  p <- HSD_spatial_dimplot(seurat_out,
                           features = "RCTD_SCT",
                           group_color = Celltype_color)
  dev.off()
  pdf("8.RCTD_SCT_CelltypeAnno_spatial.pdf",
      width = 10, height = 10)
  grid.newpage()
  grid.draw(p)
  dev.off()
  p <- DimPlot(seurat_out, group.by = "RCTD_SCT",
               cols = Celltype_color, label = F)
  p <- HSD_ggplot2_theme_1(p)
  pdf("8.RCTD_SCT_CelltypeAnno_umap.pdf",
      width = 8, height = 4.7)
  print(p)
  dev.off()
  
}

# Seurat
{
  Seurat_label <- paste0("Seurat: ",
                         levels(sce$Celltype_rename))
  Seurat_label <- gsub(pattern = "_", replacement = "-", fixed = T, x = Seurat_label)
  Seurat <- seurat_out@meta.data[,Seurat_label]
  Seurat_2 <- c()
  for (i in 1:nrow(Seurat)) {
    if (sum(is.na(Seurat[i,])) > 0) {
      Seurat_2 <- c(Seurat_2, "Unidentified")
    } else {
      Seurat_2 <- c(Seurat_2,
                    colnames(Seurat)[which.max(Seurat[i,])])
    }
  }
  seurat_out$Seurat <- Seurat_2
  Celltype_color <-  c("#004949", "#009292", "#ff6db6", "#ffb6db", 
                       "#006ddb", "#b66dff", "#6db6ff", "#924900",
                       "#db6d00", "#ffff6d", "gray")
  names(Celltype_color) <- c(Seurat_label, "Unidentified")
  p <- HSD_spatial_dimplot_highlight(seurat_obj = seurat_out,
                                     features = "Seurat",
                                     panel_width = 13,
                                     group_order = names(Celltype_color))
  dev.off()
  pdf("7.Seurat_CelltypeAnno_splited.pdf",
      width = 9.5, height = 12)
  for (i in 1:length(p)) {
    grid.newpage()
    grid.draw(p[[i]])
  }
  dev.off()
  
  p <- HSD_spatial_dimplot(seurat_out,
                           features = "Seurat",
                           group_color = Celltype_color)
  dev.off()
  pdf("8.Seurat_CelltypeAnno_spatial.pdf",
      width = 10, height = 10)
  grid.newpage()
  grid.draw(p)
  dev.off()
  p <- DimPlot(seurat_out, group.by = "Seurat",
               cols = Celltype_color, label = F)
  p <- HSD_ggplot2_theme_1(p)
  pdf("8.Seurat_CelltypeAnno_umap.pdf",
      width = 8, height = 4.7)
  print(p)
  dev.off()
  
}

# Cell2location
{
  Cell2location <- seurat_out@meta.data[,paste0("Cell2location: ",
                                                levels(sce$Celltype_rename))]
  Cell2location_2 <- c()
  for (i in 1:nrow(Cell2location)) {
    if (sum(is.na(Cell2location[i,])) > 0) {
      Cell2location_2 <- c(Cell2location_2, "Unidentified")
    } else {
      Cell2location_2 <- c(Cell2location_2,
                    colnames(Cell2location)[which.max(Cell2location[i,])])
    }
  }
  seurat_out$Cell2location <- Cell2location_2
  Celltype_color <-  c("#004949", "#009292", "#ff6db6", "#ffb6db", 
                       "#006ddb", "#b66dff", "#6db6ff", "#924900",
                       "#db6d00", "#ffff6d", "gray")
  names(Celltype_color) <- c(paste0("Cell2location: ",
                                    levels(sce$Celltype_rename)), "Unidentified")
  p <- HSD_spatial_dimplot_highlight(seurat_obj = seurat_out,
                                     features = "Cell2location",
                                     panel_width = 13,
                                     group_order = names(Celltype_color))
  dev.off()
  pdf("7.Cell2location_CelltypeAnno_splited.pdf",
      width = 9.5, height = 12)
  for (i in 1:length(p)) {
    grid.newpage()
    grid.draw(p[[i]])
  }
  dev.off()
  
  p <- HSD_spatial_dimplot(seurat_out,
                           features = "Cell2location",
                           group_color = Celltype_color)
  dev.off()
  pdf("8.Cell2location_CelltypeAnno_spatial.pdf",
      width = 10, height = 10)
  grid.newpage()
  grid.draw(p)
  dev.off()
  p <- DimPlot(seurat_out, group.by = "Cell2location",
               cols = Celltype_color, label = F)
  p <- HSD_ggplot2_theme_1(p)
  pdf("8.Cell2location_CelltypeAnno_umap.pdf",
      width = 8, height = 4.7)
  print(p)
  dev.off()
  
}

DefaultAssay(seurat_out) <- "SCT"
seurat_out@assays$SCT@data["TNC",]

genes <- c("TNC", "TNN", "RUNX2", "ACTA2", "MYLK",
           "MYL9", "CDH5", "MKI67", "TOP2A")
genes %in% rownames(seurat_out@assays$SCT@data)
plist <- lapply(genes, function(x){
  FeaturePlot(seurat_out, features = c(x), order = T,
              reduction = "umap", min.cutoff = 0) +
    scale_color_gradientn(colours = c("#10357A", "#1072A7", "#50B777",
                                      "#EDAC37", "#E2712E", "#FF0000"))
})
pdf("9.Progenitor_markers_featurePlot_umap.pdf",
    height = 10, width = 11)
cowplot::plot_grid(plotlist = plist, nrow = 3)
dev.off()


pdf("9.Progenitor_markers_featurePlot_spatial.pdf",
    width = 10, height = 10)
for (i in genes) {
  p <- HSD_spatial_featureplot_gene(seurat_obj = seurat_out,
                                    features = i,
                                    panel_width = 13,
                                    title = i,
                                    assay = "SCT")
  grid.newpage()
  grid.draw(p)
}
dev.off()

seurat_out$seurat_clusters
seurat_out$Cluster_4_15 <- ifelse(seurat_out$seurat_clusters %in% c("4", "15"),
                                  "4_15", "Other cluster")
Idents(seurat_out) <- seurat_out$Cluster_4_15
Diff_4_15 <- FindMarkers(seurat_out,
                         ident.1 = "4_15",
                         ident.2 = "Other cluster",
                         only.pos = T)
Diff_4_15 <- Diff_4_15[order(Diff_4_15$avg_log2FC,
                             decreasing = T),]
Diff_4_15$gene <- rownames(Diff_4_15)
openxlsx::write.xlsx(Diff_4_15, "Cluster_4&15_高表达基因.xlsx")

Idents(seurat_out) <- seurat_out$seurat_clusters
Diff_4 <- FindMarkers(seurat_out,
                         ident.1 = "4",
                         only.pos = T)
Diff_4 <- Diff_4[order(Diff_4$avg_log2FC,
                             decreasing = T),]
Diff_4$gene <- rownames(Diff_4)
openxlsx::write.xlsx(Diff_4, "Cluster_4_高表达基因.xlsx")

Diff_15 <- FindMarkers(seurat_out,
                      ident.1 = "15",
                      only.pos = T)
Diff_15 <- Diff_15[order(Diff_15$avg_log2FC,
                       decreasing = T),]
Diff_15$gene <- rownames(Diff_15)
openxlsx::write.xlsx(Diff_15, "Cluster_15_高表达基因.xlsx")

Diff_15_4 <- FindMarkers(seurat_out,
                         ident.1 = "15",
                         ident.2 = "4",
                         only.pos = F)
Diff_15_4 <- Diff_15_4[order(Diff_15_4$avg_log2FC,
                         decreasing = T),]
Diff_15_4$gene <- rownames(Diff_15_4)
Diff_15_4$cluster <- ifelse(Diff_15_4$avg_log2FC > 0,
                            "Cluster 15", "Cluster 4")
Diff_15_4 <- as.data.frame(Diff_15_4)
Diff_15_4 <- split.data.frame(Diff_15_4, f = list(Diff_15_4$cluster))
Diff_15_4[["Cluster 4"]] <- Diff_15_4[["Cluster 4"]][order(Diff_15_4[["Cluster 4"]]$avg_log2FC, decreasing = F),]
openxlsx::write.xlsx(Diff_15_4, "Cluster_15_vs_4_高表达基因.xlsx")

DefaultAssay(seurat_out) <- "SCT"
Spatial_cluster_diff <- FindAllMarkers(seurat_out,
                                       only.pos = T,
                                       logfc.threshold = 0.5)
Spatial_cluster_diff <- Spatial_cluster_diff[Spatial_cluster_diff$p_val_adj < 0.05,]
Spatial_cluster_diff_list <- split.data.frame(Spatial_cluster_diff,
                                              f = list(Spatial_cluster_diff$cluster))
Spatial_cluster_diff_list <- lapply(Spatial_cluster_diff_list,
                                    function(x){
                                      x[order(x$avg_log2FC, decreasing = T),]
                                    })
openxlsx::write.xlsx(Spatial_cluster_diff_list,
                     "Spatial_cluster_DEGs.xlsx")

# 空间芯片手动注释
Idents(seurat_out)
seurat_out <- RenameIdents(seurat_out,
                           "0" = "Chondrocytes",
                           "1" = "Hypertrophic chondrocytes",
                           "2" = "Endothelial cells",
                           "3" = "Proliferative_progenitor cells",
                           "4" = "Osteoblasts",
                           "5" = "Chondroclasts",
                           "6" = "Endothelial cells",
                           "7" = "Monocytes_Macrophages",
                           "8" = "Chondroclasts",
                           "9" = "AnMCs",
                           "10" = "Chondroclasts",
                           "11" = "Progenitor cells",
                           "12" = "PTCH1+ cells",
                           "13" = "Osteoblasts",
                           "14" = "Progenitor cells",
                           "15" = "Osteoblasts",
                           "16" = "Mural cells"
                           )

seurat_out$Celltype <- as.character(Idents(seurat_out))
seurat_out$Celltype <- factor(seurat_out$Celltype,
                              levels = sort(unique(seurat_out$Celltype)))
Celltype_color <- c("#A4C9DD", "#D5231E", "#399938", "#EF7C1C",
                    "#F19695", "#C6B0D2", "#F5BB6F", "#483D8B",
                    "#ADD487", "#2572A9", "#1E90FF")
names(Celltype_color) <- levels(seurat_out$Celltype)


p <- HSD_spatial_dimplot(seurat_out,
                         features = "Celltype",
                         group_color = Celltype_color)
dev.off()
pdf("10.Celltype_annotation_by_markers_spatial.pdf",
    height = 10, width = 8)
grid.newpage()
grid.draw(p)
dev.off()

p <- HSD_spatial_dimplot_highlight(seurat_obj = seurat_out,
                                   features = "Celltype",
                                   panel_width = 13,
                                   group_order = names(Celltype_color))
dev.off()
pdf("11.Celltype_annotation_by_markers_spatial_splited.pdf",
    width = 9.5, height = 12)
for (i in 1:length(p)) {
  grid.newpage()
  grid.draw(p[[i]])
}
dev.off()

p <- DimPlot(seurat_out, group.by = "Celltype",
             cols = Celltype_color, label = F, raster = TRUE)
p <- HSD_ggplot2_theme_1(p)
pdf("12.Celltype_annotation_by_markers_umap.pdf",
    width = 7, height = 4.7)
print(p)
dev.off()


p <- DimPlot(seurat_out, group.by = "Celltype", split.by = "Celltype",
             cols = Celltype_color, label = F, ncol = 3, raster = TRUE)
p <- HSD_ggplot2_theme_1(p)
pdf("12.Celltype_annotation_by_markers_splited_umap.pdf",
    width = 10, height = 10)
print(p)
dev.off()

Idents(seurat_out) <- seurat_out$Celltype
DotPlot(seurat_out,
        features = c("SFRP2", "IGF1", "PRRX1", "THY1", "TWIST2", # AnMCs
                     "ACP5", "CSTB", # Chondroclasts
                     "COMP", "COL9A3", "COL2A1", "COL9A2", "COL9A1", "IHH", # Chondrocytes
                     "PLVAP", "PECAM1", "TOP2A", "CD34", "CDH5", "ENG", # Endothelial cells
                     "COL10A1", "RUNX2", # Hypertrophic chondrocytes
                     "CD68", "CSF1R", # Monocytes_Macrophages
                     "ACTA2", "MYL9", # Mural cells
                     "IBSP", "SPP1", "PHEX", "BGLAP", "OMD", # Osteoblasts
                     "CDKN1C", "IGF1R", "TNC", "TNN", "POSTN", # Progenitor cells
                     "PTN", "CCND2", "MKI67", # Proliferative_progenitor cells
                     "PLTP", "HHIP", "PTCH1" # PTCH1+ cells
                     )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_gradientn(colours = c("white", #"#1072A7", "#50B777",
                                    "#EDAC37", "#FF0000"))


plist <- lapply(c("SFRP2", "IGF1", "PRRX1", "THY1", "TWIST2"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 3)
pdf("13.Celltype_marker_featureplot(AnMCs)_umap.pdf",
    height = 7, width = 11)
p
dev.off()

plist <- lapply(c("ACP5", "CSTB"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 2)
pdf("13.Celltype_marker_featureplot(Chondroclasts)_umap.pdf",
    height = 3.5, width = 8)
p
dev.off()

plist <- lapply(c("COMP", "COL9A3", "COL2A1",
                  "COL9A2", "COL9A1", "IHH"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 3)
pdf("13.Celltype_marker_featureplot(Chondrocytes)_umap.pdf",
    height = 7, width = 11)
p
dev.off()

plist <- lapply(c("PLVAP", "PECAM1", "TOP2A",
                  "CD34", "CDH5", "ENG"), function(x){
                    FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
                      scale_color_gradientn(colours = c("gray", "blue", "red"))
                  })
p <- cowplot::plot_grid(plotlist = plist, ncol = 3)
pdf("13.Celltype_marker_featureplot(Endothelial cells)_umap.pdf",
    height = 7, width = 11)
p
dev.off()

plist <- lapply(c("COL10A1", "RUNX2"), function(x){
                    FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
                      scale_color_gradientn(colours = c("gray", "blue", "red"))
                  })
p <- cowplot::plot_grid(plotlist = plist, ncol = 2)
pdf("13.Celltype_marker_featureplot(Hypertrophic chondrocytes)_umap.pdf",
    height = 3.5, width = 8)
p
dev.off()

plist <- lapply(c("CD68", "CSF1R"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 2)
pdf("13.Celltype_marker_featureplot(Monocytes_Macrophages)_umap.pdf",
    height = 3.5, width = 8)
p
dev.off()

plist <- lapply(c("ACTA2", "MYL9"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 2)
pdf("13.Celltype_marker_featureplot(Mural cells)_umap.pdf",
    height = 3.5, width = 8)
p
dev.off()

plist <- lapply(c("IBSP", "SPP1", "PHEX", "BGLAP", "OMD"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 3)
pdf("13.Celltype_marker_featureplot(Osteoblasts)_umap.pdf",
    height = 7, width = 11)
p
dev.off()

plist <- lapply(c("CDKN1C", "IGF1R", "TNC", "TNN", "POSTN"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 3)
pdf("13.Celltype_marker_featureplot(Progenitor cells)_umap.pdf",
    height = 7, width = 11)
p
dev.off()

plist <- lapply(c("PTN", "CCND2", "MKI67"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 3)
pdf("13.Celltype_marker_featureplot(Proliferative_progenitor cells)_umap.pdf",
    height = 3.5, width = 11)
p
dev.off()

plist <- lapply(c("PLTP", "HHIP", "PTCH1"), function(x){
  FeaturePlot(seurat_out, features = x, order = T, raster = TRUE) +
    scale_color_gradientn(colours = c("gray", "blue", "red"))
})
p <- cowplot::plot_grid(plotlist = plist, ncol = 3)
pdf("13.Celltype_marker_featureplot(PTCH1+ cells)_umap.pdf",
    height = 3.5, width = 11)
p
dev.off()


pdf("14.Celltype_marker_featureplot(PTCH1+ cells)_spatial.pdf",
    width = 10, height = 10)
for (i in c("PLTP", "HHIP", "PTCH1")) {
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

pdf("14.Celltype_marker_featureplot(Proliferative_progenitor cells)_spatial.pdf",
    width = 10, height = 10)
for (i in c("PTN", "CCND2", "MKI67")) {
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

pdf("14.Celltype_marker_featureplot(Progenitor cells)_spatial.pdf",
    width = 10, height = 10)
for (i in c("CDKN1C", "IGF1R", "TNC", "TNN", "POSTN")) {
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

pdf("14.Celltype_marker_featureplot(Osteoblasts)_spatial.pdf",
    width = 10, height = 10)
for (i in c("IBSP", "SPP1", "PHEX", "BGLAP", "OMD")) {
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

pdf("14.Celltype_marker_featureplot(Mural cells)_spatial.pdf",
    width = 10, height = 10)
for (i in c("ACTA2", "MYL9")) {
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

pdf("14.Celltype_marker_featureplot(Monocytes_Macrophages)_spatial.pdf",
    width = 10, height = 10)
for (i in c("CD68", "CSF1R")) {
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

pdf("14.Celltype_marker_featureplot(Hypertrophic chondrocytes)_spatial.pdf",
    width = 10, height = 10)
for (i in c("COL10A1", "RUNX2")) {
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

pdf("14.Celltype_marker_featureplot(Endothelial cells)_spatial.pdf",
    width = 10, height = 10)
for (i in c("PLVAP", "PECAM1", "TOP2A",
            "CD34", "CDH5", "ENG")) {
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

pdf("14.Celltype_marker_featureplot(Chondrocytes)_spatial.pdf",
    width = 10, height = 10)
for (i in c("COMP", "COL9A3", "COL2A1",
            "COL9A2", "COL9A1", "IHH")) {
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

pdf("14.Celltype_marker_featureplot(Chondroclasts)_spatial.pdf",
    width = 10, height = 10)
for (i in c("ACP5", "CSTB")) {
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

saveRDS(seurat_out, "seurat_out_annotated.rds")


setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/Stereo-seq/pythonProject")
seurat_out <- readRDS("seurat_out_annotated.rds")
seurat_out_os <- subset(seurat_out, subset = Celltype %in% c("Chondrocytes", "Hypertrophic chondrocytes",
                                                             "PTCH1+ cells", "Osteoblasts"))
##### Monocle2
library(monocle)

DefaultAssay(seurat_out) <- "Spatial"
seurat_out <- magic(seurat_out)

seurat_out



