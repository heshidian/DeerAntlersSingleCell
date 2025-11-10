library(Seurat)
library(monocle)
library(monocle3)
library(harmony)
library(tidyverse)
library(patchwork)
library(BSgenome.malu)
library(ArchR)
library(ggalt)
library(Seurat)
library(parallel)
library(clustree)
library(dplyr)
library(patchwork)
library(BSgenome.malu)
library(slingshot)
library(SingleCellExperiment)
library(qs)
library(tidyverse)
library(RColorBrewer)
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)

setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
Limb <- readRDS("./OtherDataSets/sce.0.13.integrate_umap_last.3000.RDS")
DimPlot(Limb, label = T)
levels(Limb$celltype)
Limb_sub <- subset(Limb, celltype %in% c("OCPs", "Osteoprogenitors",
                                         "PMSC", "Chondroblasts",
                                         "Chondrocytes"))
Limb_sub_color <- c("#a6cee3", "#2f83af", "#95d175",
                    "#759e50", "#399938")
names(Limb_sub_color) <- c("OCPs", "Osteoprogenitors",
                           "PMSC", "Chondroblasts",
                           "Chondrocytes")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Limb_sub <- CellCycleScoring(Limb_sub,
                              s.features = s.genes,
                              g2m.features = g2m.genes,
                              set.ident = FALSE)

table(Limb_sub$celltype,
      Limb_sub$Phase)
Limb_sub$Phase <- factor(as.character(Limb_sub$Phase),
                           levels = c("G1", "S", "G2M"))
pdf("./轨迹比较/Limb_Cellcycle_phase.pdf",
    height = 5, width = 5)
percent_bar(Limb_sub@meta.data,
            Ident = "Phase",
            Group = "celltype",
            fill_color = c("#009999", "#FFCC00", "#CC3366"),
            fill_label = "",
            flow = FALSE)
dev.off()


library(CytoTRACE2)
cytotrace2 <- function (input, species = "mouse", is_seurat = FALSE, slot_type = "counts", 
                        full_model = FALSE, batch_size = 10000, smooth_batch_size = 1000, 
                        parallelize_models = TRUE, parallelize_smoothing = TRUE, 
                        ncores = NULL, max_pcs = 200, seed = 14) 
{
  set.seed(seed)
  message("cytotrace2: Started loading data")
  if (class(input) == "Seurat" & is_seurat == FALSE) {
    stop("The input is a Seurat object. Please make sure to set is_seurat = TRUE.")
  }
  if (is.character(input) && file.exists(input)) {
    if (is_seurat == FALSE) {
      data <- loadData(input)
    } else {
      seurat <- readRDS(input)
      data <- loadData_fromSeurat(seurat, slot_type)
    }
  } else {
    if (is_seurat == FALSE) {
      data <- copy(input)
    } else {
      seurat <- copy(input)
      data <- loadData_fromSeurat(input, slot_type)
    }
  }
  if (any(duplicated(colnames(data)))) {
    stop("Please make sure the cell/sample names are unique")
  }
  if (any(duplicated(rownames(data)))) {
    stop("Please make sure the gene names are unique")
  }
  message("Dataset contains ", dim(data)[1], " genes and ", 
          dim(data)[2], " cells.")
  is_human <- sum(sapply(rownames(data), function(x) all(toupper(x) == 
                                                           x)))/nrow(data) > 0.9
  is_mouse <- sum(sapply(rownames(data), function(x) all(toupper(x) != 
                                                           x)))/nrow(data) > 0.9
  if (is_human & species == "mouse") {
    warning("Species is most likely human. Please revise the 'species' input to the function.")
  }
  if (is_mouse & species == "human") {
    warning("Species is most likely mouse. Please revise the 'species' input to the function.")
  }
  if (is.null(batch_size)) {
    batch_size <- ncol(data)
  }
  if (is.null(smooth_batch_size)) {
    smooth_batch_size <- ncol(data)
    parallelize_smoothing == FALSE
  }
  if (ncol(data) <= 1000 && parallelize_smoothing == TRUE) {
    parallelize_smoothing <- FALSE
    message("The number of cells in your dataset is less than 1000. Fast mode has been disabled.")
  }
  if (parallelize_smoothing || parallelize_models) {
    if (Sys.info()["sysname"] == "Windows") {
      ncores = 1
      message("Windows OS can run only on 1 core")
    }
    else {
      if (is.null(ncores)) {
        ncores <- parallel::detectCores(all.tests = FALSE, 
                                        logical = TRUE) - 1
      }
    }
  } else {
    if (is.null(ncores)) {
      ncores <- 1
    }
  }
  if (batch_size > ncol(data)) {
    message("The passed subsample size is greater than the number of cells in dataset.\nNow setting subsample size to ", 
            ncol(data))
    batch_size <- ncol(data)
  } else if (ncol(data) > 10000 && batch_size > 10000) {
    message("Please consider reducing the batch_size to 10000 for runtime and memory efficiency.")
  }
  init_order <- colnames(data)
  chunk <- round(ncol(data)/batch_size)
  subsamples <- split(seq_len(ncol(data)), sample(factor(seq_len(ncol(data))%%chunk)))
  sample_names <- lapply(subsamples, function(x) colnames(data)[x])
  message("cytotrace2: Running on ", chunk, " subsample(s) approximately of length ", 
          batch_size)
  if (full_model) {
    parameter_dict <- readRDS(system.file("extdata", "parameter_dict_17.rds", 
                                          package = "CytoTRACE2"))
  } else {
    parameter_dict <- readRDS(system.file("extdata", "parameter_dict_5_best.rds", 
                                          package = "CytoTRACE2"))
  }
  nc <- min(length(parameter_dict), ncores)
  subsample_processing_f <- function(subsample) {
    dt <- data[, subsample]
    message("cytotrace2: Started preprocessing.")
    ranked_data <- preprocessData(dt, species)
    gene_names <- colnames(ranked_data)
    cell_names <- rownames(ranked_data)
    message("cytotrace2: Started prediction.")
    predicted_df <- predictData(parameter_dict, ranked_data, 
                                parallelize_models, ncores = nc)
    predicted_df <- predicted_df[cell_names, ]
    num_genes <- ncol(ranked_data)
    ranked_data <- as.matrix(ranked_data)
    disp_fn <- function (x) {
      if (length(unique(x)) == 1) {
        return(0)
      }
      else {
        return(stats::var(x)/mean(x))
      }
    }
    dispersion_index <- sapply(1:num_genes, function(i) disp_fn(ranked_data[, i]))
    top_genes <- gene_names[order(dispersion_index, decreasing = TRUE)[1:min(1000, 
                                                                             num_genes)]]
    message("cytotrace2: Started postprocessing.")
    smoothScore <- smoothData(ranked_data, predicted_df, 
                              top_genes, ncores = ncores, smooth_batch_size = smooth_batch_size, 
                              parallelize_smoothing = parallelize_smoothing, seed = seed)
    smoothScore <- smoothScore[cell_names]
    predicted_df$preKNN_CytoTRACE2_Score <- smoothScore
    if (nrow(ranked_data) <= 10) {
      message("cytotrace2: Number of cells in data is less than 10. Skipping postprocessing.")
      predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
      predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
    }
    else {
      predicted_df <- binData(predicted_df)
      predicted_df <- predicted_df[cell_names, ]
      if (nrow(ranked_data) <= 100) {
        message("cytotrace2: Number of cells in data is less than 100. Skipping kNN smoothing.")
        predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
        predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
      }
      if (sd(ranked_data) == 0) {
        message("cytotrace2: Zero variance of ranked matrix. Skipping kNN smoothing.")
        predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
        predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
      }
      if (nrow(ranked_data) > 100 && sd(ranked_data) != 
          0) {
        ranked_df <- data.frame(base::t(ranked_data))
        rownames(ranked_df) <- gene_names
        colnames(ranked_df) <- cell_names
        predicted_df <- smoothDatakNN(ranked_df, predicted_df, 
                                      top_genes, max_pcs, seed)
        predicted_df <- predicted_df[cell_names, ]
      }
    }
    predicted_df
  }
  message("cytotrace2: Started running on subsample(s). This will take a few minutes.")
  results <- lapply(subsamples, subsample_processing_f)
  predicted_df <- do.call(rbind, results)
  rownames(predicted_df) <- unlist(sample_names)
  predicted_df <- predicted_df[init_order, ]
  ranked_scores <- rank(predicted_df$CytoTRACE2_Score)
  predicted_df$CytoTRACE2_Relative <- (ranked_scores - min(ranked_scores))/(max(ranked_scores) - 
                                                                              min(ranked_scores))
  predicted_df <- predicted_df[c("CytoTRACE2_Score", "CytoTRACE2_Potency", 
                                 "CytoTRACE2_Relative", "preKNN_CytoTRACE2_Score", "preKNN_CytoTRACE2_Potency")]
  if (!is_seurat) {
    message("cytotrace2: Finished")
    return(predicted_df)
  }
  else {
    predicted_df <- predicted_df[colnames(seurat), ]
    seurat <- AddMetaData(object = seurat, metadata = predicted_df)
    message("cytotrace2: Finished")
    return(seurat)
  }
}
cytotrace2_result <- cytotrace2(Limb_sub,   
                                species = "human",
                                is_seurat = TRUE,
                                slot_type = "data",
                                full_model = FALSE,
                                # batch_size = 100000,
                                # smooth_batch_size = 10000,
                                # parallelize_models = TRUE,
                                # parallelize_smoothing = TRUE,
                                ncores = 20,
                                seed = 14)

pdf("./轨迹比较/Limb_cytoTrace2.pdf", width = 5, height = 5.5)
violin_linebox(meta = cytotrace2_result@meta.data,
               group.by = "celltype",
               features = c("CytoTRACE2_Score"),
               group.color = Limb_sub_color,
               linewidth = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none")
dev.off()

dim_use <- 30
groups <- "celltype"
batch <- "orig.ident"
reduction_method <- "UMAP"
data <- GetAssayData(Limb_sub, assay = 'RNA', slot = 'counts')
cell_metadata <- Limb_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
Limb_monocle3_cds <- new_cell_data_set(data,
                                         cell_metadata = cell_metadata,
                                         gene_metadata = gene_annotation)
Limb_monocle3_cds <- preprocess_cds(Limb_monocle3_cds, num_dim = dim_use)
Limb_monocle3_cds <- align_cds(Limb_monocle3_cds, alignment_group = batch,
                                 alignment_k = dim_use)
Limb_monocle3_cds <- reduce_dimension(Limb_monocle3_cds,
                                        reduction_method = "UMAP",
                                        preprocess_method = "Aligned",
                                        umap.min_dist = 0.1,
                                        umap.n_neighbors = 30)
cds.embed <- Limb_monocle3_cds@int_colData$reducedDims$UMAP
colnames(cds.embed) <- c("UMAP_1", "UMAP_2")
cds.embed <- as.data.frame(cds.embed)
cds.embed$Celltype <- Limb_sub@meta.data[rownames(cds.embed), groups]
ggplot(data = cds.embed, aes(x = UMAP_1, y = UMAP_2, colour = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = Limb_sub_color)
cds.embed <- Limb_monocle3_cds@int_colData$reducedDims$PCA
cds.embed <- as.data.frame(cds.embed)
cds.embed$Celltype <- Limb_sub@meta.data[rownames(cds.embed), groups]
ggplot(data = cds.embed, aes(x = PC2, y = PC3, colour = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = Limb_sub_color)
int.embed <- Embeddings(Limb_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
Limb_monocle3_cds@int_colData$reducedDims$UMAP <- int.embed
Limb_monocle3_cds <- cluster_cells(Limb_monocle3_cds,
                                     reduction_method = "UMAP",
                                     k = dim_use)
Limb_monocle3_cds <- learn_graph(Limb_monocle3_cds,
                                   use_partition = F,
                                   close_loop = T,
                                   learn_graph_control = list(nn.k = 15, # 15
                                                              minimal_branch_len = 10)) # 16
umap <- Embeddings(Limb, reduction = "umap")
umap <- as.data.frame(umap)
p <- plot_cells_2(cds = Limb_monocle3_cds,
                  reduction_method = reduction_method, # 可以选择TSNE， PCA
                  color_cells_by = groups, # 展示细胞的分组，这里设置的是Sample
                  show_trajectory_graph = TRUE, trajectory_graph_color = "grey28", trajectory_graph_segment_size = 1.5, # 展示轨迹
                  group_cells_color = Limb_sub_color, # 设置细胞分组的颜色，默认是NULL，可以自己设置
                  label_cell_groups = FALSE, label_groups_by_cluster = TRUE, # 设置是否展示细胞分组
                  group_label_size = 6, # 设置分组label的大小
                  label_principal_points = FALSE, label_leaves = FALSE, graph_label_size = 0, # 设置轨迹上的分至点以及叶节点
                  cell_size = 0.8, # 设置细胞大小
                  background_cell = umap,
                  background_cell_size = 0.8
)
pdf(paste0("./轨迹比较/Limb_Monocle3_Step1_",
           reduction_method,
           "_", groups, ".pdf"),
    width = 7.2, height = 5)
print(p)
dev.off()

Limb_monocle3_cds <- Monocle3_Step2(CDS = Limb_monocle3_cds, # cds 数据对象，上一步得到的结果
                                      reduction_method = reduction_method, # 降维方法
                                      root_cells = NULL # 手动设置根细胞，默认为NULL，即用鼠标圈选。如果不想用鼠标圈，可以直接输入细胞的barcode，比如 TTTCACAGTAGTCGTT-2，可以输入多个细胞作为起点
)
#### 可视化拟时序
p <- plot_cells_2(cds = Limb_monocle3_cds,
                  reduction_method = reduction_method, # 可以选择TSNE， PCA
                  color_cells_by = "pseudotime", # 展示细胞的分组，这里设置的是Sample
                  show_trajectory_graph = TRUE, trajectory_graph_color = "grey28", trajectory_graph_segment_size = 1.5, # 展示轨迹
                  group_cells_color = NULL, # 设置细胞分组的颜色，默认是NULL，可以自己设置
                  label_cell_groups = FALSE, label_groups_by_cluster = TRUE, # 设置是否展示细胞分组
                  group_label_size = 6, # 设置分组label的大小
                  label_principal_points = FALSE, label_leaves = FALSE, graph_label_size = 0, # 设置轨迹上的分至点以及叶节点
                  cell_size = 0.8, # 设置细胞大小
                  background_cell = umap,
                  background_cell_size = 0.8
)
pdf(paste0("./轨迹比较/Limb_Monocle3_Step2_",
           reduction_method, "_pseudotime", ".pdf"),
    width = 6.8, height = 5)
print(p)
dev.off()

saveRDS(Limb_monocle3_cds,
        "3.轨迹比较.Limb_monocle3_cds.rds")


