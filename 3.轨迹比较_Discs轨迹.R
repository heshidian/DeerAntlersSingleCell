library(Seurat)
library(dplyr)
library(monocle3)
Discs <- readRDS("./OtherDataSets/disc_sce.0.3.integrate_umap_last.3000.RDS")
DimPlot(Discs, label = T)
levels(Discs$celltype)
Discs_sub <- subset(Discs, celltype %in% c("HomCs", "preHTCs", "HTCs"))
DimPlot(Discs_sub, label = T)
Discs_sub_color <- c("#A4C9DD", "#54989E", "#51A743",
                     "#ED8684", "#E46A45", "#F0841F",
                     "#AD91C0", "#C4B497", "#AA5728")
names(Discs_sub_color) <- levels(Discs$celltype)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Discs_sub <- CellCycleScoring(Discs_sub,
                              s.features = s.genes,
                              g2m.features = g2m.genes,
                              set.ident = FALSE,
                              search = TRUE)

table(Discs_sub$celltype,
      Discs_sub$Phase)
Discs_sub$Phase <- factor(as.character(Discs_sub$Phase),
                          levels = c("G1", "S", "G2M"))
pdf("./轨迹比较/Discs_Cellcycle_phase.pdf",
    height = 6, width = 5.5)
percent_bar(Discs_sub@meta.data,
            Ident = "Phase",
            Group = "celltype",
            fill_color = c("#009999", "#FFCC00", "#CC3366"),
            fill_label = "",
            flow = FALSE)
dev.off()


library(CytoTRACE2)
cytotrace2_result <- cytotrace2(Discs_sub,   
                                species = "human",
                                is_seurat = TRUE,
                                slot_type = "counts",
                                full_model = FALSE,
                                batch_size = 10000,
                                smooth_batch_size = 1000,
                                parallelize_models = FALSE,
                                parallelize_smoothing = FALSE,
                                ncores = 20,
                                seed = 14)

pdf("./轨迹比较/Discs_cytoTrace2.pdf", width = 3, height = 5.5)
violin_linebox(meta = cytotrace2_result@meta.data,
               group.by = "celltype",
               features = c("CytoTRACE2_Score"),
               group.color = Discs_sub_color,
               linewidth = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none")
dev.off()


dim_use <- 30
groups <- "celltype"
batch <- "orig.ident"
reduction_method <- "UMAP"
data <- GetAssayData(Discs_sub, assay = 'RNA', slot = 'counts')
cell_metadata <- Discs_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
Discs_monocle3_cds <- new_cell_data_set(data,
                                        cell_metadata = cell_metadata,
                                        gene_metadata = gene_annotation)
Discs_monocle3_cds <- preprocess_cds(Discs_monocle3_cds, num_dim = dim_use)
Discs_monocle3_cds <- align_cds(Discs_monocle3_cds, alignment_group = batch,
                                alignment_k = dim_use)
Discs_monocle3_cds <- reduce_dimension(Discs_monocle3_cds,
                                       reduction_method = "UMAP",
                                       preprocess_method = "Aligned",
                                       umap.min_dist = 0.1,
                                       umap.n_neighbors = 30)
cds.embed <- Discs_monocle3_cds@int_colData$reducedDims$UMAP
colnames(cds.embed) <- c("UMAP_1", "UMAP_2")
cds.embed <- as.data.frame(cds.embed)
cds.embed$Celltype <- Discs_sub@meta.data[rownames(cds.embed), groups]
ggplot(data = cds.embed, aes(x = UMAP_1, y = UMAP_2, colour = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = Discs_sub_color)
cds.embed <- Discs_monocle3_cds@int_colData$reducedDims$PCA
cds.embed <- as.data.frame(cds.embed)
cds.embed$Celltype <- Discs_sub@meta.data[rownames(cds.embed), groups]
ggplot(data = cds.embed, aes(x = PC1, y = PC3, colour = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = Discs_sub_color)
int.embed <- Embeddings(Discs_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
Discs_monocle3_cds@int_colData$reducedDims$UMAP <- int.embed
Discs_monocle3_cds <- cluster_cells(Discs_monocle3_cds,
                                    reduction_method = "UMAP",
                                    k = dim_use)
Discs_monocle3_cds <- learn_graph(Discs_monocle3_cds,
                                  use_partition = F,
                                  close_loop = T,
                                  learn_graph_control = list(nn.k = 20, # 15
                                                             minimal_branch_len = 16)) # 16
umap <- Embeddings(Discs, reduction = "umap")
umap <- as.data.frame(umap)
p <- plot_cells_2(cds = Discs_monocle3_cds,
                  reduction_method = reduction_method, # 可以选择TSNE， PCA
                  color_cells_by = groups, # 展示细胞的分组，这里设置的是Sample
                  show_trajectory_graph = TRUE, trajectory_graph_color = "grey28", trajectory_graph_segment_size = 1.5, # 展示轨迹
                  group_cells_color = Discs_sub_color, # 设置细胞分组的颜色，默认是NULL，可以自己设置
                  label_cell_groups = FALSE, label_groups_by_cluster = TRUE, # 设置是否展示细胞分组
                  group_label_size = 6, # 设置分组label的大小
                  label_principal_points = FALSE, label_leaves = FALSE, graph_label_size = 0, # 设置轨迹上的分至点以及叶节点
                  cell_size = 0.8, # 设置细胞大小
                  background_cell = umap,
                  background_cell_size = 0.8
)
pdf(paste0("./轨迹比较/Discs_Monocle3_Step1_",
           reduction_method,
           "_", groups, ".pdf"),
    width = 6.5, height = 5)
print(p)
dev.off()


Discs_monocle3_cds <- Monocle3_Step2(CDS = Discs_monocle3_cds, # cds 数据对象，上一步得到的结果
                                     reduction_method = reduction_method, # 降维方法
                                     root_cells = NULL # 手动设置根细胞，默认为NULL，即用鼠标圈选。如果不想用鼠标圈，可以直接输入细胞的barcode，比如 TTTCACAGTAGTCGTT-2，可以输入多个细胞作为起点
)
#### 可视化拟时序
p <- plot_cells_2(cds = Discs_monocle3_cds,
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
pdf(paste0("./轨迹比较/Discs_Monocle3_Step2_",
           reduction_method, "_pseudotime", ".pdf"),
    width = 6.8, height = 5)
print(p)
dev.off()

saveRDS(Discs_monocle3_cds,
        "3.轨迹比较.Discs_monocle3_cds.rds")

df <- colData(Discs_monocle3_cds)
df <- as.data.frame(df)
df$Pseudotime <- pseudotime(Discs_monocle3_cds)
df$celltype <- factor(as.character(df$celltype),
                      levels = c("HomCs", "preHTCs", "HTCs"))
p <- jitter_linebox(meta = df, group.by = "celltype",
                    features = "Pseudotime",
                    group.color = Discs_sub_color,
                    linewidth = 2, jitter_width = 0.2, stat_cor_size = 4) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10)) +
  labs(x = "", y = "Pseudotime")
pdf("./轨迹比较/Discs_Monocle3_pesudotime_celltype.pdf",
    width = 3, height = 5)
print(p)
dev.off()


