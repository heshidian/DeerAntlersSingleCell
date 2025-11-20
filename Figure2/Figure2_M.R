

seurat_out_os_PTCH1_3 <- subset(seurat_out_os_PTCH1,
                                subset = Osteoblasts_sub %in% c("Chondrocytes",
                                                                "Hypertrophic chondrocytes",
                                                                "Osteoblasts (PTCH1+ cells)",
                                                                # "Osteoblasts (Core)",
                                                                # "Osteoblasts (Endothelial cells)",
                                                                "PTCH1+ cells"))
source("../../monocle3_util.R")
library(monocle3)
#### 参数设置
dim_use <- 30
groups <- "Osteoblasts_sub"
batch <- "orig.ident"
reduction_method <- "UMAP"

data <- GetAssayData(seurat_out_os_PTCH1_3, assay = 'SCT', slot = 'counts')
cell_metadata <- seurat_out_os_PTCH1_3@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
monocle3_cds <- new_cell_data_set(data,
                                  cell_metadata = cell_metadata,
                                  gene_metadata = gene_annotation)
monocle3_cds <- preprocess_cds(monocle3_cds, num_dim = dim_use)
monocle3_cds <- align_cds(monocle3_cds,
                          alignment_k = dim_use)
monocle3_cds <- reduce_dimension(monocle3_cds,
                                 reduction_method = "UMAP",
                                 preprocess_method = "PCA",
                                 umap.min_dist = 0.5,
                                 umap.n_neighbors = 30)
cds.embed <- monocle3_cds@int_colData$reducedDims$UMAP
colnames(cds.embed) <- c("UMAP_1", "UMAP_2")
cds.embed <- as.data.frame(cds.embed)
cds.embed$Celltype <- seurat_out_os_PTCH1_3@meta.data[rownames(cds.embed), "Osteoblasts_sub"]
ggplot(data = cds.embed, aes(x = UMAP_1, y = UMAP_2, colour = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = Celltype_color_3)
int.embed <- Embeddings(seurat_out_os_PTCH1_2, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
monocle3_cds@int_colData$reducedDims$UMAP <- int.embed
monocle3_cds <- cluster_cells(monocle3_cds,
                              reduction_method = "UMAP",
                              k = dim_use)
monocle3_cds <- learn_graph(monocle3_cds,
                            use_partition = F,
                            close_loop = T,
                            learn_graph_control = list(nn.k = 20, # 20
                                                       minimal_branch_len = 6)) # 5

umap <- Embeddings(seurat_out_os_PTCH1_2, reduction = "umap")
umap <- as.data.frame(umap)
p <- plot_cells_2(cds = monocle3_cds,
                  reduction_method = reduction_method, # 可以选择TSNE， PCA
                  color_cells_by = groups, # 展示细胞的分组，这里设置的是Sample
                  show_trajectory_graph = TRUE, trajectory_graph_color = "grey28", trajectory_graph_segment_size = 1.5, # 展示轨迹
                  group_cells_color = Celltype_color_3, # 设置细胞分组的颜色，默认是NULL，可以自己设置
                  label_cell_groups = FALSE, label_groups_by_cluster = TRUE, # 设置是否展示细胞分组
                  group_label_size = 6, # 设置分组label的大小
                  label_principal_points = FALSE, label_leaves = FALSE, graph_label_size = 0, # 设置轨迹上的分至点以及叶节点
                  cell_size = 0.8, # 设置细胞大小
                  background_cell = umap,
                  background_cell_size = 0.8
)

pdf("../骨轨迹/Osteoblasts_PTCH1_traj_V2.pdf",
    width = 8.3, height = 5)
print(p)
dev.off()