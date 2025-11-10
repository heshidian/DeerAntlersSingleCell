Tumor <- readRDS("./OtherDataSets/tumor_sce.0.15.integrate_umap_last.3000.ham.RDS")
DimPlot(Tumor, label = T)
levels(Tumor$celltype)
Tumor_sub <- subset(Tumor, celltype %in% c("Tumor_Osteoblastic", "Tumor_Chondroblasts"))
DimPlot(Tumor_sub, label = T)
Tumor_sub_color <- c("#379938", "#F19695")
names(Tumor_sub_color) <- c("Tumor_Osteoblastic", "Tumor_Chondroblasts")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Tumor_sub <- CellCycleScoring(Tumor_sub,
                               s.features = s.genes,
                               g2m.features = g2m.genes,
                               set.ident = FALSE)

table(Tumor_sub$celltype,
      Tumor_sub$Phase)
Tumor_sub$Phase <- factor(as.character(Tumor_sub$Phase),
                           levels = c("G1", "S", "G2M"))
pdf("./轨迹比较/Tumor_Cellcycle_phase.pdf",
    height = 6, width = 4)
percent_bar(Tumor_sub@meta.data,
            Ident = "Phase",
            Group = "celltype",
            fill_color = c("#009999", "#FFCC00", "#CC3366"),
            fill_label = "",
            flow = FALSE)
dev.off()


library(CytoTRACE2)
cytotrace2_result <- cytotrace2(Tumor_sub,   
                                species = "human",
                                is_seurat = TRUE,
                                slot_type = "counts",
                                full_model = FALSE,
                                batch_size = 10000,
                                smooth_batch_size = 1000,
                                parallelize_models = TRUE,
                                parallelize_smoothing = TRUE,
                                ncores = 20,
                                seed = 14)

pdf("./轨迹比较/Tumor_cytoTrace2.pdf", width = 3, height = 5.5)
violin_linebox(meta = cytotrace2_result@meta.data,
               group.by = "celltype",
               features = c("CytoTRACE2_Score"),
               group.color = Tumor_sub_color,
               linewidth = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none")
dev.off()


dim_use <- 30
groups <- "celltype"
batch <- "orig.ident"
reduction_method <- "UMAP"
data <- GetAssayData(Tumor_sub, assay = 'RNA', slot = 'counts')
cell_metadata <- Tumor_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
Tumor_monocle3_cds <- new_cell_data_set(data,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_annotation)
Tumor_monocle3_cds <- preprocess_cds(Tumor_monocle3_cds, num_dim = dim_use)
Tumor_monocle3_cds <- align_cds(Tumor_monocle3_cds, alignment_group = batch,
                               alignment_k = dim_use)
Tumor_monocle3_cds <- reduce_dimension(Tumor_monocle3_cds,
                                      reduction_method = "UMAP",
                                      preprocess_method = "Aligned",
                                      umap.min_dist = 0.1,
                                      umap.n_neighbors = 30)
cds.embed <- Tumor_monocle3_cds@int_colData$reducedDims$UMAP
colnames(cds.embed) <- c("UMAP_1", "UMAP_2")
cds.embed <- as.data.frame(cds.embed)
cds.embed$Celltype <- Tumor_sub@meta.data[rownames(cds.embed), groups]
ggplot(data = cds.embed, aes(x = UMAP_1, y = UMAP_2, colour = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = Tumor_sub_color)
cds.embed <- Tumor_monocle3_cds@int_colData$reducedDims$PCA
cds.embed <- as.data.frame(cds.embed)
cds.embed$Celltype <- Tumor_sub@meta.data[rownames(cds.embed), groups]
ggplot(data = cds.embed, aes(x = PC1, y = PC3, colour = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = Tumor_sub_color)
int.embed <- Embeddings(Tumor_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
Tumor_monocle3_cds@int_colData$reducedDims$UMAP <- int.embed
Tumor_monocle3_cds <- cluster_cells(Tumor_monocle3_cds,
                                   reduction_method = "UMAP",
                                   k = dim_use)
Tumor_monocle3_cds <- learn_graph(Tumor_monocle3_cds,
                                 use_partition = F,
                                 close_loop = T,
                                 learn_graph_control = list(nn.k = 15, # 15
                                                            minimal_branch_len = 10)) # 16
umap <- Embeddings(Tumor, reduction = "umap")
umap <- as.data.frame(umap)
p <- plot_cells_2(cds = Tumor_monocle3_cds,
                  reduction_method = reduction_method, # 可以选择TSNE， PCA
                  color_cells_by = groups, # 展示细胞的分组，这里设置的是Sample
                  show_trajectory_graph = TRUE, trajectory_graph_color = "grey28", trajectory_graph_segment_size = 1.5, # 展示轨迹
                  group_cells_color = Tumor_sub_color, # 设置细胞分组的颜色，默认是NULL，可以自己设置
                  label_cell_groups = FALSE, label_groups_by_cluster = TRUE, # 设置是否展示细胞分组
                  group_label_size = 6, # 设置分组label的大小
                  label_principal_points = FALSE, label_leaves = FALSE, graph_label_size = 0, # 设置轨迹上的分至点以及叶节点
                  cell_size = 0.8, # 设置细胞大小
                  background_cell = umap,
                  background_cell_size = 0.8
)
pdf(paste0("./轨迹比较/Tumor_Monocle3_Step1_",
           reduction_method,
           "_", groups, ".pdf"),
    width = 7.5, height = 5)
print(p)
dev.off()

Tumor_monocle3_cds@colData@listData[["CytoTRACE2_Score"]] <- cytotrace2_result@meta.data[colnames(Tumor_monocle3_cds), "CytoTRACE2_Score"]
pdf("./轨迹比较/Tumor_Monocle3_CytoTRACE2_Score.pdf",
    height = 5, width = 7.5)
plot_cells_2(cds = Tumor_monocle3_cds,
             reduction_method = reduction_method, # 可以选择TSNE， PCA
             color_cells_by = "CytoTRACE2_Score", # 展示细胞的分组，这里设置的是Sample
             show_trajectory_graph = TRUE, trajectory_graph_color = "grey28", trajectory_graph_segment_size = 1.5, # 展示轨迹
             group_cells_color = NULL, # 设置细胞分组的颜色，默认是NULL，可以自己设置
             label_cell_groups = FALSE, label_groups_by_cluster = TRUE, # 设置是否展示细胞分组
             group_label_size = 6, # 设置分组label的大小
             label_principal_points = FALSE, label_leaves = FALSE, graph_label_size = 0, # 设置轨迹上的分至点以及叶节点
             cell_size = 0.8, # 设置细胞大小
             background_cell = umap,
             background_cell_size = 0.8
)
dev.off()

Tumor_monocle3_cds <- Monocle3_Step2(CDS = Tumor_monocle3_cds, # cds 数据对象，上一步得到的结果
                                    reduction_method = reduction_method, # 降维方法
                                    root_cells = NULL # 手动设置根细胞，默认为NULL，即用鼠标圈选。如果不想用鼠标圈，可以直接输入细胞的barcode，比如 TTTCACAGTAGTCGTT-2，可以输入多个细胞作为起点
)
#### 可视化拟时序
p <- plot_cells_2(cds = Tumor_monocle3_cds,
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
pdf(paste0("./轨迹比较/Tumor_Monocle3_Step2_",
           reduction_method, "_pseudotime", ".pdf"),
    width = 6.8, height = 5)
print(p)
dev.off()

saveRDS(Tumor_monocle3_cds,
        "3.轨迹比较.Tumor_monocle3_cds.rds")






