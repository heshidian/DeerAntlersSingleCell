
RNA_singlet_EC_1 <- subset(RNA_singlet_EC,
                           subset = seurat_clusters %in% c(0, 1))
DimPlot(RNA_singlet_EC_1)

DimPlot(RNA_singlet_EC_1, label = TRUE, group.by = "seurat_clusters")
DimPlot(RNA_singlet_EC_1, label = TRUE, group.by = "Phase")
DimPlot(RNA_singlet_EC_1, label = TRUE)
DefaultAssay(RNA_singlet_EC_1) <- "RNA"
FeaturePlot(RNA_singlet_EC_1,
            features = c("MKI67"), order = TRUE)
VlnPlot(RNA_singlet_EC_1, features = "MKI67",
        group.by = "sample")
VlnPlot(RNA_singlet_EC_1, features = "CytoTrace2",
        group.by = "sample")

pdf("(EC)ATAC_3.CytoTRACE2.pdf", height = 4, width = 3.6)
violin_linebox(RNA_singlet_EC_1@meta.data, group.by = "sample",
               features = c("CytoTrace2"), linewidth = 2,
               group.color = RNA_group_color) +
  # 添加上方红色色块
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = 0.3, ymax = Inf, alpha = 0.15, fill = "red") +
  # 添加下方灰色色块
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = 0.3, alpha = 0.15, fill = "gray") +
  annotate("text", x = 0.5, y = 0.57, label = "EC-active", hjust = 0, vjust = 1, size = 5) +
  annotate("text", x = 0.5, y = -0.03, label = "EC-dormant", hjust = 0, vjust = 0, size = 5) +
  NoLegend() +
  labs(title = "", y = "CytoTRACE2 score") +
  geom_hline(yintercept = 0.3, color = "red",
             linetype = "dashed", linewidth = 1)
dev.off()




cell_prolife <- read.gmt("CELL_PROLIFERATION_GO_0008283.v2024.1.Hs.gmt")

cell_prolife_AUC <- AUCell_calcAUC(geneSets = list("Cell proliferation" = cell_prolife$gene),
                                   rankings = ranking,
                                   nCores = 20)
RNA_singlet_Mural$CELL_PROLIFERATION <- cell_prolife_AUC@assays@data$AUC[,colnames(RNA_singlet_Mural)]

RNA_singlet_Mural_1 <- subset(RNA_singlet_Mural,
                              subset = seurat_clusters == 1)

DimPlot(RNA_singlet_Mural_1, label = TRUE, group.by = "seurat_clusters")
DimPlot(RNA_singlet_Mural_1, label = TRUE, group.by = "Phase")
DimPlot(RNA_singlet_Mural_1, label = TRUE)
DefaultAssay(RNA_singlet_Mural_1) <- "RNA"
FeaturePlot(RNA_singlet_Mural_1,
            features = c("MKI67"), order = TRUE)
VlnPlot(RNA_singlet_Mural_1, features = "MKI67",
        group.by = "sample")
VlnPlot(RNA_singlet_Mural_1, features = "CytoTrace2",
        group.by = "sample")

pdf("ATAC_3.CytoTRACE2.pdf", height = 4, width = 3.6)
violin_linebox(RNA_singlet_Mural_1@meta.data, group.by = "sample",
               features = c("CytoTrace2"), linewidth = 2,
               group.color = RNA_group_color) +
  # 添加上方红色色块
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = 0.3, ymax = Inf, alpha = 0.15, fill = "red") +
  # 添加下方灰色色块
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = 0.3, alpha = 0.15, fill = "gray") +
  annotate("text", x = 0.5, y = 0.57, label = "Mural-active", hjust = 0, vjust = 1, size = 5) +
  annotate("text", x = 0.5, y = -0.03, label = "Mural-dormant", hjust = 0, vjust = 0, size = 5) +
  NoLegend() +
  labs(title = "", y = "CytoTRACE2 score") +
  geom_hline(yintercept = 0.3, color = "red",
             linetype = "dashed", linewidth = 1)
dev.off()

