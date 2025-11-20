Osteoblasts_sub <- subset(seurat_out,
                          subset = Osteoblasts_sub %in% c("Osteoblasts (PTCH1+ cells)",
                                                          "Osteoblasts (Core)"))
DefaultAssay(Osteoblasts_sub) <- "Spatial"
Osteoblasts_sub <- SCTransform(Osteoblasts_sub, assay = "Spatial",
                               verbose = TRUE)
Osteoblasts_sub <- RunPCA(Osteoblasts_sub, assay = "SCT", verbose = TRUE)
Osteoblasts_sub <- FindNeighbors(Osteoblasts_sub, reduction = "pca", dims = 1:30)
Osteoblasts_sub <- FindClusters(Osteoblasts_sub, resolution = 0.17)
Osteoblasts_sub <- RunUMAP(Osteoblasts_sub, reduction = "pca", dims = 1:30,
                           min.dist = 0.6)
pdf("../骨轨迹/Osteoblasts_sub_umap.pdf",
    width = 6.5, height = 4)
DimPlot(Osteoblasts_sub, group.by = "Osteoblasts_sub",
        cols = Celltype_color_3) + ggtitle("")
dev.off()
