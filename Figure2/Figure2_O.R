
monocle3_pseudotime <- pseudotime(monocle3_cds)
monocle3_pseudotime <- data.frame(Pseudotime = monocle3_pseudotime)

sum(rownames(monocle3_pseudotime) %in% colnames(seurat_out_os_PTCH1_3))
DefaultAssay(seurat_out_os_PTCH1_3) <- "Spatial"
seurat_out_os_PTCH1_3 <- NormalizeData(seurat_out_os_PTCH1_3)

monocle3_pseudotime$COL10A1 <- seurat_out_os_PTCH1_3@assays$Spatial@data["COL10A1",
                                                                         rownames(monocle3_pseudotime)]
monocle3_pseudotime$COL11A2 <- seurat_out_os_PTCH1_3@assays$Spatial@data["COL11A2",
                                                                         rownames(monocle3_pseudotime)]
monocle3_pseudotime$COL1A2 <- seurat_out_os_PTCH1_3@assays$Spatial@data["COL1A2",
                                                                        rownames(monocle3_pseudotime)]
monocle3_pseudotime$IBSP <- seurat_out_os_PTCH1_3@assays$Spatial@data["IBSP",
                                                                      rownames(monocle3_pseudotime)]
monocle3_pseudotime$MMP13 <- seurat_out_os_PTCH1_3@assays$Spatial@data["MMP13",
                                                                       rownames(monocle3_pseudotime)]
monocle3_pseudotime$MGP <- seurat_out_os_PTCH1_3@assays$Spatial@data["MGP",
                                                                     rownames(monocle3_pseudotime)]
monocle3_pseudotime$COL2A1 <- seurat_out_os_PTCH1_3@assays$Spatial@data["COL2A1",
                                                                        rownames(monocle3_pseudotime)]
monocle3_pseudotime$PHEX <- seurat_out_os_PTCH1_3@assays$Spatial@data["PHEX",
                                                                      rownames(monocle3_pseudotime)]
monocle3_pseudotime$COMP <- seurat_out_os_PTCH1_3@assays$Spatial@data["COMP",
                                                                      rownames(monocle3_pseudotime)]
monocle3_pseudotime$BGLAP <- seurat_out_os_PTCH1_3@assays$Spatial@data["BGLAP",
                                                                       rownames(monocle3_pseudotime)]
monocle3_pseudotime$OMD <- seurat_out_os_PTCH1_3@assays$Spatial@data["OMD",
                                                                     rownames(monocle3_pseudotime)]
monocle3_pseudotime$VCAN <- seurat_out_os_PTCH1_3@assays$Spatial@data["VCAN",
                                                                      rownames(monocle3_pseudotime)]
monocle3_pseudotime$ASPN <- seurat_out_os_PTCH1_3@assays$Spatial@data["ASPN",
                                                                      rownames(monocle3_pseudotime)]
monocle3_pseudotime$CHAD <- seurat_out_os_PTCH1_3@assays$Spatial@data["CHAD",
                                                                      rownames(monocle3_pseudotime)]
monocle3_pseudotime$SP7 <- seurat_out_os_PTCH1_3@assays$Spatial@data["SP7",
                                                                     rownames(monocle3_pseudotime)]
monocle3_pseudotime$SPP1 <- seurat_out_os_PTCH1_3@assays$Spatial@data["SPP1",
                                                                      rownames(monocle3_pseudotime)]

# for (i in 2:ncol(monocle3_pseudotime)) {
#   monocle3_pseudotime[,i] <- (monocle3_pseudotime[,i] - min(monocle3_pseudotime[,i])) / (max(monocle3_pseudotime[,i]) - min(monocle3_pseudotime[,i]))
# }
monocle3_pseudotime_df <- reshape::melt.data.frame(monocle3_pseudotime,
                                                   id.vars = "Pseudotime")

monocle3_pseudotime_df <- monocle3_pseudotime_df[monocle3_pseudotime_df$variable %in% c("COL10A1",
                                                                                        "IBSP",
                                                                                        "BGLAP",
                                                                                        "SP7",
                                                                                        "SPP1"),]
pdf("../骨轨迹/Osteoblasts_PTCH1_traj_pseudo_gene_curve.pdf",
    width = 5, height = 2.5)
ggplot() +
  geom_smooth(data = monocle3_pseudotime_df,
              aes(x = Pseudotime, y = value, color = variable),
              method = "loess", se = FALSE) +
  geom_segment(aes(x = min(monocle3_pseudotime_df$Pseudotime),
                   xend = max(monocle3_pseudotime_df$Pseudotime),
                   y = 0, yend = 0),
               arrow = arrow(length = unit(0.2, "cm"),
                             ends = "last", type = "open"),
               color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10, colour = "black")) +
  labs(color = "Gene", x = "Pseudotime", y = "Normalized expression level")
dev.off()
