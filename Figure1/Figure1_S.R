
PC_targets <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/打分补充/PC_MYC_targets_peaks.txt",
                       header = FALSE)
common_genes <- intersect(unique(PC_targets[,1]),
                          unique(MYC_peaks$nearestGene))
common_genes <- MYC_peaks[MYC_peaks$nearestGene %in% common_genes,]
p <- ggplot() +
  geom_point(data = MYC_peaks,
             aes(x = peak_gene_cor, y = MYC_target_cor),
             size = 1.5, color = "gray") +
  geom_point(data = MYC_peaks_target,
             aes(x = peak_gene_cor, y = MYC_target_cor),
             size = 2.5, color = "green") +
  geom_point(data = common_genes,
             aes(x = peak_gene_cor, y = MYC_target_cor),
             size = 2, color = "black",
             shape = 21, stroke = 1#,
             # color = "#F19695"
  ) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  geom_hline(yintercept = 0.3, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, colour = "black")) +
  labs(x = "Pearson correlation coefficient between\npeak accessibility and target expression",
       y = "Pearson correlation coefficient between\nMYC expression and target expression")
p <- egg::set_panel_size(p, width = unit(8, "cm"),
                         height = unit(8, "cm"))
dev.off()
pdf("./Tumor_MYC_targets_PC_peaks.pdf", width = 8, height = 8)
grid::grid.draw(p)
dev.off()



Tumor_peaks <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/Public_tumor_data/GSE215758/Tumor_MYC_targets_peaks.txt",
                        header = FALSE, na.strings = "")
common_genes <- intersect(unique(Tumor_peaks[,1]),
                          unique(MYC_peaks$nearestGene))
common_genes <- MYC_peaks[MYC_peaks$nearestGene %in% common_genes,]
p <- ggplot() +
  geom_point(data = MYC_peaks,
             aes(x = peak_gene_cor, y = MYC_target_cor),
             size = 1.5, color = "gray") +
  geom_point(data = MYC_peaks_target,
             aes(x = peak_gene_cor, y = MYC_target_cor),
             size = 2.5, color = "red") +
  geom_point(data = common_genes,
             aes(x = peak_gene_cor, y = MYC_target_cor),
             size = 2, color = "black",
             shape = 21, stroke = 1#,
             # color = "#F19695"
  ) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  geom_hline(yintercept = 0.3, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10, colour = "black")) +
  labs(x = "Pearson correlation coefficient between\npeak accessibility and target expression",
       y = "Pearson correlation coefficient between\nMYC expression and target expression")
p <- egg::set_panel_size(p, width = unit(8, "cm"),
                         height = unit(8, "cm"))
dev.off()
pdf("./打分补充/PC_MYC_targets_Tumor_peaks.pdf", width = 8, height = 8)
grid::grid.draw(p)
dev.off()
