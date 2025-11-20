
plot1 <- openxlsx::read.xlsx("AnPC_vs_Tumor_Chondroblasts_Osteoblastic_selected.xlsx",
                             sheet = "Tumor_Osteoblastic")
# plot1 <- plot1[plot1$P < 0.05,]
plot1$pathway <- factor(plot1$pathway,
                        levels = rev(plot1$pathway))
plot1$log10P <- -log10(plot1$P)
plot1$log10P[plot1$log10P == Inf] <- max(plot1$log10P[plot1$log10P != Inf])

plot1$fdr <- p.adjust(plot1$P, method = "fdr")

pdf("PCs_vs_Osteoblastic.pdf", width = 7, height = 6)
p <- ggplot() +
  geom_bar(data = plot1, aes(x = log2FC, y = pathway),
           stat = "identity", width = 0.1) +
  geom_point(data = plot1, aes(x = log2FC, y = pathway,
                               size = abs(log2FC),
                               color = log10P)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_gradientn(limits = c(0, 100),
                        colors = c("#006400", "#FFD700", "#F4A460", "#FA8072"),
                        # colors = viridis(n = 100, option = "D"),
                        oob = scales::squish) +
  scale_size_continuous(
    range = c(1, 4),         # 设置点的最小/最大半径
    labels = c(0.1, 0.4, 0.8, 1.2),     # 图例中显示的参考值
    breaks = c(0.1, 0.4, 0.8, 1.2),
    name = "abs(log2FC)"     # 图例标题
  ) +
  labs(x = "log2(Fold Change)", y = "",
       size = "Abs(log2FC)", color = "-log10(P-value)") +
  ggtitle("PCs vs. Osteoblastic")
p <- egg::set_panel_size(p,
                         width = unit(5, "cm"),
                         height = unit(9, "cm"))
grid::grid.draw(p)
dev.off()

