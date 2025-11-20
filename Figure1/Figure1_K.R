RNA_singlet <- readRDS("RNA_singlet.rds")
temp <- subset(RNA_singlet,
               subset = Celltype_rename %in% c("Progenitor cells",
                                               "Proliferative_progenitor cells"))

temp$sample <- unlist(
  lapply(as.character(temp$sample),
         function(x){
           switch (x,
                   "1.RM" = "RM",
                   "2.PC" = "PC",
                   "3.TZ" = "TZ",
                   "4.CA" = "CA",
                   "5.MC" = "MC"
           )
         })
)
temp$sample <- factor(temp$sample,
                      levels = c("RM", "PC", "TZ", "CA", "MC"))
cytotrace2_result <- cytotrace2(temp,   
                                species = "human",
                                is_seurat = TRUE,
                                slot_type = "counts",
                                batch_size = 10000,
                                smooth_batch_size = 1000,
                                parallelize_models = TRUE,
                                parallelize_smoothing = TRUE,
                                ncores = 10,
                                seed = 14) 

temp$CytoTrace2 <- cytotrace2_result@meta.data[colnames(temp), "CytoTRACE2_Score"]

RNA_group_color <- c("#E7211A", "#EFEA3C", "#72C8D5", "#6AB82D", "#18499E")
names(RNA_group_color) <- c("RM", "PC", "TZ", "CA", "MC")

pdf("CytoTRACE2_PC_ProlifPC.pdf",
    height = 5, width = 2)
violin_linebox(temp@meta.data, group.by = "Celltype_rename",
               features = c("CytoTrace2"), linewidth = 2,
               group.color = RNA_Celltype_color) +
  # 添加上方红色色块
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = 0.3, ymax = Inf, alpha = 0.15, fill = "red") +
  # 添加下方灰色色块
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = 0.3, alpha = 0.15, fill = "gray") +
  annotate("text", x = 0.5, y = 0.7, label = "active", hjust = 0, vjust = 1, size = 5) +
  annotate("text", x = 0.5, y = 0, label = "dormant", hjust = 0, vjust = 0, size = 5) +
  NoLegend() +
  labs(title = "", y = "CytoTRACE2 score") +
  geom_hline(yintercept = 0.3, color = "red",
             linetype = "dashed", linewidth = 1)
dev.off()
