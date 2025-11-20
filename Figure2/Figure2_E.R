##################


ATAC_Mural_1 <- addBgdPeaks(ATAC_Mural_1, force = TRUE)
ATAC_Mural_1 <- addMotifAnnotations(ArchRProj = ATAC_Mural_1,
                                    motifPWMs = Motif_pwm,
                                    name = "Motif",
                                    cutOff = 1e-05,
                                    force = TRUE)
ATAC_Mural_1 <- addDeviationsMatrix(
  ArchRProj = ATAC_Mural_1,
  peakAnnotation = "Motif",
  force = TRUE
)

markersGS <- getMarkerFeatures(
  ArchRProj = ATAC_Mural_1, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "predictedGroup_Un",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersMotif <- getMarkerFeatures(
  ArchRProj = ATAC_Mural_1, 
  useMatrix = "MotifMatrix", 
  groupBy = "predictedGroup_Un",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersExp <- getMarkerFeatures(
  ArchRProj = ATAC_Mural_1, 
  useMatrix = "GeneIntegrationMatrix_2", 
  groupBy = "predictedGroup_Un",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markerListGS <- getMarkers(markersGS,
                           cutOff = "FDR <= 1 & Log2FC >= -Inf")
markerListGS <- as.data.frame(markerListGS@listData$`Mural-active`)
rownames(markerListGS) <- markerListGS$name

markerListMotif <- getMarkers_custom(markersMotif,
                                     cutOff = "FDR <= 1 & AUC >= -Inf")
markerListMotif <- as.data.frame(markerListMotif@listData$`Mural-active`)
rownames(markerListMotif) <- markerListMotif$name

markerListExp <- getMarkers(markersExp,
                            cutOff = "FDR <= 1 & Log2FC >= -Inf")
markerListExp <- as.data.frame(markerListExp@listData$`Mural-active`)
rownames(markerListExp) <- markerListExp$name

Idents(RNA_singlet_Mural_1) <- RNA_singlet_Mural_1$CytoTrace2_group
markerListExp_RNA <- FindMarkers(RNA_singlet_Mural_1, logfc.threshold = -Inf,
                                 min.pct = 0.05,
                                 ident.1 = "Mural-active", ident.2 = "Mural-dormant")
markerListExp_RNA$Log2FC <- markerListExp_RNA$avg_log2FC

TFs <- read.table("/media/heshidian/RAID5_42TB/1.ReferenceFiles/pySCENIC/Human/allTFs_hg38.txt",
                  header = FALSE, sep = "\t")
TFs <- TFs[TFs$V1 %in% rownames(RNA_singlet_Mural_1@assays$RNA@data),]

GS_Exp <- data.frame(Gene = TFs,
                     GS_Log2FC = markerListGS[TFs, "Log2FC"],
                     # Exp_Log2FC = markerListExp[rownames(markerListGS), "Log2FC"],
                     Exp_Log2FC = markerListExp[TFs, "Log2FC"]
)
GS_Exp <- as.data.frame(na.omit(GS_Exp))

r <- c()
for (i in GS_Exp$Gene) {
  r <- c(r,
         cor(RNA_singlet_Mural_1@assays$RNA@data[i,],
             RNA_singlet_Mural_1$CytoTrace2,
             method = "pearson"))
}
GS_Exp$r <- r

ATAC_motif <- .getFeatureDF(getArrowFiles(ATAC_singlet),
                            "MotifMatrix", threads = 300)

GS_Exp <- GS_Exp[GS_Exp$Gene %in% ATAC_motif$name,]
ggplot(data = GS_Exp,
       aes(x = GS_Log2FC, y = Exp_Log2FC, color = r),
       size = 3) +
  geom_point() +
  scale_color_gradientn(colours = c("blue", "white", "red"),
                        limit = c(-0.2, 0.2),
                        oob = scales::squish
  ) +
  geom_vline(xintercept = 0.25) +
  geom_vline(xintercept = -0.25) +
  geom_hline(yintercept = 0.25) +
  geom_hline(yintercept = -0.25)

GS_Exp_active <- GS_Exp[GS_Exp$GS_Log2FC > 0 & GS_Exp$Exp_Log2FC > 0,]
GS_Exp_active$Group <- "Mural-active"
GS_Exp_dormant <- GS_Exp[GS_Exp$GS_Log2FC < -0 & GS_Exp$Exp_Log2FC < -0,]
GS_Exp_dormant$Group <- "Mural-dormant"

GS_Exp_active <- GS_Exp_active[order(GS_Exp_active$r,
                                     decreasing = TRUE),]
active_top5 <- GS_Exp_active[1:5,]
GS_Exp_dormant <- GS_Exp_dormant[order(GS_Exp_dormant$r,
                                       decreasing = FALSE),]
dormant_top5 <- GS_Exp_dormant[1:5,]

GS_Exp <- GS_Exp[order(abs(GS_Exp$r),
                       decreasing = FALSE),]
pdf("ATAC_6.positive_TFs_V2.pdf",
    width = 6.4, height = 4)
ggplot() +
  geom_point(data = GS_Exp,
             aes(x = Exp_Log2FC, y = GS_Log2FC, color = r),
             size = 3, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size = 12, color = "black",
                                  face = "bold", hjust = 0.5)) +
  # coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0.3, 0.7)) +
  labs(x = expression(log[2]*"FC of TF (RNA expression)"),
       y = expression(log[2]*"FC of TF (ATAC gene activity)"),
       title = "Mural-active vs. Mural-dormant",
       color = "Pearson correlation\n(Expression and cytoTrace score)") +
  scale_color_gradientn(colours = c("blue", "white", "red"),
                        limit = c(-0.2, 0.2),
                        oob = scales::squish
  ) +
  geom_vline(xintercept = 0, colour = "#696969", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "#696969", linetype = "dashed") +
  ggrepel::geom_text_repel(data = dormant_top5, aes(x = Exp_Log2FC, y = GS_Log2FC, 
                                                    label = Gene),
                           size = 3.3, fontface = "italic",
                           box.padding = unit(1, "lines"),
                           point.padding = unit(0, "lines"), 
                           min.segment.length = 0,
                           segment.color = "black",
                           colour="#000000",
                           show.legend = FALSE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 10)) +
  ggrepel::geom_text_repel(data = active_top5, aes(x = Exp_Log2FC, y = GS_Log2FC, 
                                                   label = Gene),
                           size = 3.3, fontface = "italic",
                           box.padding = unit(1.8, "lines"),
                           point.padding = unit(0, "lines"), 
                           min.segment.length = 0,
                           segment.color = "black",
                           colour="#000000",
                           show.legend = FALSE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 10))
dev.off()




###############


ATAC_EC_1 <- addBgdPeaks(ATAC_EC_1, force = TRUE)
ATAC_EC_1 <- addMotifAnnotations(ArchRProj = ATAC_EC_1,
                                 motifPWMs = Motif_pwm,
                                 name = "Motif",
                                 cutOff = 1e-05,
                                 force = TRUE)
ATAC_EC_1 <- addDeviationsMatrix(
  ArchRProj = ATAC_EC_1,
  peakAnnotation = "Motif",
  force = TRUE
)

markersGS <- getMarkerFeatures(
  ArchRProj = ATAC_EC_1, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "predictedGroup_Un",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersMotif <- getMarkerFeatures(
  ArchRProj = ATAC_EC_1, 
  useMatrix = "MotifMatrix", 
  groupBy = "predictedGroup_Un",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersExp <- getMarkerFeatures(
  ArchRProj = ATAC_EC_1, 
  useMatrix = "GeneIntegrationMatrix_2", 
  groupBy = "predictedGroup_Un",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markerListGS <- getMarkers(markersGS,
                           cutOff = "FDR <= 1 & Log2FC >= -Inf")
markerListGS <- as.data.frame(markerListGS@listData$`EC-active`)
rownames(markerListGS) <- markerListGS$name

markerListMotif <- getMarkers_custom(markersMotif,
                                     cutOff = "FDR <= 1 & AUC >= -Inf")
markerListMotif <- as.data.frame(markerListMotif@listData$`EC-active`)
rownames(markerListMotif) <- markerListMotif$name

markerListExp <- getMarkers(markersExp,
                            cutOff = "FDR <= 1 & Log2FC >= -Inf")
markerListExp <- as.data.frame(markerListExp@listData$`EC-active`)
rownames(markerListExp) <- markerListExp$name

markerListExp_RNA <- FindMarkers(RNA_singlet_EC_1, logfc.threshold = -Inf,
                                 min.pct = 0.05,
                                 ident.1 = "EC-active", ident.2 = "EC-dormant")
markerListExp_RNA$Log2FC <- markerListExp_RNA$avg_log2FC

TFs <- read.table("/media/heshidian/RAID5_42TB/1.ReferenceFiles/pySCENIC/Human/allTFs_hg38.txt",
                  header = FALSE, sep = "\t")
TFs <- TFs[TFs$V1 %in% rownames(RNA_singlet_EC_1@assays$RNA@data),]

GS_Exp <- data.frame(Gene = TFs,
                     GS_Log2FC = markerListGS[TFs, "Log2FC"],
                     # Exp_Log2FC = markerListExp[rownames(markerListGS), "Log2FC"],
                     Exp_Log2FC = markerListExp[TFs, "Log2FC"]
)
GS_Exp <- as.data.frame(na.omit(GS_Exp))

r <- c()
for (i in GS_Exp$Gene) {
  r <- c(r,
         cor(RNA_singlet_EC_1@assays$RNA@data[i,],
             RNA_singlet_EC_1$CytoTrace2,
             method = "spearman"))
}
GS_Exp$r <- r

ATAC_motif <- .getFeatureDF(getArrowFiles(ATAC_singlet),
                            "MotifMatrix", threads = 300)

GS_Exp <- GS_Exp[GS_Exp$Gene %in% ATAC_motif$name,]
ggplot(data = GS_Exp,
       aes(x = GS_Log2FC, y = Exp_Log2FC, color = r),
       size = 3) +
  geom_point() +
  scale_color_gradientn(colours = c("blue", "white", "red"),
                        limit = c(-0.2, 0.2),
                        oob = scales::squish
  ) +
  geom_vline(xintercept = 0.25) +
  geom_vline(xintercept = -0.25) +
  geom_hline(yintercept = 0.25) +
  geom_hline(yintercept = -0.25)

GS_Exp_active <- GS_Exp[GS_Exp$GS_Log2FC > 0 & GS_Exp$Exp_Log2FC > 0,]
GS_Exp_active$Group <- "EC-active"
GS_Exp_dormant <- GS_Exp[GS_Exp$GS_Log2FC < -0 & GS_Exp$Exp_Log2FC < -0,]
GS_Exp_dormant$Group <- "EC-dormant"
GS_Exp_active <- GS_Exp_active[order(GS_Exp_active$r,
                                     decreasing = TRUE),]
active_top5 <- GS_Exp_active[1:5,]
GS_Exp_dormant <- GS_Exp_dormant[order(GS_Exp_dormant$r,
                                       decreasing = FALSE),]
dormant_top5 <- GS_Exp_dormant[1:5,]

GS_Exp <- GS_Exp[order(abs(GS_Exp$r),
                       decreasing = FALSE),]

library(ggfx)
pdf("(EC)ATAC_6.positive_TFs_V2.pdf",
    width = 6.4, height = 4)
ggplot() +
  geom_point(data = GS_Exp,
             aes(x = Exp_Log2FC, y = GS_Log2FC, color = r),
             size = 3, alpha = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size = 12, color = "black",
                                  face = "bold", hjust = 0.5)) +
  # coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0.3, 0.7)) +
  labs(x = expression(log[2]*"FC of TF (RNA expression)"),
       y = expression(log[2]*"FC of TF (ATAC gene activity)"),
       title = "EC-active vs. EC-dormant",
       color = "Pearson correlation\n(Expression and cytoTrace score)") +
  scale_color_gradientn(colours = c("blue", "white", "red"),
                        limit = c(-0.2, 0.2),
                        oob = scales::squish
  ) +
  geom_vline(xintercept = 0, colour = "#696969", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "#696969", linetype = "dashed") +
  ggrepel::geom_text_repel(data = dormant_top5, aes(x = Exp_Log2FC, y = GS_Log2FC, 
                                                    label = Gene),
                           size = 3.3, fontface = "italic",
                           box.padding = unit(1, "lines"),
                           point.padding = unit(0, "lines"), 
                           min.segment.length = 0,
                           segment.color = "black",
                           colour="#000000",
                           show.legend = FALSE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 10)) +
  ggrepel::geom_text_repel(data = active_top5, aes(x = Exp_Log2FC, y = GS_Log2FC, 
                                                   label = Gene),
                           size = 3.3, fontface = "italic",
                           box.padding = unit(1.8, "lines"),
                           point.padding = unit(0, "lines"), 
                           min.segment.length = 0,
                           segment.color = "black",
                           colour="#000000",
                           show.legend = FALSE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 10))
dev.off()
