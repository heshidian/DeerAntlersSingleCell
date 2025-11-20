ATAC_singlet <- readRDS("ATAC_singlet.rds")
ATAC_Singlet_GS <- getMatrixFromProject(ATAC_singlet, useMatrix = "GeneScoreMatrix")
ATAC_Singlet_GS@assays@data@listData[["GeneScoreMatrix"]]@Dimnames[[1]] <- ATAC_Singlet_GS@elementMetadata@listData[["name"]]
ATAC_Singlet_GS_log <- log2(ATAC_Singlet_GS@assays@data@listData[["GeneScoreMatrix"]] + 1)
ATAC_meta <- getCellColData(ATAC_singlet)
ATAC_meta <- as.data.frame(ATAC_meta)

cells <- rownames(ATAC_meta)[ATAC_meta$Celltype %in% c("AnMCs",
                                                       "Progenitor cells",
                                                       "Proliferative_progenitor cells",
                                                       "Chondrocytes",
                                                       "Hypertrophic chondrocytes")]
ATAC_Singlet_GS_log <- ATAC_Singlet_GS_log[,cells]

p <- DotPlot_2(Exp = ATAC_Singlet_GS_log,
               group = as.character(ATAC_meta[colnames(ATAC_Singlet_GS_log),
                                              "Celltype"]),
               genes = c("TP53", "MYC"),
               genes_order = c("TP53", "MYC"),
               group_order = unique(as.character(ATAC_meta[colnames(ATAC_Singlet_GS_log),
                                                           "Celltype"])),
               fill_label = "Gene Activity", size_label = "Percent",
               high_color = "red", low_color = "blue")
p <- egg::set_panel_size(p, width = unit(3, "cm"), height = unit(5, "cm"))
dev.off()
pdf("imeta_反修_ATAC_MYC_TP53.pdf", width = 7, height = 7)
grid.draw(p)
dev.off()