smoothMat_PK <- getTrajectory_custom(ATAC_Chondrocyte,
                                     name = "Pseudotime",
                                     break_n = 60,
                                     smoothWindow = 10,
                                     log2Norm = TRUE,
                                     scaleTo = 10000,
                                     useMatrix = "PeakMatrix")
peaks <- as.data.frame(smoothMat_PK[["seTrajectory"]]@elementMetadata@listData)
peaks <- paste0(peaks$seqnames, "_", peaks$start, "_", peaks$end)
smoothMat_PK_matrix <- smoothMat_PK[["seTrajectory"]]@assays@data@listData[["smoothMat"]]
rownames(smoothMat_PK_matrix) <- peaks

smoothMat_MM <- getTrajectory_custom(ATAC_Chondrocyte,
                                     name = "Pseudotime",
                                     break_n = 60, 
                                     smoothWindow = 10,
                                     log2Norm = FALSE, # 因为有负值
                                     scaleTo = NULL, # 因为有负值
                                     useMatrix = "MotifMatrix")
rownames(smoothMat_MM[["seTrajectory"]]@assays@data@listData[["smoothMat"]]) <- smoothMat_MM[["seTrajectory"]]@NAMES
smoothMat_MM_matrix <- smoothMat_MM[["seTrajectory"]]@assays@data@listData[["smoothMat"]]
rownames(smoothMat_MM_matrix) <- unlist(lapply(rownames(smoothMat_MM_matrix), function(x){
  unlist(strsplit(x, ":", fixed = TRUE))[2]
}))

smoothMat_GS <- getTrajectory_custom(ATAC_Chondrocyte,
                                     name = "Pseudotime",
                                     break_n = 60,
                                     smoothWindow = 10,
                                     log2Norm = TRUE,
                                     scaleTo = 10000,
                                     useMatrix = "GeneScoreMatrix")
smoothMat_GS_matrix <- smoothMat_GS[["seTrajectory"]]@assays@data@listData[["smoothMat"]]
rownames(smoothMat_GS_matrix) <- smoothMat_GS[["seTrajectory"]]@NAMES
rownames(smoothMat_GS_matrix) <- unlist(lapply(rownames(smoothMat_GS_matrix), function(x){
  unlist(strsplit(x, ":", fixed = TRUE))[2]
}))

smoothMat_Exp <- getTrajectory_custom(ATAC_Chondrocyte,
                                      name = "Pseudotime",
                                      break_n = 60,
                                      smoothWindow = 10,
                                      log2Norm = TRUE,
                                      scaleTo = 10000,
                                      useMatrix = "GeneIntegrationMatrix_Co_harmony_Annotation_Group")
smoothMat_Exp_matrix <- smoothMat_Exp[["seTrajectory"]]@assays@data@listData[["smoothMat"]]
rownames(smoothMat_Exp_matrix) <- smoothMat_Exp[["seTrajectory"]]@NAMES
rownames(smoothMat_Exp_matrix) <- unlist(lapply(rownames(smoothMat_Exp_matrix), function(x){
  unlist(strsplit(x, ":", fixed = TRUE))[2]
}))

RNA_monocle3_dynamic_genes <- readRDS("RNA_monocle3_dynamic_genes.rds")
RNA_monocle3_dynamic_genes_cluster <- readRDS("RNA_monocle3_dynamic_genes_cluster.rds")

sum(rownames(RNA_monocle3_dynamic_genes) %in% rownames(smoothMat_GS_matrix))
sum(rownames(RNA_monocle3_dynamic_genes) %in% rownames(smoothMat_Exp_matrix))
nrow(RNA_monocle3_dynamic_genes)

smoothMat_Exp_matrix_RNA <- smoothMat_Exp_matrix[rownames(RNA_monocle3_dynamic_genes),]
smoothMat_Exp_matrix_RNA <- t(smoothMat_Exp_matrix_RNA)
for (i in 1:ncol(smoothMat_Exp_matrix_RNA)) {
  smoothMat_Exp_matrix_RNA[,i] <- (smoothMat_Exp_matrix_RNA[,i] - mean(smoothMat_Exp_matrix_RNA[,i])) / sd(smoothMat_Exp_matrix_RNA[,i])
}
smoothMat_Exp_matrix_RNA <- t(smoothMat_Exp_matrix_RNA)
smoothMat_Exp_matrix_RNA[smoothMat_Exp_matrix_RNA > 2] <- 2
smoothMat_Exp_matrix_RNA[smoothMat_Exp_matrix_RNA < -2] <- -2
row_anno <- rowAnnotation(
  GeneCluster = as.character(RNA_monocle3_dynamic_genes_cluster[rownames(smoothMat_Exp_matrix_RNA), "Cluster"]), 
  col = list(GeneCluster = c("1" = "#008080", "2" = "#DA70D6",
                          "3" = "#778899", "4" = "#FB8072",
                          "5" = "#B0E0E6", "6" = "#FDB462",
                          "7" = "#00CED1", "8" = "#FFC0CB")),
  show_annotation_name = FALSE
)

# RNA_monocle3_dynamic_genes[RNA_monocle3_dynamic_genes > 2] <- 2
# RNA_monocle3_dynamic_genes[RNA_monocle3_dynamic_genes < -2] <- -2
h1 <- Heatmap(RNA_monocle3_dynamic_genes,
              show_column_names = FALSE, show_row_names = FALSE,
              cluster_rows = FALSE, cluster_columns = FALSE,
              col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
              row_names_side = "left",
              name = "RNA\nRow scaled gene expression",
              width = unit(2, "null"),
              heatmap_legend_param = list(
                # title_position = "topcenter",  # 图例标题居中
                legend_direction = "horizontal"  # 图例方向为水平
              ),
              left_annotation = row_anno
)
h1

h2 <- Heatmap(smoothMat_Exp_matrix_RNA,
        show_column_names = FALSE, show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,
        col = paletteContinuous(set = "beach", n = 100),
        row_names_side = "left",
        name = "ATAC\nRow scaled gene expression",
        width = unit(2, "null"),
        heatmap_legend_param = list(
          # title_position = "topcenter",  # 图例标题居中
          legend_direction = "horizontal"  # 图例方向为水平
        )#,
        # left_annotation = row_anno
)

pdf("ATAC_Monocle3_RNA_ATAC_genes_heatmap.pdf",
    height = 4.5, width = 6)
h1 + h2
dev.off()

peakset_df_pd_genes <- peakset_df_pd[peakset_df_pd$nearestGene %in% rownames(smoothMat_Exp_matrix_RNA),]
length(unique(peakset_df_pd_genes$nearestGene))
peakset_df_pd_genes$Gene_Cluster <- RNA_monocle3_dynamic_genes_cluster[peakset_df_pd_genes$nearestGene, "Cluster"]
peakset_df_pd_genes$Gene_Cluster <- as.character(peakset_df_pd_genes$Gene_Cluster)
length(unique(peakset_df_pd_genes$nearestGene))

Motif_peak_match <- getPeakAnnotation(ATAC_Chondrocyte, name = "Motif")
Motif_peak_match <- readRDS(Motif_peak_match[["Matches"]])
peakset_df <- Motif_peak_match@rowRanges
peakset_df@ranges@NAMES <- NULL
peakset_df <- as.data.frame(peakset_df)
rownames(peakset_df) <- paste0(peakset_df$seqnames, "_",
                               peakset_df$start, "_",
                               peakset_df$end)
Motif_peak_match <- Motif_peak_match@assays@data@listData[["matches"]]
rownames(Motif_peak_match) <- rownames(peakset_df)
Motif_peak_match <- as.data.frame(Motif_peak_match)
Motif_peak_match <- Motif_peak_match[colnames(Motif_peak_match) %in% rownames(smoothMat_Exp_matrix)]

TF_peak_gene_cor <- data.frame(Target = peakset_df_pd_genes$nearestGene,
                               Gene_Cluster = peakset_df_pd_genes$Gene_Cluster,
                               Peak = peakset_df_pd_genes$Peak_id
                               )
TF_peak_gene_cor_list <- c()
for (i in 1:nrow(TF_peak_gene_cor)) {
  print(i)
  temp <- Motif_peak_match[TF_peak_gene_cor$Peak[i],]
  tf <- colnames(temp)[temp[1,] == TRUE]
  if (length(tf) > 0) {
    temp <- data.frame(Target = TF_peak_gene_cor$Target[i],
                       Gene_Cluster = TF_peak_gene_cor$Gene_Cluster[i],
                       Peak = TF_peak_gene_cor$Peak[i],
                       TF = tf)
    TF_peak_gene_cor_list <- c(TF_peak_gene_cor_list,
                               list(temp))
  }
}
TF_peak_gene_cor <- data.table::rbindlist(TF_peak_gene_cor_list)
length(unique(TF_peak_gene_cor$Target))
length(unique(TF_peak_gene_cor$TF))

TF_motif_TF_exp_cor <- c()
TF_motif_Gene_exp_cor <- c()
TF_exp_Gene_exp_cor <- c()
TF_motif_TF_exp_cor_P <- c()
TF_motif_Gene_exp_cor_P <- c()
TF_exp_Gene_exp_cor_P <- c()

for (i in 1:nrow(TF_peak_gene_cor)) {
  print(i)
  tf <- TF_peak_gene_cor$TF[i]
  # peak <- TF_peak_gene_cor$Peak[i]
  gene <- TF_peak_gene_cor$Target[i]
  # TF_motif_TF_exp_cor
  test <- cor.test(smoothMat_MM_matrix[tf,],
                   smoothMat_Exp_matrix[tf,])
  TF_motif_TF_exp_cor <- c(TF_motif_TF_exp_cor,
                           test$estimate)
  TF_motif_TF_exp_cor_P <- c(TF_motif_TF_exp_cor_P,
                           test$p.value)
  # TF_motif_Gene_exp_cor
  test <- cor.test(smoothMat_MM_matrix[tf,],
                   smoothMat_Exp_matrix[gene,])
  TF_motif_Gene_exp_cor <- c(TF_motif_Gene_exp_cor,
                             test$estimate)
  TF_motif_Gene_exp_cor_P <- c(TF_motif_Gene_exp_cor_P,
                               test$p.value)
  # TF_exp_Gene_exp_cor
  test <- cor.test(smoothMat_Exp_matrix[tf,],
                   smoothMat_Exp_matrix[gene,])
  TF_exp_Gene_exp_cor <- c(TF_exp_Gene_exp_cor,
                         test$estimate)
  TF_exp_Gene_exp_cor_P <- c(TF_exp_Gene_exp_cor_P,
                           test$p.value)
}

TF_peak_gene_cor$TF_motif_TF_exp_cor <- TF_motif_TF_exp_cor
TF_peak_gene_cor$TF_motif_TF_exp_cor_P <- TF_motif_TF_exp_cor_P
TF_peak_gene_cor$TF_motif_Gene_exp_cor <- TF_motif_Gene_exp_cor
TF_peak_gene_cor$TF_motif_Gene_exp_cor_P <- TF_motif_Gene_exp_cor_P
TF_peak_gene_cor$TF_exp_Gene_exp_cor <- TF_exp_Gene_exp_cor
TF_peak_gene_cor$TF_exp_Gene_exp_cor_P <- TF_exp_Gene_exp_cor_P

cor_threshold <- 0.6
TF_peak_gene_cor$select_positive <- (TF_peak_gene_cor$TF_motif_TF_exp_cor > cor_threshold & TF_peak_gene_cor$TF_motif_TF_exp_cor_P < 0.05) & # 有正有负
  (TF_peak_gene_cor$TF_motif_Gene_exp_cor > cor_threshold & TF_peak_gene_cor$TF_motif_Gene_exp_cor_P < 0.05) & # 有正有负
  (TF_peak_gene_cor$TF_exp_Gene_exp_cor > cor_threshold & TF_peak_gene_cor$TF_exp_Gene_exp_cor_P < 0.05)
table(TF_peak_gene_cor$select_positive)

TF_peak_gene_cor$select_negative <- (TF_peak_gene_cor$TF_motif_TF_exp_cor > cor_threshold & TF_peak_gene_cor$TF_motif_TF_exp_cor_P < 0.05) & # 有正有负
  (TF_peak_gene_cor$TF_motif_Gene_exp_cor < -cor_threshold & TF_peak_gene_cor$TF_motif_Gene_exp_cor_P < 0.05) & # 有正有负
  (TF_peak_gene_cor$TF_exp_Gene_exp_cor < -cor_threshold & TF_peak_gene_cor$TF_exp_Gene_exp_cor_P < 0.05)
table(TF_peak_gene_cor$select_negative)

TF_peak_gene_cor$select <- TF_peak_gene_cor$select_positive + TF_peak_gene_cor$select_negative
table(TF_peak_gene_cor$select)

TF_peak_gene_cor_select <- TF_peak_gene_cor[TF_peak_gene_cor$select == 1,]
length(unique(TF_peak_gene_cor_select$Target))
length(unique(TF_peak_gene_cor_select$TF))
table(TF_peak_gene_cor_select$TF, TF_peak_gene_cor_select$Gene_Cluster)

TF_peak_gene_cor_select_positive <- TF_peak_gene_cor_select[TF_peak_gene_cor_select$select_positive == 1,]
TF_peak_gene_cor_select_negative <- TF_peak_gene_cor_select[TF_peak_gene_cor_select$select_negative == 1,]
openxlsx::write.xlsx(TF_peak_gene_cor_select_positive,
                     "TF_gene_positive.xlsx")
openxlsx::write.xlsx(TF_peak_gene_cor_select_negative,
                     "TF_gene_negative.xlsx")

TF_cluster_gene_enrichment_p_positive <- matrix(1, nrow = length(unique(TF_peak_gene_cor_select$TF)), ncol = 8)
rownames(TF_cluster_gene_enrichment_p_positive) <- unique(TF_peak_gene_cor_select$TF)
colnames(TF_cluster_gene_enrichment_p_positive) <- paste0("C", 1:8)
for (i in unique(TF_peak_gene_cor_select$TF)) {
  temp <- TF_peak_gene_cor_select_positive[TF_peak_gene_cor_select_positive$TF == i,]
  if (nrow(temp) == 0) {
    TF_cluster_gene_enrichment_p_positive[i, paste0("C",j)] <- 1
  } else {
    total_targets <- unique(temp$Target)
    for (j in c("1","2","3","4","5","6","7","8")) {
      total_cluster_genes <- rownames(RNA_monocle3_dynamic_genes_cluster)[as.character(RNA_monocle3_dynamic_genes_cluster$Cluster) == j]
      total_cluster_genes <- total_cluster_genes[total_cluster_genes %in% unique(TF_peak_gene_cor_select$Target)]
      common_genes <- intersect(total_targets, total_cluster_genes)
      p_value <- 1 - phyper(length(common_genes) - 1,
                            length(total_cluster_genes),
                            902 - length(total_cluster_genes),
                            length(total_targets), lower.tail = TRUE)
      if (p_value == 0) {
        p_value <- 1e-30
      }
      TF_cluster_gene_enrichment_p_positive[i, paste0("C",j)] <- p_value
    }
  }
}
TF_cluster_gene_enrichment_p_negative <- matrix(1, nrow = length(unique(TF_peak_gene_cor_select$TF)), ncol = 8)
rownames(TF_cluster_gene_enrichment_p_negative) <- unique(TF_peak_gene_cor_select$TF)
colnames(TF_cluster_gene_enrichment_p_negative) <- paste0("C", 1:8)
for (i in unique(TF_peak_gene_cor_select$TF)) {
  temp <- TF_peak_gene_cor_select_negative[TF_peak_gene_cor_select_negative$TF == i,]
  if (nrow(temp) == 0) {
    TF_cluster_gene_enrichment_p_negative[i, paste0("C",j)] <- 1
  } else {
    total_targets <- unique(temp$Target)
    for (j in c("1","2","3","4","5","6","7","8")) {
      total_cluster_genes <- rownames(RNA_monocle3_dynamic_genes_cluster)[as.character(RNA_monocle3_dynamic_genes_cluster$Cluster) == j]
      total_cluster_genes <- total_cluster_genes[total_cluster_genes %in% unique(TF_peak_gene_cor_select$Target)]
      common_genes <- intersect(total_targets, total_cluster_genes)
      p_value <- 1 - phyper(length(common_genes) - 1,
                            length(total_cluster_genes),
                            902 - length(total_cluster_genes),
                            length(total_targets), lower.tail = TRUE)
      if (p_value == 0) {
        p_value <- 1e-30
      }
      TF_cluster_gene_enrichment_p_negative[i, paste0("C",j)] <- p_value
    }
  }
}



{
  TF_selected <- unique(TF_peak_gene_cor_select$TF)
  smoothMat_MM_matrix_selected <- smoothMat_MM_matrix[TF_selected,]
  smoothMat_MM_matrix_selected <- t(smoothMat_MM_matrix_selected)
  TF_reorder <- c()
  for (i in 1:ncol(smoothMat_MM_matrix_selected)) {
    TF_reorder <- c(TF_reorder,
                    which.max(smoothMat_MM_matrix_selected[,i]))
  }
  names(TF_reorder) <- colnames(smoothMat_MM_matrix_selected)
  TF_reorder <- sort(TF_reorder, decreasing = FALSE)
  TF_reorder <- names(TF_reorder)
  
  for (i in 1:ncol(smoothMat_MM_matrix_selected)) {
    smoothMat_MM_matrix_selected[,i] <- (smoothMat_MM_matrix_selected[,i] - mean(smoothMat_MM_matrix_selected[,i])) / sd(smoothMat_MM_matrix_selected[,i])
  }
  smoothMat_MM_matrix_selected <- t(smoothMat_MM_matrix_selected)
  smoothMat_MM_matrix_selected[smoothMat_MM_matrix_selected > 2] <- 2
  smoothMat_MM_matrix_selected[smoothMat_MM_matrix_selected < -2] <- -2
  smoothMat_MM_matrix_selected <- smoothMat_MM_matrix_selected[TF_reorder,]
  h1 <- Heatmap(smoothMat_MM_matrix_selected,
                show_column_names = FALSE,
                cluster_rows = FALSE, cluster_columns = FALSE,
                col = paletteContinuous(set = "solarExtra", n = 100),
                row_names_side = "left",
                name = "Row scaled TF deviation",
                width = unit(2.5, "null"),
                heatmap_legend_param = list(
                  # title_position = "topcenter",  # 图例标题居中
                  legend_direction = "horizontal"  # 图例方向为水平
                )
  )
  h1
  
  smoothMat_Exp_matrix_selected <- smoothMat_Exp_matrix[rownames(smoothMat_MM_matrix_selected),]
  smoothMat_Exp_matrix_selected <- t(smoothMat_Exp_matrix_selected)
  for (i in 1:ncol(smoothMat_Exp_matrix_selected)) {
    smoothMat_Exp_matrix_selected[,i] <- (smoothMat_Exp_matrix_selected[,i] - mean(smoothMat_Exp_matrix_selected[,i])) / sd(smoothMat_Exp_matrix_selected[,i])
  }
  smoothMat_Exp_matrix_selected <- t(smoothMat_Exp_matrix_selected)
  smoothMat_Exp_matrix_selected[smoothMat_Exp_matrix_selected > 2] <- 2
  smoothMat_Exp_matrix_selected[smoothMat_Exp_matrix_selected < -2] <- -2
  # smoothMat_GS_matrix_selected <- smoothMat_GS_matrix_selected[TF_reorder,]
  h2 <- Heatmap(smoothMat_Exp_matrix_selected,
                show_column_names = FALSE,
                cluster_rows = FALSE, cluster_columns = FALSE,
                col = paletteContinuous(set = "beach", n = 100),
                row_names_side = "left",
                name = "Row scaled TF expression",
                width = unit(2.5, "null"),
                heatmap_legend_param = list(
                  # title_position = "topcenter",  # 图例标题居中
                  legend_direction = "horizontal"  # 图例方向为水平
                )
  )
  h2
  
  TF_cluster_gene_enrichment_p_positive_selected <- TF_cluster_gene_enrichment_p_positive[rownames(smoothMat_MM_matrix_selected),]
  TF_cluster_gene_enrichment_p_positive_selected <- -log2(TF_cluster_gene_enrichment_p_positive_selected)
  # for (i in 1:nrow(TF_cluster_gene_enrichment_p_positive_selected)) {
  #   TF_cluster_gene_enrichment_p_positive_selected[i,] <- (TF_cluster_gene_enrichment_p_positive_selected[i,] - mean(TF_cluster_gene_enrichment_p_positive_selected[i,])) / sd(TF_cluster_gene_enrichment_p_positive_selected[i,])
  # }
  temp <- TF_cluster_gene_enrichment_p_positive[rownames(smoothMat_MM_matrix_selected),]
  TF_cluster_gene_enrichment_p_positive_selected[temp >= 0.05] <- NA
  TF_cluster_gene_enrichment_p_positive_selected[TF_cluster_gene_enrichment_p_positive_selected > 10] <- 10
  h3 <- Heatmap(TF_cluster_gene_enrichment_p_positive_selected[,c("C3", "C6", "C7",
                                            "C1", "C8", "C5",
                                            "C4", "C2")],
                cluster_rows = FALSE, cluster_columns = FALSE,
                col = paletteContinuous(set = "coolwarm", n = 100),
                row_names_side = "left",
                width = unit(1.5, "null"),
                name = "Positive regulation\nEnrichment significance\n-log2(P-value)",
                heatmap_legend_param = list(
                  # title_position = "topcenter",  # 图例标题居中
                  legend_direction = "horizontal"  # 图例方向为水平
                )
  )
  h3
  
  TF_cluster_gene_enrichment_p_negative_selected <- TF_cluster_gene_enrichment_p_negative[rownames(smoothMat_MM_matrix_selected),]
  TF_cluster_gene_enrichment_p_negative_selected <- -log2(TF_cluster_gene_enrichment_p_negative_selected)
  # for (i in 1:nrow(TF_cluster_gene_enrichment_p_positive_selected)) {
  #   TF_cluster_gene_enrichment_p_positive_selected[i,] <- (TF_cluster_gene_enrichment_p_positive_selected[i,] - mean(TF_cluster_gene_enrichment_p_positive_selected[i,])) / sd(TF_cluster_gene_enrichment_p_positive_selected[i,])
  # }
  temp <- TF_cluster_gene_enrichment_p_negative[rownames(smoothMat_MM_matrix_selected),]
  TF_cluster_gene_enrichment_p_negative_selected[temp >= 0.05] <- NA
  TF_cluster_gene_enrichment_p_negative_selected[TF_cluster_gene_enrichment_p_negative_selected > 10] <- 10
  h4 <- Heatmap(TF_cluster_gene_enrichment_p_negative_selected[,c("C3", "C6", "C7",
                                                                  "C1", "C8", "C5",
                                                                  "C4", "C2")],
                cluster_rows = FALSE, cluster_columns = FALSE,
                col = paletteContinuous(set = "greenBlue", n = 100),
                row_names_side = "left",
                width = unit(1.5, "null"),
                name = "Negative regulation\nEnrichment significance\n-log2(P-value)",
                heatmap_legend_param = list(
                  # title_position = "topcenter",  # 图例标题居中
                  legend_direction = "horizontal"  # 图例方向为水平
                )
  )
  h4
  pdf("ATAC_Monocle3_RNA_ATAC_TF_selected_Heatmap.pdf",
      width = 12, height = 15)
  h1 + h2 + h3 + h4
  dev.off()
}

temp <- RNA_monocle3_dynamic_genes_cluster[TF_selected,,drop = FALSE]
temp <- as.data.frame(na.omit(temp))
TF_selected %in% rownames(RNA_monocle3_dynamic_genes_cluster)
