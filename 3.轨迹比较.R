percent_bar <- function(meta, Ident, Group,
                        fill_color = NULL, fill_label = "",
                        flow = FALSE) {
  if (flow == TRUE) {
    temp <- data.frame(Ident = meta[,Ident],
                       Group = meta[,Group])
    temp <- as.data.frame.array(table(temp$Ident, temp$Group))
    for (i in 1:ncol(temp)) {
      temp[,i] <- temp[,i] / sum(temp[,i])
    }
    temp <- data.frame(Ident = rownames(temp),
                       temp)
    dat <- reshape::melt(temp, id = "Ident")
    dat$Ident <- factor(dat$Ident, levels = sort(unique(meta[,Ident])))
    
    p <- ggplot(dat, aes(x = variable, y = value, fill = Ident, 
                         stratum = Ident, alluvium = Ident)) +
      ggalluvial::geom_stratum(width = 0.55) +  #代替 geom_col() 绘制堆叠柱形图
      ggalluvial::geom_flow(width = 0.55, alpha = 1) +  #绘制同类别之间的连接线\
      scale_y_continuous(labels = scales::percent) +
      theme_bw() +
      labs(x = "", y = "Percentage", fill = fill_label) +
      theme(axis.text = element_text(size = 13, color = "black"),
            axis.title = element_text(size = 13, color = "black"),
            legend.text = element_text(size = 13, color = "black"),
            legend.title = element_text(size = 13, color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
  } else {
    temp <- data.frame(Ident = meta[,Ident],
                       Group = meta[,Group])
    temp$Ident <- factor(temp$Ident,
                         levels = sort(unique((meta[,Ident]))))
    temp$Group <- factor(temp$Group,
                         levels = sort(unique((meta[,Group]))))
    p <- ggplot(data = temp, aes(x = Group, fill = Ident)) +
      geom_bar(stat = "count", position = "fill", width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      theme_bw() +
      labs(x = "", y = "Percentage", fill = fill_label) +
      theme(axis.text = element_text(size = 13, color = "black"),
            axis.title = element_text(size = 13, color = "black"),
            legend.text = element_text(size = 13, color = "black"),
            legend.title = element_text(size = 13, color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }
  
  if (!is.null(fill_color)) {
    p <- p + scale_fill_manual(values = fill_color)
  }
  p
}
jitter_linebox <- function(meta, group.by = "Day",
                           features = c("nFeature_RNA"),
                           group.color = group_color,
                           linewidth = 2, jitter_width,
                           stat_cor_size, font_size = 10) {
  temp <- lapply(unique(meta[,group.by]), function(x){
    temp <- meta[which(meta[,group.by] == as.character(x)),]
    quan <- quantile(temp[, features])
    IQR <- quan[4] - quan[2]
    Q1 <- quan[2]
    Q2 <- quan[3]
    Q3 <- quan[4]
    low <- Q1-1.5*IQR
    if (low <= min(temp[, features])) {
      low <- min(temp[, features])
    }
    high <- Q3+1.5*IQR
    if (high >= max(temp[, features])) {
      high <- max(temp[, features])
    }
    data.frame(group = as.character(x),
               low = low,
               Q1 = Q1,
               Q2 = Q2,
               Q3 = Q3,
               high = high
    )
  })
  temp <- rlist::list.rbind(temp)
  temp <- as.data.frame(temp)
  temp$group <- factor(temp$group,
                       levels = levels(meta[, group.by]))
  g <- ggplot() +
    geom_jitter(data = meta, aes(x = meta[, group.by], y = meta[, features],
                                 color = meta[, group.by]),
                width = jitter_width) +
    scale_color_manual(values = group.color) +
    geom_segment(aes(x = group, y = low, xend = group, yend = high), linewidth = 0.5, colour = "black", data = temp) +
    geom_segment(aes(x = group, y = Q1, xend = group, yend = Q3), linewidth = linewidth, colour = "black", data = temp) +
    geom_point(aes(x = group, y = Q2), data = temp, colour = "white", size = linewidth) +
    geom_smooth(aes(x = as.numeric(group), y = Q2), data = temp,
                method = "lm", colour = "red", se = FALSE) +
    ggpubr::stat_cor(aes(x = as.numeric(group), y = Q2), data = temp,
                     method = "pearson", size = stat_cor_size) +
    theme_bw() +
    ggtitle(features) +
    theme(axis.text = element_text(size = font_size, colour = "black"),
          axis.title = element_text(size = font_size, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(size = font_size, face = "bold", hjust = 0.5),
          panel.grid = element_blank()) +
    labs(x = "", y = features, fill = group.by)
  return(g)
}
violin_linebox <- function(meta, group.by = "Day", features = c("nFeature_RNA"), group.color = group_color, linewidth = 2) {
  temp <- lapply(unique(meta[,group.by]), function(x){
    temp <- meta[which(meta[,group.by] == as.character(x)),]
    quan <- quantile(temp[, features])
    IQR <- quan[4] - quan[2]
    Q1 <- quan[2]
    Q2 <- quan[3]
    Q3 <- quan[4]
    low <- Q1-1.5*IQR
    if (low <= min(temp[, features])) {
      low <- min(temp[, features])
    }
    high <- Q3+1.5*IQR
    if (high >= max(temp[, features])) {
      high <- max(temp[, features])
    }
    data.frame(group = as.character(x),
               low = low,
               Q1 = Q1,
               Q2 = Q2,
               Q3 = Q3,
               high = high
    )
  })
  temp <- rlist::list.rbind(temp)
  g <- ggplot() +
    geom_violin(data = meta, aes(x = meta[,group.by], y = meta[, features], fill = meta[, group.by]),
                scale = "width") +
    scale_fill_manual(values = group.color) +
    geom_segment(aes(x = group, y = low, xend = group, yend = high), linewidth = 0.5, colour = "black", data = temp) +
    geom_segment(aes(x = group, y = Q1, xend = group, yend = Q3), linewidth = linewidth, colour = "black", data = temp) +
    geom_point(aes(x = group, y = Q2), data = temp, colour = "white", size = linewidth) +
    theme_bw() +
    ggtitle(features) +
    theme(axis.text = element_text(size = 14, colour = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          panel.grid = element_blank()) +
    labs(x = "", y = features, fill = group.by)
  return(g)
}

setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
library(Seurat)
library(monocle)
library(monocle3)
library(harmony)
library(tidyverse)
library(patchwork)
library(BSgenome.malu)
library(ArchR)
library(ggalt)
library(parallel)
library(clustree)
library(dplyr)
library(patchwork)
library(BSgenome.malu)
library(slingshot)
library(SingleCellExperiment)
library(qs)
library(tidyverse)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(pheatmap)

rscripts <- list.files(path = "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/ArchR_R", pattern = ".R$", recursive = F, full.names = F)
for (i in rscripts) {
  source(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/ArchR_R/",i))
}
rscripts_fixed <- list.files(path = "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/ArchR_R_bug_fixed", pattern = ".R$", recursive = F, full.names = F)
for (i in rscripts_fixed) {
  source(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/ArchR_R_bug_fixed/",i))
}

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.malu)

gene_file <- "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/wx_create_gene_malu.txt"
tss_file <- "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/wx_create_tss_malu.txt"
exon_file <- "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/wx_create_exon_malu.txt"

tss <- read.csv(tss_file, sep = "\t", header = F)
colnames(tss) <- c("chr", "start", "end", "strand", "tx_id", "ts_name")
tss$end <- tss$start
tss1 <- tss[-grep("chrMT", tss$chr, fixed = T),]
tss1 <- makeGRangesFromDataFrame(tss1, keep.extra.columns = T)

gene <- read.csv(gene_file, sep = "\t", header = F)
colnames(gene) <- c("chr", "start", "end", "strand", "gene_id", "symbol")
gene$symbol <- make.unique(gene$symbol)
gene1 <- gene[-grep("chrMT", gene$chr, fixed = T),]
gene1 <- makeGRangesFromDataFrame(gene1, keep.extra.columns = T)

exon <- read.csv(exon_file, sep = "\t", header = F)
colnames(exon) <- c("chr", "start", "end", "strand", "gene_id", "symbol")
exon$symbol <- make.unique(exon$symbol)
exon1 <- exon[-grep("chrMT", exon$chr, fixed = T),]
exon1 <- makeGRangesFromDataFrame(exon1, keep.extra.columns = T)

geneAnnotation <- createGeneAnnotation(TSS = tss1,
                                       exons = exon1,
                                       genes = gene1)
addArchRThreads(threads = 30)

RNA_Celltype_color <- c("#A4C9DD", "#2572A9", "#ADD487", "#399938",
                        "#F19695", "#D5231E", "#F5BB6F", "#EF7C1C",
                        "#C6B0D2", "#653B90")
names(RNA_Celltype_color) <- c("AnMCs", "Proliferative_progenitor cells", "Progenitor cells",
                               "Chondrocytes", "Hypertrophic chondrocytes", "Chondroclasts",
                               "Mural cells", "Endothelial cells", "Monocytes_Macrophages",
                               "Mast cells")
Discs_sub_color <- c("#A4C9DD", "#54989E", "#51A743",
                     "#ED8684", "#E46A45", "#F0841F",
                     "#AD91C0", "#C4B497", "#AA5728")
names(Discs_sub_color) <- c("NPPCs", "Stroma", "HomCs",
                            "preHTCs", "HTCs", "Notochord",
                            "Pericyte", "EC", "blood")
Limb_sub_color <- c("#a6cee3", "#2f83af", "#95d175",
                    "#759e50", "#f06161")
names(Limb_sub_color) <- c("OCPs", "Osteoprogenitors",
                           "PMSC", "Chondroblasts",
                           "Chondrocytes")
Tumor_sub_color <- c("#379938", "#F19695")
names(Tumor_sub_color) <- c("Tumor_Osteoblastic", "Tumor_Chondroblasts")

RNA_monocle3 <- readRDS("RNA_monocle3_cds.rds")
Limb_monocle3 <- readRDS("3.轨迹比较.Limb_monocle3_cds.rds")
Tumor_monocle3 <- readRDS("3.轨迹比较.Tumor_monocle3_cds.rds")
Discs_monocle3 <- readRDS("3.轨迹比较.Discs_monocle3_cds.rds")
Limb <- readRDS("./OtherDataSets/sce.0.13.integrate_umap_last.3000.RDS")
Limb_sub <- subset(Limb, celltype %in% c("OCPs", "Osteoprogenitors",
                                         "PMSC", "Chondroblasts",
                                         "Chondrocytes"))
Tumor <- readRDS("./OtherDataSets/tumor_sce.0.15.integrate_umap_last.3000.ham.RDS")
Tumor_sub <- subset(Tumor, celltype %in% c("Tumor_Osteoblastic", "Tumor_Chondroblasts"))

Discs <- readRDS("./OtherDataSets/disc_sce.0.3.integrate_umap_last.3000.RDS")
Discs_sub <- subset(Discs, celltype %in% c("HomCs", "preHTCs", "HTCs"))

RNA_singlet_sub <- readRDS("./RNA_singlet_sub.rds")
RNA_singlet_sub <- NormalizeData(RNA_singlet_sub)
RNA_singlet_sub <- FindVariableFeatures(RNA_singlet_sub)
RNA_singlet_sub_pseudotime <- pseudotime(RNA_monocle3)
RNA_singlet_sub_pseudotime <- (RNA_singlet_sub_pseudotime - min(RNA_singlet_sub_pseudotime)) / (max(RNA_singlet_sub_pseudotime) - min(RNA_singlet_sub_pseudotime))
RNA_singlet_sub$pseudotime <- RNA_singlet_sub_pseudotime[colnames(RNA_singlet_sub)]

Discs_sub_pseudotime <- pseudotime(Discs_monocle3)
Discs_sub_pseudotime <- (Discs_sub_pseudotime - min(Discs_sub_pseudotime)) / (max(Discs_sub_pseudotime) - min(Discs_sub_pseudotime))
Discs_sub$pseudotime <- Discs_sub_pseudotime[colnames(Discs_sub)]
colnames(Discs_sub)

Limb_sub_pseudotime <- pseudotime(Limb_monocle3)
Limb_sub_pseudotime <- (Limb_sub_pseudotime - min(Limb_sub_pseudotime)) / (max(Limb_sub_pseudotime) - min(Limb_sub_pseudotime))
Limb_sub$pseudotime <- Limb_sub_pseudotime[colnames(Limb_sub)]

Tumor_sub_pseudotime <- pseudotime(Tumor_monocle3)
Tumor_sub_pseudotime <- (Tumor_sub_pseudotime - min(Tumor_sub_pseudotime)) / (max(Tumor_sub_pseudotime) - min(Tumor_sub_pseudotime))
Tumor_sub$pseudotime <- Tumor_sub_pseudotime[colnames(Tumor_sub)]

# bin_width <- max(RNA_singlet_sub$pseudotime) / 100
# Pseudotime_bins <- cut(RNA_singlet_sub$pseudotime,
#                        breaks = seq(0, max(RNA_singlet_sub$pseudotime), by = bin_width),
#                        include.lowest = TRUE)
# RNA_singlet_sub$Pseudotime_bin <- Pseudotime_bins
# sort(table(RNA_singlet_sub$Pseudotime_bin))

#####
# CellAlign
library(cellAlign)
CellAlign_Global_fun <- function(interGlobal_ref, interGlobal_query,
                          ref_Traj, query_Traj, numPts = 200,
                          sharedMarkers
                          ){
  interScaledGlobal_ref <- cellAlign::scaleInterpolate(interGlobal_ref)
  interScaledGlobal_query <- cellAlign::scaleInterpolate(interGlobal_query)
  selected_ref <- interScaledGlobal_ref$scaledData[sharedMarkers,]
  selected_query <- interScaledGlobal_query$scaledData[sharedMarkers,]
  ref_query_alignment <- globalAlign(y = selected_ref,
                                     x = selected_query,
                                     # scores = list(query = selected_query$traj, 
                                     #               ref = selected_ref$traj),
                                     sigCalc = T, numPerm = 50)
  ref_query_mapping <- mapRealDataGlobal(ref_query_alignment,
                                         intTrajQuery = interScaledGlobal_query$traj,
                                         realTrajQuery = query_Traj,
                                         intTrajRef = interScaledGlobal_ref$traj,
                                         realTrajRef = ref_Traj)
  return(list("ref_query_alignment" = ref_query_alignment,
              "ref_query_mapping" = ref_query_mapping))
}
CellAlign_Local_fun <- function(interGlobal_ref,
                                interGlobal_query,
                                ref_Traj, query_Traj,
                                numPts = 200, Thresh = 0.2,
                                sharedMarkers) {
  interScaledGlobal_ref <- cellAlign::scaleInterpolate(interGlobal_ref)
  interScaledGlobal_query <- cellAlign::scaleInterpolate(interGlobal_query)
  selected_ref <- interScaledGlobal_ref$scaledData[sharedMarkers,]
  selected_query <- interScaledGlobal_query$scaledData[sharedMarkers,]
  A <- calcDistMat(x = selected_query,
                   y = selected_query,
                   dist.method = 'Euclidean')
  A[A > 10*Thresh] <- max(A)
  ref_query_alignment <- localAlign(x = selected_query,
                          y = selected_ref,
                          threshPercent = Thresh)
  
  ref_query_mapping <- mapRealDataLocal(alignment = ref_query_alignment,
                   intTrajQuery = interScaledGlobal_query$traj,
                   realTrajQuery = query_Traj,
                   intTrajRef = interScaledGlobal_ref$traj,
                   realTrajRef = ref_Traj)
  return(list("A" = A,
              "ref_query_alignment" = ref_query_alignment,
              "ref_query_mapping" = ref_query_mapping))
}
metaNodesMapping <- function(alignment, intTrajQuery, realTrajQuery,
                             intTrajRef, realTrajRef){
  alignMat = data.frame(queryInt = alignment$align[[1]]$index1, refInt = alignment$align[[1]]$index2)
  alignMat$ptQueryInt = intTrajQuery[alignMat$queryInt]
  alignMat$ptRefInt = intTrajRef[alignMat$refInt]
  alignMat$ptQuery = sapply(alignMat$ptQueryInt, function(ptQInt){
    return(realTrajQuery[which.min(abs(realTrajQuery - ptQInt))])})
  alignMat$ptRef = sapply(alignMat$ptRefInt, function(ptRInt){
    return(realTrajRef[which.min(abs(realTrajRef - ptRInt))])})
  alignMat$cellIDQuery = sapply(alignMat$ptQueryInt, function(ptQInt){
    return(names(realTrajQuery)[which.min(abs(realTrajQuery - ptQInt))])})
  alignMat$cellIDRef = sapply(alignMat$ptRefInt, function(ptRInt){
    return(names(realTrajRef)[which.min(abs(realTrajRef - ptRInt))])})
  alignMat = unique(alignMat)
  
  return(alignMat)
}
library(cellAlign)
###### RNA & Discs
{
  numPts <- 200
  # interGlobal_Discs_sub <- cellAlign::interWeights(expDataBatch = Discs_sub@assays$RNA@data,
  #                                                  trajCond = Discs_sub_pseudotime,
  #                                                  winSz = 0.1,
  #                                                  numPts = numPts)
  # saveRDS(interGlobal_Discs_sub, "./轨迹比较/interGlobal_Discs_sub_3.rds")
  interGlobal_RNA_sub <- readRDS("./轨迹比较/interGlobal_RNA_sub.rds")
  interGlobal_RNA_sub$interpolatedVals <- interGlobal_RNA_sub$interpolatedVals[rowSums(interGlobal_RNA_sub$interpolatedVals) > 0,]
  interGlobal_RNA_sub$error <- interGlobal_RNA_sub$error[rownames(interGlobal_RNA_sub$interpolatedVals),]
  interGlobal_Discs_sub <- readRDS("./轨迹比较/interGlobal_Discs_sub_3.rds")
  interGlobal_Discs_sub$interpolatedVals <- interGlobal_Discs_sub$interpolatedVals[rowSums(interGlobal_Discs_sub$interpolatedVals) > 0,]
  interGlobal_Discs_sub$error <- interGlobal_Discs_sub$error[rownames(interGlobal_Discs_sub$interpolatedVals),]
  
  sharedMarkers <- intersect(rownames(interGlobal_Discs_sub$interpolatedVals),
                             rownames(interGlobal_RNA_sub$interpolatedVals))
  length(sharedMarkers)
  
  RNA_Discs_CellAlign_Global <- CellAlign_Global_fun(interGlobal_ref = interGlobal_RNA_sub,
                                                     interGlobal_query = interGlobal_Discs_sub,
                                                     ref_Traj = RNA_singlet_sub$pseudotime,
                                                     query_Traj = Discs_sub$pseudotime,
                                                     numPts = 200, sharedMarkers = sharedMarkers)
  
  plotAlign(RNA_Discs_CellAlign_Global$ref_query_alignment)
  
  plotMapping(RNA_Discs_CellAlign_Global$ref_query_mapping)
  
  Global_node_mapping <- RNA_Discs_CellAlign_Global[["ref_query_mapping"]][["metaNodesPt"]]
  Global_node_mapping <- Global_node_mapping[,-5]
  Global_node_mapping$index <-paste0(Global_node_mapping$ptRef,
                                     "_", Global_node_mapping$ptQuery)
  
  RNA_Discs_CellAlign_Local <- CellAlign_Local_fun(interGlobal_ref = interGlobal_RNA_sub,
                                                   interGlobal_query = interGlobal_Discs_sub,
                                                   ref_Traj = RNA_singlet_sub$pseudotime,
                                                   query_Traj = Discs_sub$pseudotime,
                                                   numPts = 200, Thresh = 0.2,
                                                   sharedMarkers = sharedMarkers)
  
  plotAlign(RNA_Discs_CellAlign_Local$ref_query_alignment)
  
  plotMapping(RNA_Discs_CellAlign_Local$ref_query_mapping)
  
  Local_node_mapping <- RNA_Discs_CellAlign_Local[["ref_query_mapping"]][["metaNodesPt"]]
  Local_node_mapping$index <- paste0(Local_node_mapping$ptRef,
                                     "_", Local_node_mapping$ptQuery)
  sum(Local_node_mapping$index %in% Global_node_mapping$index)
  Global_node_mapping$Mapping <- ifelse(Global_node_mapping$index %in% Local_node_mapping$index,
                                        "Conserved", "Not")
  Global_node_mapping$rowIndex <- 1:nrow(Global_node_mapping)
  
  Global_ref_node_cells <- Global_node_mapping[,c("metaNodeRef",
                                                  "ptRef")]
  Global_ref_node_cells <- Global_ref_node_cells[!duplicated(Global_ref_node_cells),]
  Global_ref_node_cells$ptRef <- paste0("ptRef", "_", Global_ref_node_cells$ptRef)
  rownames(Global_ref_node_cells) <- Global_ref_node_cells$ptRef
  
  Global_query_node_cells <- Global_node_mapping[,c("metaNodeQuery",
                                                    "ptQuery")]
  Global_query_node_cells <- Global_query_node_cells[!duplicated(Global_query_node_cells),]
  Global_query_node_cells$ptQuery <- paste0("ptQuery", "_", Global_query_node_cells$ptQuery)
  rownames(Global_query_node_cells) <- Global_query_node_cells$ptQuery
  node_cells <- data.frame(node = c(Global_ref_node_cells$ptRef,
                                    Global_query_node_cells$ptQuery),
                           cells = c(Global_ref_node_cells$metaNodeRef,
                                     Global_query_node_cells$metaNodeQuery))
  rownames(node_cells) <- node_cells$node
  
  Global_node_mapping_df <- melt(Global_node_mapping[, c("ptQuery", "ptRef", "rowIndex")],
                                 id.vars = c("rowIndex"))
  Global_node_mapping_df$cells <- paste0(Global_node_mapping_df$variable,
                                         "_", Global_node_mapping_df$value)
  Global_node_mapping_df$cells <- node_cells[Global_node_mapping_df$cells, "cells"]
  cells <- c(RNA_Discs_CellAlign_Global$ref_query_mapping$refAssign,
             RNA_Discs_CellAlign_Global$ref_query_mapping$queryAssign)
  main_celltype <- c()
  for (i in 1:nrow(Global_node_mapping_df)) {
    if (Global_node_mapping_df$variable[i] == "ptQuery") {
      temp <- table(Discs_sub@meta.data[cells[[Global_node_mapping_df$cells[i]]], "celltype"])
      temp_celltype <- names(temp)[which.max(temp)]
      main_celltype <- c(main_celltype,
                         temp_celltype)
    }
    if (Global_node_mapping_df$variable[i] == "ptRef") {
      temp <- table(RNA_singlet_sub@meta.data[cells[[Global_node_mapping_df$cells[i]]], "Celltype_rename"])
      temp_celltype <- names(temp)[which.max(temp)]
      main_celltype <- c(main_celltype,
                         temp_celltype)
    }
  }
  Global_node_mapping_df$main_celltype <- main_celltype
  Global_node_mapping_df_Conserved <- melt(Global_node_mapping[Global_node_mapping$Mapping == "Conserved",
                                                               c("ptQuery", "ptRef", "rowIndex")],
                                           id.vars = c("rowIndex"))
  Global_node_mapping_df_Not <- melt(Global_node_mapping[Global_node_mapping$Mapping == "Not",
                                                         c("ptQuery", "ptRef", "rowIndex")],
                                     id.vars = c("rowIndex"))
  pdf("./轨迹比较/CellAlign_1.point_alignment_Antler_Discs.pdf",
      height = 2.5, width = 7.5)
  ggplot() + 
    theme_bw() +
    geom_line(data = Global_node_mapping_df_Conserved,
              aes(x = variable, y = value, group = rowIndex),
              color = "red") +
    geom_line(data = Global_node_mapping_df_Not,
              aes(x = variable, y = value, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df,
               aes(x = variable, y = value,
                   colour = main_celltype),
               shape = 15) +
    scale_color_manual(values = c(RNA_Celltype_color,
                                  Discs_sub_color)) +
    coord_flip() + ggtitle("Ref_Antler vs. Query_Discs") +
    labs(color = "Celltype") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          legend.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  dev.off()
  
  costMat <- RNA_Discs_CellAlign_Global$ref_query_alignment$localCostMatrix
  costMat <- t(apply(costMat, 1, function(x) {
    return(as.numeric(x))
  })) # 列 ref, 行 query
  linearInd <- sub2ind(nrow(costMat), RNA_Discs_CellAlign_Global$ref_query_alignment$align[[1]]$index1, 
                       RNA_Discs_CellAlign_Global$ref_query_alignment$align[[1]]$index2)
  costMat[linearInd] <- NA
  costMat <- as.data.frame(costMat)
  
  metaNodesMap <- metaNodesMapping(RNA_Discs_CellAlign_Global$ref_query_alignment,
                                   intTrajQuery = interGlobal_Discs_sub$traj, 
                                   realTrajQuery = Discs_sub$pseudotime,
                                   intTrajRef = interGlobal_RNA_sub$traj,
                                   realTrajRef = RNA_singlet_sub$pseudotime)
  metaNodesMap_Query <- metaNodesMap[!duplicated(metaNodesMap$queryInt),]
  metaNodesMap_Query <- metaNodesMap_Query[order(metaNodesMap_Query$queryInt,
                                                 decreasing = F),]
  metaNodesMap_Ref <- metaNodesMap[!duplicated(metaNodesMap$refInt),]
  metaNodesMap_Ref <- metaNodesMap_Ref[order(metaNodesMap_Ref$refInt,
                                             decreasing = F),]
  metaNodesMap_Ref$ptRef <- paste0("ptRef_",
                                   metaNodesMap_Ref$ptRef)
  metaNodesMap_Ref <- metaNodesMap_Ref[!duplicated(metaNodesMap_Ref$ptRef),]
  metaNodesMap_Query$ptQuery <- paste0("ptQuery_",
                                       metaNodesMap_Query$ptQuery)
  metaNodesMap_Query <- metaNodesMap_Query[!duplicated(metaNodesMap_Query$ptQuery),]
  
  costMat <- costMat[,metaNodesMap_Ref$refInt]
  costMat <- costMat[metaNodesMap_Query$queryInt,]
  colnames(costMat) <- metaNodesMap_Ref$ptRef # 列 ref
  rownames(costMat) <- metaNodesMap_Query$ptQuery # 行 query
  
  col_celltype <- c()
  for (i in 1:ncol(costMat)) {
    temp <- table(RNA_singlet_sub@meta.data[cells[[metaNodesMap_Ref$cellIDRef[i]]], "Celltype_rename"])
    temp_celltype <- names(temp)[which.max(temp)]
    col_celltype <- c(col_celltype,
                      temp_celltype)
  }
  row_celltype <- c()
  for (i in 1:nrow(costMat)) {
    temp <- table(Discs_sub@meta.data[cells[[metaNodesMap_Query$cellIDQuery[i]]], "celltype"])
    temp_celltype <- names(temp)[which.max(temp)]
    row_celltype <- c(row_celltype,
                      temp_celltype)
  }
  
  Global_node_mapping_df_Conserved$node <- paste0(Global_node_mapping_df_Conserved$variable, "_",
                                                  Global_node_mapping_df_Conserved$value)
  annotCols <- data.frame(Ref_celltype = col_celltype,
                          Local_alignment = "Not conserved",
                          row.names = colnames(costMat))
  annotCols[Global_node_mapping_df_Conserved$node[Global_node_mapping_df_Conserved$variable == "ptRef"], "Local_alignment"] <- "Conserved"
  annotRows <- data.frame(Query_celltype = row_celltype,
                          Local_alignment = "Not conserved",
                          row.names = rownames(costMat))
  annotRows[Global_node_mapping_df_Conserved$node[Global_node_mapping_df_Conserved$variable == "ptQuery"], "Local_alignment"] <- "Conserved"
  
  annotColors <- list(Ref_celltype = RNA_Celltype_color,
                      Query_celltype = Discs_sub_color,
                      Local_alignment = c("Not conserved" = "gray", "Conserved" = "red"))
  pdf("./轨迹比较/CellAlign_2.alignment_pheatmap_Antler_Discs.pdf",
      width = 8, height = 8)
  pheatmap(costMat, na_col = "purple",
           show_rownames = F, show_colnames = F,
           annotation_row = annotRows,
           annotation_col = annotCols,
           annotation_colors = annotColors,
           cluster_rows = F, cluster_cols = F,
           cellheight = 1.5, cellwidth = 1.5)
  dev.off()
  
  interScaledLocal_RNA_sub <- cellAlign::scaleInterpolate(interGlobal_RNA_sub)
  interScaledLocal_Discs_sub <- cellAlign::scaleInterpolate(interGlobal_Discs_sub)
  ref_colMeans <- colMeans(interScaledLocal_RNA_sub$scaledData[sharedMarkers,])
  ref_colMeans <- ref_colMeans[metaNodesMap_Ref$refInt]
  names(ref_colMeans) <- colnames(costMat)
  query_colMeans <- colMeans(interScaledLocal_Discs_sub$scaledData[sharedMarkers,])
  query_colMeans <- query_colMeans[metaNodesMap_Query$queryInt]
  names(query_colMeans) <- rownames(costMat)
  
  metaNodesMap_Ref_exp <- interScaledLocal_RNA_sub$scaledData[sharedMarkers, metaNodesMap_Ref$refInt]
  metaNodesMap_Ref_exp <- data.frame(t(metaNodesMap_Ref_exp),
                                     check.rows = F, check.names = F)
  rownames(metaNodesMap_Ref_exp) <- colnames(costMat)
  metaNodesMap_Query_exp <- interScaledLocal_Discs_sub$scaledData[sharedMarkers, metaNodesMap_Query$queryInt]
  metaNodesMap_Query_exp <- data.frame(t(metaNodesMap_Query_exp),
                                       check.rows = F, check.names = F)
  rownames(metaNodesMap_Query_exp) <- rownames(costMat)
  
  Global_node_mapping_df_ref <- Global_node_mapping_df[Global_node_mapping_df$variable == "ptRef",]
  Global_node_mapping_df_query <- Global_node_mapping_df[Global_node_mapping_df$variable == "ptQuery",]
  Global_node_mapping_df_ref$mean_expression <- ref_colMeans[paste0(Global_node_mapping_df_ref$variable,
                                                                    "_", Global_node_mapping_df_ref$value)]
  Global_node_mapping_df_query$mean_expression <- query_colMeans[paste0(Global_node_mapping_df_query$variable,
                                                                        "_", Global_node_mapping_df_query$value)]
  Global_node_mapping_df_ref <- as.data.frame(cbind(Global_node_mapping_df_ref,
                                                    metaNodesMap_Ref_exp[paste0(Global_node_mapping_df_ref$variable,
                                                                                "_", Global_node_mapping_df_ref$value),]))
  Global_node_mapping_df_query <- as.data.frame(cbind(Global_node_mapping_df_query,
                                                      metaNodesMap_Query_exp[paste0(Global_node_mapping_df_query$variable,
                                                                                    "_", Global_node_mapping_df_query$value),]))
  Global_node_mapping_df_2 <- as.data.frame(rbind(Global_node_mapping_df_ref,
                                                  Global_node_mapping_df_query))
  View(Global_node_mapping)
  sharedMarkers_cor <- c()
  for (i in sharedMarkers) {
    temp_ref_exp <- metaNodesMap_Ref_exp[paste0("ptRef_", Global_node_mapping$ptRef),i]
    temp_query_exp <- metaNodesMap_Query_exp[paste0("ptQuery_", Global_node_mapping$ptQuery),i]
    temp_cor_test <- cor.test(temp_ref_exp,
                              temp_query_exp)
    sharedMarkers_cor <- as.data.frame(rbind(sharedMarkers_cor,
                                             data.frame(Gene = i,
                                                        Cor = temp_cor_test$estimate,
                                                        P = temp_cor_test$p.value)))
  }
  
  sharedMarkers_cor$FDR <- p.adjust(sharedMarkers_cor$P,
                                    method = "fdr")
  ggplot(data = sharedMarkers_cor,
         aes(x = Cor)) +
    geom_histogram()
  
  selected_sharedMarkers_cor <- sharedMarkers_cor[abs(sharedMarkers_cor$Cor) > 0.8,]
  selected_sharedMarkers_cor_positive <- selected_sharedMarkers_cor[selected_sharedMarkers_cor$Cor > 0,]
  selected_sharedMarkers_cor_negative <- selected_sharedMarkers_cor[selected_sharedMarkers_cor$Cor < 0,]
  
  lines_cosine <- function(exp) {
    x <- 1:nrow(exp)
    y1 <- seq(from = 0, to = 1,
              length.out = length(x)) # 单调递增
    y2 <- seq(from = 1, to = 0,
              length.out = length(x)) # 单调递减
    y3 <- ((1/100)*x-1)^2 # 先减后增
    y4 <- 1-y3 # 先增后减
    cosine_df <- c()
    for (i in colnames(exp)) {
      cosine_1 <- lsa::cosine(cbind(exp[,i], y1))[1,2]
      cosine_2 <- lsa::cosine(cbind(exp[,i], y2))[1,2]
      cosine_3 <- lsa::cosine(cbind(exp[,i], y3))[1,2]
      cosine_4 <- lsa::cosine(cbind(exp[,i], y4))[1,2]
      
      cosine_df <- as.data.frame(rbind(cosine_df,
                                       data.frame(Gene = i,
                                                  y1 = cosine_1,
                                                  y2 = cosine_2,
                                                  y3 = cosine_3,
                                                  y4 = cosine_4)))
    }
    return(cosine_df)
  }
  selected_sharedMarkers_cor_positive <- selected_sharedMarkers_cor_positive[order(selected_sharedMarkers_cor_positive$Cor,
                                                                                   decreasing = T),]
  metaNodesMap_Ref_exp_positive <- metaNodesMap_Ref_exp[,selected_sharedMarkers_cor_positive$Gene]
  metaNodesMap_Query_exp_positive <- metaNodesMap_Query_exp[,selected_sharedMarkers_cor_positive$Gene]
  cosine_df <- lines_cosine(metaNodesMap_Ref_exp_positive)
  Type1 <- cosine_df$Gene[cosine_df$y1 >= 0.9 & cosine_df$y2 < 0.9]
  Type1 <- unlist(apply(metaNodesMap_Ref_exp_positive[,Type1,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type1 <- Type1[order(Type1, decreasing = F)]
  Type2 <- cosine_df$Gene[cosine_df$y1 < 0.9 & cosine_df$y2 >= 0.9]
  Type2 <- unlist(apply(metaNodesMap_Ref_exp_positive[,Type2,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type2 <- Type2[order(Type2, decreasing = F)]
  
  col_fun1 <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  col_anno_ref <- HeatmapAnnotation(
    Ref_Celltype = col_celltype,
    col = list(Ref_Celltype = RNA_Celltype_color)
  )
  col_anno_query <- HeatmapAnnotation(
    Query_Celltype = row_celltype,
    col = list(Query_Celltype = Discs_sub_color)
  )
  row_anno_gene <- rowAnnotation(
    Gene_type = c(rep("Type2", length(Type2)),
                  rep("Type1", length(Type1))),
    col = list(Gene_type = c("Type1" = "#F4A460",
                             "Type2" = "#1874CD"))
  )
  h1 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Ref_exp_positive[,c(names(Type2), names(Type1))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun1(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof RefTraj",
                                top_annotation = col_anno_ref,
                                left_annotation = row_anno_gene)
  col_fun2 <- circlize::colorRamp2(c(0, 0.5, 1), c("#009ACD", "#FFFACD", "#551A8B"))
  h2 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Query_exp_positive[,c(names(Type2), names(Type1))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun2(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof QueryTraj",
                                top_annotation = col_anno_query,
  )
  pdf("./轨迹比较/CellAlign_3.gene_heatmap_positive_Antler_Discs.pdf",
      width = 9, height = 6)
  h1 + h2
  dev.off()
  
  {
    ggplot() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.title = element_text(colour = "black", size = 10),
            legend.text = element_text(colour = "black", size = 10),
            axis.text = element_text(colour = "black", size = 10),
            axis.title = element_text(colour = "black", size = 10)) +
      geom_line(data = Global_node_mapping_df_2,
                aes(x = value, y = Type1_mean, group = rowIndex),
                color = "gray") +
      geom_point(data = Global_node_mapping_df_2,
                 aes(x = value, y = Type1_mean, color = variable)) +
      labs(x = "Pseudotime", y = "Mean interpolated expression of type1 gene set",
           color = "Trajectory") +
      scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                    "ptRef" = "#3A5FCD"))
  }
  
  
  Global_node_mapping_df_2$Type1_mean <- rowMeans(Global_node_mapping_df_2[,names(Type1)])
  Global_node_mapping_df_2$Type2_mean <- rowMeans(Global_node_mapping_df_2[,names(Type2)])
  
  pdf("./轨迹比较/CellAlign_4.positive_type2_mean_exp_Antler_Discs.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type2_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type2_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type2 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  pdf("./轨迹比较/CellAlign_4.positive_type1_mean_exp_Antler_Discs.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type1_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type1_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type1 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  
  metaNodesMap_Ref_exp_negative <- metaNodesMap_Ref_exp[,selected_sharedMarkers_cor_negative$Gene]
  metaNodesMap_Query_exp_negative <- metaNodesMap_Query_exp[,selected_sharedMarkers_cor_negative$Gene]
  cosine_df <- lines_cosine(metaNodesMap_Ref_exp_negative)
  Type3 <- cosine_df$Gene[cosine_df$y1 >= 0.9 & cosine_df$y2 < 0.9]
  Type3 <- unlist(apply(metaNodesMap_Ref_exp_negative[,Type3,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type3 <- Type3[order(Type3, decreasing = F)]
  Type4 <- cosine_df$Gene[cosine_df$y1 < 0.9 & cosine_df$y2 >= 0.9]
  Type4 <- unlist(apply(metaNodesMap_Ref_exp_negative[,Type4,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type4 <- Type4[order(Type4, decreasing = F)]
  
  col_fun1 <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  col_anno_ref <- HeatmapAnnotation(
    Ref_Celltype = col_celltype,
    col = list(Ref_Celltype = RNA_Celltype_color)
  )
  col_anno_query <- HeatmapAnnotation(
    Query_Celltype = row_celltype,
    col = list(Query_Celltype = Discs_sub_color)
  )
  row_anno_gene <- rowAnnotation(
    Gene_type = c(rep("Type4", length(Type4)),
                  rep("Type3", length(Type3))),
    col = list(Gene_type = c("Type3" = "#836FFF",
                             "Type4" = "#D02090"))
  )
  h3 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Ref_exp_negative[,c(names(Type4), names(Type3))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun1(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof RefTraj",
                                top_annotation = col_anno_ref,
                                left_annotation = row_anno_gene)
  col_fun2 <- circlize::colorRamp2(c(0, 0.5, 1), c("#009ACD", "#FFFACD", "#551A8B"))
  h4 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Query_exp_negative[,c(names(Type4), names(Type3))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun2(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof QueryTraj",
                                top_annotation = col_anno_query,
  )
  pdf("./轨迹比较/CellAlign_5.gene_heatmap_negative_Antler_Discs.pdf",
      width = 9, height = 5)
  h3 + h4
  dev.off()
  
  Global_node_mapping_df_2$Type3_mean <- rowMeans(Global_node_mapping_df_2[,names(Type3)])
  Global_node_mapping_df_2$Type4_mean <- rowMeans(Global_node_mapping_df_2[,names(Type4)])
  
  pdf("./轨迹比较/CellAlign_6.negative_type3_mean_exp_Antler_Discs.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type3_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type3_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type3 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  pdf("./轨迹比较/CellAlign_6.positive_type4_mean_exp_Antler_Discs.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type4_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type4_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type4 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  
  library(clusterProfiler)
  GeneList <- read.gmt("msigdb.v2023.2.Hs.symbols.gmt")
  GeneList <- GeneList[GeneList$gene %in% geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]],]
  Type1_enrich <- enricher(gene = names(Type1), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type1_enrich@result <- Type1_enrich@result[Type1_enrich@result$pvalue < 0.05,]
  Type1_enrich_result <- Type1_enrich@result
  Type1_enrich_result$ID <- unlist(lapply(Type1_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type1_enrich_result$Description <- unlist(lapply(Type1_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type1_enrich_result) <- Type1_enrich_result$ID
  Type1_enrich@result <- Type1_enrich_result
  openxlsx::write.xlsx(Type1_enrich@result,
                       "./轨迹比较/CellAlign_7.type1_genes_enrichment_Antler_Discs.xlsx")
  Type2_enrich <- enricher(gene = names(Type2), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type2_enrich@result <- Type2_enrich@result[Type2_enrich@result$pvalue < 0.05,]
  Type2_enrich_result <- Type2_enrich@result
  Type2_enrich_result$ID <- unlist(lapply(Type2_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type2_enrich_result$Description <- unlist(lapply(Type2_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type2_enrich_result) <- Type2_enrich_result$ID
  Type2_enrich@result <- Type2_enrich_result
  openxlsx::write.xlsx(Type2_enrich@result,
                       "./轨迹比较/CellAlign_7.type2_genes_enrichment_Antler_Discs.xlsx")
  Type3_enrich <- enricher(gene = names(Type3), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type3_enrich@result <- Type3_enrich@result[Type3_enrich@result$pvalue < 0.05,]
  Type3_enrich_result <- Type3_enrich@result
  Type3_enrich_result$ID <- unlist(lapply(Type3_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type3_enrich_result$Description <- unlist(lapply(Type3_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type3_enrich_result) <- Type3_enrich_result$ID
  Type3_enrich@result <- Type3_enrich_result
  openxlsx::write.xlsx(Type3_enrich@result,
                       "./轨迹比较/CellAlign_7.type3_genes_enrichment_Antler_Discs.xlsx")
  
  Type4_enrich <- enricher(gene = names(Type4), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type4_enrich@result <- Type4_enrich@result[Type4_enrich@result$pvalue < 0.05,]
  Type4_enrich_result <- Type4_enrich@result
  Type4_enrich_result$ID <- unlist(lapply(Type4_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type4_enrich_result$Description <- unlist(lapply(Type4_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type4_enrich_result) <- Type4_enrich_result$ID
  Type4_enrich@result <- Type4_enrich_result
  openxlsx::write.xlsx(Type4_enrich@result,
                       "./轨迹比较/CellAlign_7.type4_genes_enrichment_Antler_Discs.xlsx")
  
  Gene_type_df <- data.frame(Type = c(rep("Type1", length(Type1)),
                                      rep("Type2", length(Type2)),
                                      rep("Type3", length(Type3)),
                                      rep("Type4", length(Type4))),
                             Gene = c(names(Type1), names(Type2),
                                      names(Type3), names(Type4)))
  Gene_type_df_list <- split.data.frame(Gene_type_df, f = list(Gene_type_df$Type))
  openxlsx::write.xlsx(Gene_type_df_list,
                       "./轨迹比较/CellAlign_8.type_genes_list_Antler_Discs.xlsx")
  
}

###### RNA & Tumor
{
  interGlobal_RNA_sub <- readRDS("./轨迹比较/interGlobal_RNA_sub.rds")
  interGlobal_RNA_sub$interpolatedVals <- interGlobal_RNA_sub$interpolatedVals[rowSums(interGlobal_RNA_sub$interpolatedVals) > 0,]
  interGlobal_RNA_sub$error <- interGlobal_RNA_sub$error[rownames(interGlobal_RNA_sub$interpolatedVals),]
  interGlobal_Tumor_sub <- readRDS("./轨迹比较/interGlobal_Tumor_sub.rds")
  interGlobal_Tumor_sub$interpolatedVals <- interGlobal_Tumor_sub$interpolatedVals[rowSums(interGlobal_Tumor_sub$interpolatedVals) > 0,]
  interGlobal_Tumor_sub$error <- interGlobal_Tumor_sub$error[rownames(interGlobal_Tumor_sub$interpolatedVals),]
  
  sharedMarkers <- intersect(rownames(interGlobal_Tumor_sub$interpolatedVals),
                             rownames(interGlobal_RNA_sub$interpolatedVals))
  length(sharedMarkers)
  
  RNA_Tumor_CellAlign_Global <- CellAlign_Global_fun(interGlobal_ref = interGlobal_RNA_sub,
                                                     interGlobal_query = interGlobal_Tumor_sub,
                                                     ref_Traj = RNA_singlet_sub$pseudotime,
                                                     query_Traj = Tumor_sub$pseudotime,
                                                     numPts = 200, sharedMarkers = sharedMarkers)
  
  plotAlign(RNA_Tumor_CellAlign_Global$ref_query_alignment)
  
  plotMapping(RNA_Tumor_CellAlign_Global$ref_query_mapping)
  
  Global_node_mapping <- RNA_Tumor_CellAlign_Global[["ref_query_mapping"]][["metaNodesPt"]]
  Global_node_mapping <- Global_node_mapping[,-5]
  Global_node_mapping$index <-paste0(Global_node_mapping$ptRef,
                                     "_", Global_node_mapping$ptQuery)
  
  RNA_Tumor_CellAlign_Local <- CellAlign_Local_fun(interGlobal_ref = interGlobal_RNA_sub,
                                                   interGlobal_query = interGlobal_Tumor_sub,
                                                   ref_Traj = RNA_singlet_sub$pseudotime,
                                                   query_Traj = Tumor_sub$pseudotime,
                                                   numPts = 200, Thresh = 0.2,
                                                   sharedMarkers = sharedMarkers)
  
  plotAlign(RNA_Tumor_CellAlign_Local$ref_query_alignment)
  
  plotMapping(RNA_Tumor_CellAlign_Local$ref_query_mapping)
  
  Local_node_mapping <- RNA_Tumor_CellAlign_Local[["ref_query_mapping"]][["metaNodesPt"]]
  Local_node_mapping$index <- paste0(Local_node_mapping$ptRef,
                                     "_", Local_node_mapping$ptQuery)
  sum(Local_node_mapping$index %in% Global_node_mapping$index)
  Global_node_mapping$Mapping <- ifelse(Global_node_mapping$index %in% Local_node_mapping$index,
                                        "Conserved", "Not")
  Global_node_mapping$rowIndex <- 1:nrow(Global_node_mapping)
  
  Global_ref_node_cells <- Global_node_mapping[,c("metaNodeRef",
                                                  "ptRef")]
  Global_ref_node_cells <- Global_ref_node_cells[!duplicated(Global_ref_node_cells),]
  Global_ref_node_cells$ptRef <- paste0("ptRef", "_", Global_ref_node_cells$ptRef)
  rownames(Global_ref_node_cells) <- Global_ref_node_cells$ptRef
  
  Global_query_node_cells <- Global_node_mapping[,c("metaNodeQuery",
                                                    "ptQuery")]
  Global_query_node_cells <- Global_query_node_cells[!duplicated(Global_query_node_cells),]
  Global_query_node_cells$ptQuery <- paste0("ptQuery", "_", Global_query_node_cells$ptQuery)
  rownames(Global_query_node_cells) <- Global_query_node_cells$ptQuery
  node_cells <- data.frame(node = c(Global_ref_node_cells$ptRef,
                                    Global_query_node_cells$ptQuery),
                           cells = c(Global_ref_node_cells$metaNodeRef,
                                     Global_query_node_cells$metaNodeQuery))
  rownames(node_cells) <- node_cells$node
  
  Global_node_mapping_df <- melt(Global_node_mapping[, c("ptQuery", "ptRef", "rowIndex")],
                                 id.vars = c("rowIndex"))
  Global_node_mapping_df$cells <- paste0(Global_node_mapping_df$variable,
                                         "_", Global_node_mapping_df$value)
  Global_node_mapping_df$cells <- node_cells[Global_node_mapping_df$cells, "cells"]
  cells <- c(RNA_Tumor_CellAlign_Global$ref_query_mapping$refAssign,
             RNA_Tumor_CellAlign_Global$ref_query_mapping$queryAssign)
  main_celltype <- c()
  for (i in 1:nrow(Global_node_mapping_df)) {
    if (Global_node_mapping_df$variable[i] == "ptQuery") {
      temp <- table(Tumor_sub@meta.data[cells[[Global_node_mapping_df$cells[i]]], "celltype"])
      temp_celltype <- names(temp)[which.max(temp)]
      main_celltype <- c(main_celltype,
                         temp_celltype)
    }
    if (Global_node_mapping_df$variable[i] == "ptRef") {
      temp <- table(RNA_singlet_sub@meta.data[cells[[Global_node_mapping_df$cells[i]]], "Celltype_rename"])
      temp_celltype <- names(temp)[which.max(temp)]
      main_celltype <- c(main_celltype,
                         temp_celltype)
    }
  }
  Global_node_mapping_df$main_celltype <- main_celltype
  Global_node_mapping_df_Conserved <- melt(Global_node_mapping[Global_node_mapping$Mapping == "Conserved",
                                                               c("ptQuery", "ptRef", "rowIndex")],
                                           id.vars = c("rowIndex"))
  Global_node_mapping_df_Not <- melt(Global_node_mapping[Global_node_mapping$Mapping == "Not",
                                                         c("ptQuery", "ptRef", "rowIndex")],
                                     id.vars = c("rowIndex"))
  pdf("./轨迹比较/CellAlign_1.point_alignment_Antler_Tumor.pdf",
      height = 2.5, width = 7.5)
  ggplot() + 
    theme_bw() +
    geom_line(data = Global_node_mapping_df_Conserved,
              aes(x = variable, y = value, group = rowIndex),
              color = "red") +
    geom_line(data = Global_node_mapping_df_Not,
              aes(x = variable, y = value, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df,
               aes(x = variable, y = value,
                   colour = main_celltype),
               shape = 15) +
    scale_color_manual(values = c(RNA_Celltype_color,
                                  Tumor_sub_color)) +
    coord_flip() + ggtitle("Ref_Antler vs. Query_Tumor") +
    labs(color = "Celltype") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          legend.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  dev.off()
  
  costMat <- RNA_Tumor_CellAlign_Global$ref_query_alignment$localCostMatrix
  costMat <- t(apply(costMat, 1, function(x) {
    return(as.numeric(x))
  })) # 列 ref, 行 query
  linearInd <- sub2ind(nrow(costMat), RNA_Tumor_CellAlign_Global$ref_query_alignment$align[[1]]$index1, 
                       RNA_Tumor_CellAlign_Global$ref_query_alignment$align[[1]]$index2)
  costMat[linearInd] <- NA
  costMat <- as.data.frame(costMat)
  
  metaNodesMap <- metaNodesMapping(RNA_Tumor_CellAlign_Global$ref_query_alignment,
                                   intTrajQuery = interGlobal_Tumor_sub$traj, 
                                   realTrajQuery = Tumor_sub$pseudotime,
                                   intTrajRef = interGlobal_RNA_sub$traj,
                                   realTrajRef = RNA_singlet_sub$pseudotime)
  metaNodesMap_Query <- metaNodesMap[!duplicated(metaNodesMap$queryInt),]
  metaNodesMap_Query <- metaNodesMap_Query[order(metaNodesMap_Query$queryInt,
                                                 decreasing = F),]
  metaNodesMap_Ref <- metaNodesMap[!duplicated(metaNodesMap$refInt),]
  metaNodesMap_Ref <- metaNodesMap_Ref[order(metaNodesMap_Ref$refInt,
                                             decreasing = F),]
  metaNodesMap_Ref$ptRef <- paste0("ptRef_",
                                   metaNodesMap_Ref$ptRef)
  metaNodesMap_Ref <- metaNodesMap_Ref[!duplicated(metaNodesMap_Ref$ptRef),]
  metaNodesMap_Query$ptQuery <- paste0("ptQuery_",
                                       metaNodesMap_Query$ptQuery)
  metaNodesMap_Query <- metaNodesMap_Query[!duplicated(metaNodesMap_Query$ptQuery),]
  
  costMat <- costMat[,metaNodesMap_Ref$refInt]
  costMat <- costMat[metaNodesMap_Query$queryInt,]
  colnames(costMat) <- metaNodesMap_Ref$ptRef # 列 ref
  rownames(costMat) <- metaNodesMap_Query$ptQuery # 行 query
  
  col_celltype <- c()
  for (i in 1:ncol(costMat)) {
    temp <- table(RNA_singlet_sub@meta.data[cells[[metaNodesMap_Ref$cellIDRef[i]]], "Celltype_rename"])
    temp_celltype <- names(temp)[which.max(temp)]
    col_celltype <- c(col_celltype,
                      temp_celltype)
  }
  row_celltype <- c()
  for (i in 1:nrow(costMat)) {
    temp <- table(Tumor_sub@meta.data[cells[[metaNodesMap_Query$cellIDQuery[i]]], "celltype"])
    temp_celltype <- names(temp)[which.max(temp)]
    row_celltype <- c(row_celltype,
                      temp_celltype)
  }
  
  Global_node_mapping_df_Conserved$node <- paste0(Global_node_mapping_df_Conserved$variable, "_",
                                                  Global_node_mapping_df_Conserved$value)
  annotCols <- data.frame(Ref_celltype = col_celltype,
                          Local_alignment = "Not conserved",
                          row.names = colnames(costMat))
  annotCols[Global_node_mapping_df_Conserved$node[Global_node_mapping_df_Conserved$variable == "ptRef"], "Local_alignment"] <- "Conserved"
  annotRows <- data.frame(Query_celltype = row_celltype,
                          Local_alignment = "Not conserved",
                          row.names = rownames(costMat))
  annotRows[Global_node_mapping_df_Conserved$node[Global_node_mapping_df_Conserved$variable == "ptQuery"], "Local_alignment"] <- "Conserved"
  
  annotColors <- list(Ref_celltype = RNA_Celltype_color,
                      Query_celltype = Tumor_sub_color,
                      Local_alignment = c("Not conserved" = "gray", "Conserved" = "red"))
  pdf("./轨迹比较/CellAlign_2.alignment_pheatmap_Antler_Tumor.pdf",
      width = 8, height = 8)
  pheatmap(costMat, na_col = "purple",
           show_rownames = F, show_colnames = F,
           annotation_row = annotRows,
           annotation_col = annotCols,
           annotation_colors = annotColors,
           cluster_rows = F, cluster_cols = F,
           cellheight = 1.5, cellwidth = 1.5)
  dev.off()
  
  interScaledLocal_RNA_sub <- cellAlign::scaleInterpolate(interGlobal_RNA_sub)
  interScaledLocal_Tumor_sub <- cellAlign::scaleInterpolate(interGlobal_Tumor_sub)
  ref_colMeans <- colMeans(interScaledLocal_RNA_sub$scaledData[sharedMarkers,])
  ref_colMeans <- ref_colMeans[metaNodesMap_Ref$refInt]
  names(ref_colMeans) <- colnames(costMat)
  query_colMeans <- colMeans(interScaledLocal_Tumor_sub$scaledData[sharedMarkers,])
  query_colMeans <- query_colMeans[metaNodesMap_Query$queryInt]
  names(query_colMeans) <- rownames(costMat)
  
  metaNodesMap_Ref_exp <- interScaledLocal_RNA_sub$scaledData[sharedMarkers, metaNodesMap_Ref$refInt]
  metaNodesMap_Ref_exp <- data.frame(t(metaNodesMap_Ref_exp),
                                     check.rows = F, check.names = F)
  rownames(metaNodesMap_Ref_exp) <- colnames(costMat)
  metaNodesMap_Query_exp <- interScaledLocal_Tumor_sub$scaledData[sharedMarkers, metaNodesMap_Query$queryInt]
  metaNodesMap_Query_exp <- data.frame(t(metaNodesMap_Query_exp),
                                       check.rows = F, check.names = F)
  rownames(metaNodesMap_Query_exp) <- rownames(costMat)
  
  Global_node_mapping_df_ref <- Global_node_mapping_df[Global_node_mapping_df$variable == "ptRef",]
  Global_node_mapping_df_query <- Global_node_mapping_df[Global_node_mapping_df$variable == "ptQuery",]
  Global_node_mapping_df_ref$mean_expression <- ref_colMeans[paste0(Global_node_mapping_df_ref$variable,
                                                                    "_", Global_node_mapping_df_ref$value)]
  Global_node_mapping_df_query$mean_expression <- query_colMeans[paste0(Global_node_mapping_df_query$variable,
                                                                        "_", Global_node_mapping_df_query$value)]
  Global_node_mapping_df_ref <- as.data.frame(cbind(Global_node_mapping_df_ref,
                                                    metaNodesMap_Ref_exp[paste0(Global_node_mapping_df_ref$variable,
                                                                                "_", Global_node_mapping_df_ref$value),]))
  Global_node_mapping_df_query <- as.data.frame(cbind(Global_node_mapping_df_query,
                                                      metaNodesMap_Query_exp[paste0(Global_node_mapping_df_query$variable,
                                                                                    "_", Global_node_mapping_df_query$value),]))
  Global_node_mapping_df_2 <- as.data.frame(rbind(Global_node_mapping_df_ref,
                                                  Global_node_mapping_df_query))
  View(Global_node_mapping)
  sharedMarkers_cor <- c()
  for (i in sharedMarkers) {
    temp_ref_exp <- metaNodesMap_Ref_exp[paste0("ptRef_", Global_node_mapping$ptRef),i]
    temp_query_exp <- metaNodesMap_Query_exp[paste0("ptQuery_", Global_node_mapping$ptQuery),i]
    temp_cor_test <- cor.test(temp_ref_exp,
                              temp_query_exp)
    sharedMarkers_cor <- as.data.frame(rbind(sharedMarkers_cor,
                                             data.frame(Gene = i,
                                                        Cor = temp_cor_test$estimate,
                                                        P = temp_cor_test$p.value)))
  }
  
  sharedMarkers_cor$FDR <- p.adjust(sharedMarkers_cor$P,
                                    method = "fdr")
  ggplot(data = sharedMarkers_cor,
         aes(x = Cor)) +
    geom_histogram()
  
  selected_sharedMarkers_cor <- sharedMarkers_cor[abs(sharedMarkers_cor$Cor) > 0.8,]
  selected_sharedMarkers_cor_positive <- selected_sharedMarkers_cor[selected_sharedMarkers_cor$Cor > 0,]
  selected_sharedMarkers_cor_negative <- selected_sharedMarkers_cor[selected_sharedMarkers_cor$Cor < 0,]
  
  lines_cosine <- function(exp) {
    x <- 1:nrow(exp)
    y1 <- seq(from = 0, to = 1,
              length.out = length(x)) # 单调递增
    y2 <- seq(from = 1, to = 0,
              length.out = length(x)) # 单调递减
    y3 <- ((1/100)*x-1)^2 # 先减后增
    y4 <- 1-y3 # 先增后减
    cosine_df <- c()
    for (i in colnames(exp)) {
      cosine_1 <- lsa::cosine(cbind(exp[,i], y1))[1,2]
      cosine_2 <- lsa::cosine(cbind(exp[,i], y2))[1,2]
      cosine_3 <- lsa::cosine(cbind(exp[,i], y3))[1,2]
      cosine_4 <- lsa::cosine(cbind(exp[,i], y4))[1,2]
      
      cosine_df <- as.data.frame(rbind(cosine_df,
                                       data.frame(Gene = i,
                                                  y1 = cosine_1,
                                                  y2 = cosine_2,
                                                  y3 = cosine_3,
                                                  y4 = cosine_4)))
    }
    return(cosine_df)
  }
  selected_sharedMarkers_cor_positive <- selected_sharedMarkers_cor_positive[order(selected_sharedMarkers_cor_positive$Cor,
                                                                                   decreasing = T),]
  metaNodesMap_Ref_exp_positive <- metaNodesMap_Ref_exp[,selected_sharedMarkers_cor_positive$Gene]
  metaNodesMap_Query_exp_positive <- metaNodesMap_Query_exp[,selected_sharedMarkers_cor_positive$Gene]
  cosine_df <- lines_cosine(metaNodesMap_Ref_exp_positive)
  Type1 <- cosine_df$Gene[cosine_df$y1 >= 0.9 & cosine_df$y2 < 0.9]
  Type1 <- unlist(apply(metaNodesMap_Ref_exp_positive[,Type1,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type1 <- Type1[order(Type1, decreasing = F)]
  Type2 <- cosine_df$Gene[cosine_df$y1 < 0.9 & cosine_df$y2 >= 0.9]
  Type2 <- unlist(apply(metaNodesMap_Ref_exp_positive[,Type2,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type2 <- Type2[order(Type2, decreasing = F)]
  
  library(ComplexHeatmap)
  col_fun1 <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  col_anno_ref <- HeatmapAnnotation(
    Ref_Celltype = col_celltype,
    col = list(Ref_Celltype = RNA_Celltype_color)
  )
  col_anno_query <- HeatmapAnnotation(
    Query_Celltype = row_celltype,
    col = list(Query_Celltype = Tumor_sub_color)
  )
  row_anno_gene <- rowAnnotation(
    Gene_type = c(rep("Type2", length(Type2)),
                  rep("Type1", length(Type1))),
    col = list(Gene_type = c("Type1" = "#F4A460",
                             "Type2" = "#1874CD"))
  )
  h1 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Ref_exp_positive[,c(names(Type2), names(Type1))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun1(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof RefTraj",
                                top_annotation = col_anno_ref,
                                left_annotation = row_anno_gene)
  col_fun2 <- circlize::colorRamp2(c(0, 0.5, 1), c("#009ACD", "#FFFACD", "#551A8B"))
  h2 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Query_exp_positive[,c(names(Type2), names(Type1))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun2(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof QueryTraj",
                                top_annotation = col_anno_query,
  )
  pdf("./轨迹比较/CellAlign_3.gene_heatmap_positive_Antler_Tumor.pdf",
      width = 9, height = 6)
  h1 + h2
  dev.off()
  
  {
    ggplot() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.title = element_text(colour = "black", size = 10),
            legend.text = element_text(colour = "black", size = 10),
            axis.text = element_text(colour = "black", size = 10),
            axis.title = element_text(colour = "black", size = 10)) +
      geom_line(data = Global_node_mapping_df_2,
                aes(x = value, y = Type1_mean, group = rowIndex),
                color = "gray") +
      geom_point(data = Global_node_mapping_df_2,
                 aes(x = value, y = Type1_mean, color = variable)) +
      labs(x = "Pseudotime", y = "Mean interpolated expression of type1 gene set",
           color = "Trajectory") +
      scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                    "ptRef" = "#3A5FCD"))
  }
  
  
  Global_node_mapping_df_2$Type1_mean <- rowMeans(Global_node_mapping_df_2[,names(Type1)])
  Global_node_mapping_df_2$Type2_mean <- rowMeans(Global_node_mapping_df_2[,names(Type2)])
  
  pdf("./轨迹比较/CellAlign_4.positive_type2_mean_exp_Antler_Tumor.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type2_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type2_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type2 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  pdf("./轨迹比较/CellAlign_4.positive_type1_mean_exp_Antler_Tumor.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type1_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type1_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type1 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  
  metaNodesMap_Ref_exp_negative <- metaNodesMap_Ref_exp[,selected_sharedMarkers_cor_negative$Gene]
  metaNodesMap_Query_exp_negative <- metaNodesMap_Query_exp[,selected_sharedMarkers_cor_negative$Gene]
  cosine_df <- lines_cosine(metaNodesMap_Ref_exp_negative)
  Type3 <- cosine_df$Gene[cosine_df$y1 >= 0.9 & cosine_df$y2 < 0.9]
  Type3 <- unlist(apply(metaNodesMap_Ref_exp_negative[,Type3,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type3 <- Type3[order(Type3, decreasing = F)]
  Type4 <- cosine_df$Gene[cosine_df$y1 < 0.9 & cosine_df$y2 >= 0.9]
  Type4 <- unlist(apply(metaNodesMap_Ref_exp_negative[,Type4,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type4 <- Type4[order(Type4, decreasing = F)]
  
  col_fun1 <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  col_anno_ref <- HeatmapAnnotation(
    Ref_Celltype = col_celltype,
    col = list(Ref_Celltype = RNA_Celltype_color)
  )
  col_anno_query <- HeatmapAnnotation(
    Query_Celltype = row_celltype,
    col = list(Query_Celltype = Tumor_sub_color)
  )
  row_anno_gene <- rowAnnotation(
    Gene_type = c(rep("Type4", length(Type4)),
                  rep("Type3", length(Type3))),
    col = list(Gene_type = c("Type3" = "#836FFF",
                             "Type4" = "#D02090"))
  )
  h3 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Ref_exp_negative[,c(names(Type4), names(Type3))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun1(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof RefTraj",
                                top_annotation = col_anno_ref,
                                left_annotation = row_anno_gene)
  col_fun2 <- circlize::colorRamp2(c(0, 0.5, 1), c("#009ACD", "#FFFACD", "#551A8B"))
  h4 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Query_exp_negative[,c(names(Type4), names(Type3))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun2(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof QueryTraj",
                                top_annotation = col_anno_query,
  )
  pdf("./轨迹比较/CellAlign_5.gene_heatmap_negative_Antler_Tumor.pdf",
      width = 9, height = 5)
  h3 + h4
  dev.off()
  
  Global_node_mapping_df_2$Type3_mean <- rowMeans(Global_node_mapping_df_2[,names(Type3)])
  Global_node_mapping_df_2$Type4_mean <- rowMeans(Global_node_mapping_df_2[,names(Type4)])
  
  pdf("./轨迹比较/CellAlign_6.negative_type3_mean_exp_Antler_Tumor.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type3_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type3_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type3 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  pdf("./轨迹比较/CellAlign_6.positive_type4_mean_exp_Antler_Tumor.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type4_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type4_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type4 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  
  library(clusterProfiler)
  GeneList <- read.gmt("msigdb.v2023.2.Hs.symbols.gmt")
  GeneList <- GeneList[GeneList$gene %in% geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]],]
  Type1_enrich <- enricher(gene = names(Type1), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type1_enrich@result <- Type1_enrich@result[Type1_enrich@result$pvalue < 0.05,]
  Type1_enrich_result <- Type1_enrich@result
  Type1_enrich_result$ID <- unlist(lapply(Type1_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type1_enrich_result$Description <- unlist(lapply(Type1_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type1_enrich_result) <- Type1_enrich_result$ID
  Type1_enrich@result <- Type1_enrich_result
  openxlsx::write.xlsx(Type1_enrich@result,
                       "./轨迹比较/CellAlign_7.type1_genes_enrichment_Antler_Tumor.xlsx")
  Type2_enrich <- enricher(gene = names(Type2), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type2_enrich@result <- Type2_enrich@result[Type2_enrich@result$pvalue < 0.05,]
  Type2_enrich_result <- Type2_enrich@result
  Type2_enrich_result$ID <- unlist(lapply(Type2_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type2_enrich_result$Description <- unlist(lapply(Type2_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type2_enrich_result) <- Type2_enrich_result$ID
  Type2_enrich@result <- Type2_enrich_result
  openxlsx::write.xlsx(Type2_enrich@result,
                       "./轨迹比较/CellAlign_7.type2_genes_enrichment_Antler_Tumor.xlsx")
  Type3_enrich <- enricher(gene = names(Type3), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type3_enrich@result <- Type3_enrich@result[Type3_enrich@result$pvalue < 0.05,]
  Type3_enrich_result <- Type3_enrich@result
  Type3_enrich_result$ID <- unlist(lapply(Type3_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type3_enrich_result$Description <- unlist(lapply(Type3_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type3_enrich_result) <- Type3_enrich_result$ID
  Type3_enrich@result <- Type3_enrich_result
  openxlsx::write.xlsx(Type3_enrich@result,
                       "./轨迹比较/CellAlign_7.type3_genes_enrichment_Antler_Tumor.xlsx")
  
  Type4_enrich <- enricher(gene = names(Type4), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type4_enrich@result <- Type4_enrich@result[Type4_enrich@result$pvalue < 0.05,]
  Type4_enrich_result <- Type4_enrich@result
  Type4_enrich_result$ID <- unlist(lapply(Type4_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type4_enrich_result$Description <- unlist(lapply(Type4_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type4_enrich_result) <- Type4_enrich_result$ID
  Type4_enrich@result <- Type4_enrich_result
  openxlsx::write.xlsx(Type4_enrich@result,
                       "./轨迹比较/CellAlign_7.type4_genes_enrichment_Antler_Tumor.xlsx")
  
  Gene_type_df <- data.frame(Type = c(rep("Type1", length(Type1)),
                                      rep("Type2", length(Type2)),
                                      rep("Type3", length(Type3)),
                                      rep("Type4", length(Type4))),
                             Gene = c(names(Type1), names(Type2),
                                      names(Type3), names(Type4)))
  Gene_type_df_list <- split.data.frame(Gene_type_df, f = list(Gene_type_df$Type))
  openxlsx::write.xlsx(Gene_type_df_list,
                       "./轨迹比较/CellAlign_8.type_genes_list_Antler_Tumor.xlsx")
  
}

###### RNA & Limb
{
  interGlobal_RNA_sub <- readRDS("./轨迹比较/interGlobal_RNA_sub.rds")
  interGlobal_RNA_sub$interpolatedVals <- interGlobal_RNA_sub$interpolatedVals[rowSums(interGlobal_RNA_sub$interpolatedVals) > 0,]
  interGlobal_RNA_sub$error <- interGlobal_RNA_sub$error[rownames(interGlobal_RNA_sub$interpolatedVals),]
  interGlobal_Limb_sub <- readRDS("./轨迹比较/interGlobal_Limb_sub.rds")
  interGlobal_Limb_sub$interpolatedVals <- interGlobal_Limb_sub$interpolatedVals[rowSums(interGlobal_Limb_sub$interpolatedVals) > 0,]
  interGlobal_Limb_sub$error <- interGlobal_Limb_sub$error[rownames(interGlobal_Limb_sub$interpolatedVals),]
  
  sharedMarkers <- intersect(rownames(interGlobal_Limb_sub$interpolatedVals),
                             rownames(interGlobal_RNA_sub$interpolatedVals))
  length(sharedMarkers)
  
  RNA_Limb_CellAlign_Global <- CellAlign_Global_fun(interGlobal_ref = interGlobal_RNA_sub,
                                                    interGlobal_query = interGlobal_Limb_sub,
                                                    ref_Traj = RNA_singlet_sub$pseudotime,
                                                    query_Traj = Limb_sub$pseudotime,
                                                    numPts = 200, sharedMarkers = sharedMarkers)
  
  plotAlign(RNA_Limb_CellAlign_Global$ref_query_alignment)
  
  plotMapping(RNA_Limb_CellAlign_Global$ref_query_mapping)
  
  Global_node_mapping <- RNA_Limb_CellAlign_Global[["ref_query_mapping"]][["metaNodesPt"]]
  Global_node_mapping <- Global_node_mapping[,-5]
  Global_node_mapping$index <-paste0(Global_node_mapping$ptRef,
                                     "_", Global_node_mapping$ptQuery)
  
  RNA_Limb_CellAlign_Local <- CellAlign_Local_fun(interGlobal_ref = interGlobal_RNA_sub,
                                                  interGlobal_query = interGlobal_Limb_sub,
                                                  ref_Traj = RNA_singlet_sub$pseudotime,
                                                  query_Traj = Limb_sub$pseudotime,
                                                  numPts = 200, Thresh = 0.2,
                                                  sharedMarkers = sharedMarkers)
  
  plotAlign(RNA_Limb_CellAlign_Local$ref_query_alignment)
  
  plotMapping(RNA_Limb_CellAlign_Local$ref_query_mapping)
  
  Local_node_mapping <- RNA_Limb_CellAlign_Local[["ref_query_mapping"]][["metaNodesPt"]]
  Local_node_mapping$index <- paste0(Local_node_mapping$ptRef,
                                     "_", Local_node_mapping$ptQuery)
  sum(Local_node_mapping$index %in% Global_node_mapping$index)
  Global_node_mapping$Mapping <- ifelse(Global_node_mapping$index %in% Local_node_mapping$index,
                                        "Conserved", "Not")
  Global_node_mapping$rowIndex <- 1:nrow(Global_node_mapping)
  
  Global_ref_node_cells <- Global_node_mapping[,c("metaNodeRef",
                                                  "ptRef")]
  Global_ref_node_cells <- Global_ref_node_cells[!duplicated(Global_ref_node_cells),]
  Global_ref_node_cells$ptRef <- paste0("ptRef", "_", Global_ref_node_cells$ptRef)
  rownames(Global_ref_node_cells) <- Global_ref_node_cells$ptRef
  
  Global_query_node_cells <- Global_node_mapping[,c("metaNodeQuery",
                                                    "ptQuery")]
  Global_query_node_cells <- Global_query_node_cells[!duplicated(Global_query_node_cells),]
  Global_query_node_cells$ptQuery <- paste0("ptQuery", "_", Global_query_node_cells$ptQuery)
  rownames(Global_query_node_cells) <- Global_query_node_cells$ptQuery
  node_cells <- data.frame(node = c(Global_ref_node_cells$ptRef,
                                    Global_query_node_cells$ptQuery),
                           cells = c(Global_ref_node_cells$metaNodeRef,
                                     Global_query_node_cells$metaNodeQuery))
  rownames(node_cells) <- node_cells$node
  
  Global_node_mapping_df <- melt(Global_node_mapping[, c("ptQuery", "ptRef", "rowIndex")],
                                 id.vars = c("rowIndex"))
  Global_node_mapping_df$cells <- paste0(Global_node_mapping_df$variable,
                                         "_", Global_node_mapping_df$value)
  Global_node_mapping_df$cells <- node_cells[Global_node_mapping_df$cells, "cells"]
  cells <- c(RNA_Limb_CellAlign_Global$ref_query_mapping$refAssign,
             RNA_Limb_CellAlign_Global$ref_query_mapping$queryAssign)
  main_celltype <- c()
  for (i in 1:nrow(Global_node_mapping_df)) {
    if (Global_node_mapping_df$variable[i] == "ptQuery") {
      temp <- table(Limb_sub@meta.data[cells[[Global_node_mapping_df$cells[i]]], "celltype"])
      temp_celltype <- names(temp)[which.max(temp)]
      main_celltype <- c(main_celltype,
                         temp_celltype)
    }
    if (Global_node_mapping_df$variable[i] == "ptRef") {
      temp <- table(RNA_singlet_sub@meta.data[cells[[Global_node_mapping_df$cells[i]]], "Celltype_rename"])
      temp_celltype <- names(temp)[which.max(temp)]
      main_celltype <- c(main_celltype,
                         temp_celltype)
    }
  }
  Global_node_mapping_df$main_celltype <- main_celltype
  Global_node_mapping_df_Conserved <- melt(Global_node_mapping[Global_node_mapping$Mapping == "Conserved",
                                                               c("ptQuery", "ptRef", "rowIndex")],
                                           id.vars = c("rowIndex"))
  Global_node_mapping_df_Not <- melt(Global_node_mapping[Global_node_mapping$Mapping == "Not",
                                                         c("ptQuery", "ptRef", "rowIndex")],
                                     id.vars = c("rowIndex"))
  pdf("./轨迹比较/CellAlign_1.point_alignment_Antler_Limb.pdf",
      height = 2.5, width = 7.5)
  ggplot() + 
    theme_bw() +
    geom_line(data = Global_node_mapping_df_Conserved,
              aes(x = variable, y = value, group = rowIndex),
              color = "red") +
    geom_line(data = Global_node_mapping_df_Not,
              aes(x = variable, y = value, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df,
               aes(x = variable, y = value,
                   colour = main_celltype),
               shape = 15) +
    scale_color_manual(values = c(RNA_Celltype_color,
                                  Limb_sub_color)) +
    coord_flip() + ggtitle("Ref_Antler vs. Query_Limb") +
    labs(color = "Celltype") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.ticks.y = element_line(colour = "black"),
          legend.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  dev.off()
  
  costMat <- RNA_Limb_CellAlign_Global$ref_query_alignment$localCostMatrix
  costMat <- t(apply(costMat, 1, function(x) {
    return(as.numeric(x))
  })) # 列 ref, 行 query
  linearInd <- sub2ind(nrow(costMat), RNA_Limb_CellAlign_Global$ref_query_alignment$align[[1]]$index1, 
                       RNA_Limb_CellAlign_Global$ref_query_alignment$align[[1]]$index2)
  costMat[linearInd] <- NA
  costMat <- as.data.frame(costMat)
  
  metaNodesMap <- metaNodesMapping(RNA_Limb_CellAlign_Global$ref_query_alignment,
                                   intTrajQuery = interGlobal_Limb_sub$traj, 
                                   realTrajQuery = Limb_sub$pseudotime,
                                   intTrajRef = interGlobal_RNA_sub$traj,
                                   realTrajRef = RNA_singlet_sub$pseudotime)
  metaNodesMap_Query <- metaNodesMap[!duplicated(metaNodesMap$queryInt),]
  metaNodesMap_Query <- metaNodesMap_Query[order(metaNodesMap_Query$queryInt,
                                                 decreasing = F),]
  metaNodesMap_Ref <- metaNodesMap[!duplicated(metaNodesMap$refInt),]
  metaNodesMap_Ref <- metaNodesMap_Ref[order(metaNodesMap_Ref$refInt,
                                             decreasing = F),]
  metaNodesMap_Ref$ptRef <- paste0("ptRef_",
                                   metaNodesMap_Ref$ptRef)
  metaNodesMap_Ref <- metaNodesMap_Ref[!duplicated(metaNodesMap_Ref$ptRef),]
  metaNodesMap_Query$ptQuery <- paste0("ptQuery_",
                                       metaNodesMap_Query$ptQuery)
  metaNodesMap_Query <- metaNodesMap_Query[!duplicated(metaNodesMap_Query$ptQuery),]
  
  costMat <- costMat[,metaNodesMap_Ref$refInt]
  costMat <- costMat[metaNodesMap_Query$queryInt,]
  colnames(costMat) <- metaNodesMap_Ref$ptRef # 列 ref
  rownames(costMat) <- metaNodesMap_Query$ptQuery # 行 query
  
  col_celltype <- c()
  for (i in 1:ncol(costMat)) {
    temp <- table(RNA_singlet_sub@meta.data[cells[[metaNodesMap_Ref$cellIDRef[i]]], "Celltype_rename"])
    temp_celltype <- names(temp)[which.max(temp)]
    col_celltype <- c(col_celltype,
                      temp_celltype)
  }
  row_celltype <- c()
  for (i in 1:nrow(costMat)) {
    temp <- table(Limb_sub@meta.data[cells[[metaNodesMap_Query$cellIDQuery[i]]], "celltype"])
    temp_celltype <- names(temp)[which.max(temp)]
    row_celltype <- c(row_celltype,
                      temp_celltype)
  }
  
  Global_node_mapping_df_Conserved$node <- paste0(Global_node_mapping_df_Conserved$variable, "_",
                                                  Global_node_mapping_df_Conserved$value)
  annotCols <- data.frame(Ref_celltype = col_celltype,
                          Local_alignment = "Not conserved",
                          row.names = colnames(costMat))
  annotCols[Global_node_mapping_df_Conserved$node[Global_node_mapping_df_Conserved$variable == "ptRef"], "Local_alignment"] <- "Conserved"
  annotRows <- data.frame(Query_celltype = row_celltype,
                          Local_alignment = "Not conserved",
                          row.names = rownames(costMat))
  annotRows[Global_node_mapping_df_Conserved$node[Global_node_mapping_df_Conserved$variable == "ptQuery"], "Local_alignment"] <- "Conserved"
  
  annotColors <- list(Ref_celltype = RNA_Celltype_color,
                      Query_celltype = Limb_sub_color,
                      Local_alignment = c("Not conserved" = "gray", "Conserved" = "red"))
  pdf("./轨迹比较/CellAlign_2.alignment_pheatmap_Antler_Limb.pdf",
      width = 8, height = 8)
  pheatmap(costMat, na_col = "purple",
           show_rownames = F, show_colnames = F,
           annotation_row = annotRows,
           annotation_col = annotCols,
           annotation_colors = annotColors,
           cluster_rows = F, cluster_cols = F,
           cellheight = 1.5, cellwidth = 1.5)
  dev.off()
  
  interScaledLocal_RNA_sub <- cellAlign::scaleInterpolate(interGlobal_RNA_sub)
  interScaledLocal_Limb_sub <- cellAlign::scaleInterpolate(interGlobal_Limb_sub)
  ref_colMeans <- colMeans(interScaledLocal_RNA_sub$scaledData[sharedMarkers,])
  ref_colMeans <- ref_colMeans[metaNodesMap_Ref$refInt]
  names(ref_colMeans) <- colnames(costMat)
  query_colMeans <- colMeans(interScaledLocal_Limb_sub$scaledData[sharedMarkers,])
  query_colMeans <- query_colMeans[metaNodesMap_Query$queryInt]
  names(query_colMeans) <- rownames(costMat)
  
  metaNodesMap_Ref_exp <- interScaledLocal_RNA_sub$scaledData[sharedMarkers, metaNodesMap_Ref$refInt]
  metaNodesMap_Ref_exp <- data.frame(t(metaNodesMap_Ref_exp),
                                     check.rows = F, check.names = F)
  rownames(metaNodesMap_Ref_exp) <- colnames(costMat)
  metaNodesMap_Query_exp <- interScaledLocal_Limb_sub$scaledData[sharedMarkers, metaNodesMap_Query$queryInt]
  metaNodesMap_Query_exp <- data.frame(t(metaNodesMap_Query_exp),
                                       check.rows = F, check.names = F)
  rownames(metaNodesMap_Query_exp) <- rownames(costMat)
  
  Global_node_mapping_df_ref <- Global_node_mapping_df[Global_node_mapping_df$variable == "ptRef",]
  Global_node_mapping_df_query <- Global_node_mapping_df[Global_node_mapping_df$variable == "ptQuery",]
  Global_node_mapping_df_ref$mean_expression <- ref_colMeans[paste0(Global_node_mapping_df_ref$variable,
                                                                    "_", Global_node_mapping_df_ref$value)]
  Global_node_mapping_df_query$mean_expression <- query_colMeans[paste0(Global_node_mapping_df_query$variable,
                                                                        "_", Global_node_mapping_df_query$value)]
  Global_node_mapping_df_ref <- as.data.frame(cbind(Global_node_mapping_df_ref,
                                                    metaNodesMap_Ref_exp[paste0(Global_node_mapping_df_ref$variable,
                                                                                "_", Global_node_mapping_df_ref$value),]))
  Global_node_mapping_df_query <- as.data.frame(cbind(Global_node_mapping_df_query,
                                                      metaNodesMap_Query_exp[paste0(Global_node_mapping_df_query$variable,
                                                                                    "_", Global_node_mapping_df_query$value),]))
  Global_node_mapping_df_2 <- as.data.frame(rbind(Global_node_mapping_df_ref,
                                                  Global_node_mapping_df_query))
  View(Global_node_mapping)
  sharedMarkers_cor <- c()
  for (i in sharedMarkers) {
    temp_ref_exp <- metaNodesMap_Ref_exp[paste0("ptRef_", Global_node_mapping$ptRef),i]
    temp_query_exp <- metaNodesMap_Query_exp[paste0("ptQuery_", Global_node_mapping$ptQuery),i]
    temp_cor_test <- cor.test(temp_ref_exp,
                              temp_query_exp)
    sharedMarkers_cor <- as.data.frame(rbind(sharedMarkers_cor,
                                             data.frame(Gene = i,
                                                        Cor = temp_cor_test$estimate,
                                                        P = temp_cor_test$p.value)))
  }
  
  sharedMarkers_cor$FDR <- p.adjust(sharedMarkers_cor$P,
                                    method = "fdr")
  ggplot(data = sharedMarkers_cor,
         aes(x = Cor)) +
    geom_histogram()
  
  selected_sharedMarkers_cor <- sharedMarkers_cor[abs(sharedMarkers_cor$Cor) > 0.8,]
  selected_sharedMarkers_cor_positive <- selected_sharedMarkers_cor[selected_sharedMarkers_cor$Cor > 0,]
  selected_sharedMarkers_cor_negative <- selected_sharedMarkers_cor[selected_sharedMarkers_cor$Cor < 0,]
  
  lines_cosine <- function(exp) {
    x <- 1:nrow(exp)
    y1 <- seq(from = 0, to = 1,
              length.out = length(x)) # 单调递增
    y2 <- seq(from = 1, to = 0,
              length.out = length(x)) # 单调递减
    y3 <- ((1/100)*x-1)^2 # 先减后增
    y4 <- 1-y3 # 先增后减
    cosine_df <- c()
    for (i in colnames(exp)) {
      cosine_1 <- lsa::cosine(cbind(exp[,i], y1))[1,2]
      cosine_2 <- lsa::cosine(cbind(exp[,i], y2))[1,2]
      cosine_3 <- lsa::cosine(cbind(exp[,i], y3))[1,2]
      cosine_4 <- lsa::cosine(cbind(exp[,i], y4))[1,2]
      
      cosine_df <- as.data.frame(rbind(cosine_df,
                                       data.frame(Gene = i,
                                                  y1 = cosine_1,
                                                  y2 = cosine_2,
                                                  y3 = cosine_3,
                                                  y4 = cosine_4)))
    }
    return(cosine_df)
  }
  selected_sharedMarkers_cor_positive <- selected_sharedMarkers_cor_positive[order(selected_sharedMarkers_cor_positive$Cor,
                                                                                   decreasing = T),]
  metaNodesMap_Ref_exp_positive <- metaNodesMap_Ref_exp[,selected_sharedMarkers_cor_positive$Gene]
  metaNodesMap_Query_exp_positive <- metaNodesMap_Query_exp[,selected_sharedMarkers_cor_positive$Gene]
  cosine_df <- lines_cosine(metaNodesMap_Ref_exp_positive)
  Type1 <- cosine_df$Gene[cosine_df$y1 >= 0.9 & cosine_df$y2 < 0.9]
  Type1 <- unlist(apply(metaNodesMap_Ref_exp_positive[,Type1,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type1 <- Type1[order(Type1, decreasing = F)]
  Type2 <- cosine_df$Gene[cosine_df$y1 < 0.9 & cosine_df$y2 >= 0.9]
  Type2 <- unlist(apply(metaNodesMap_Ref_exp_positive[,Type2,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type2 <- Type2[order(Type2, decreasing = F)]
  
  library(ComplexHeatmap)
  col_fun1 <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  col_anno_ref <- HeatmapAnnotation(
    Ref_Celltype = col_celltype,
    col = list(Ref_Celltype = RNA_Celltype_color)
  )
  col_anno_query <- HeatmapAnnotation(
    Query_Celltype = row_celltype,
    col = list(Query_Celltype = Limb_sub_color)
  )
  row_anno_gene <- rowAnnotation(
    Gene_type = c(rep("Type2", length(Type2)),
                  rep("Type1", length(Type1))),
    col = list(Gene_type = c("Type1" = "#F4A460",
                             "Type2" = "#1874CD"))
  )
  h1 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Ref_exp_positive[,c(names(Type2), names(Type1))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun1(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof RefTraj",
                                top_annotation = col_anno_ref,
                                left_annotation = row_anno_gene)
  col_fun2 <- circlize::colorRamp2(c(0, 0.5, 1), c("#009ACD", "#FFFACD", "#551A8B"))
  h2 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Query_exp_positive[,c(names(Type2), names(Type1))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun2(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof QueryTraj",
                                top_annotation = col_anno_query,
  )
  pdf("./轨迹比较/CellAlign_3.gene_heatmap_positive_Antler_Limb.pdf",
      width = 9, height = 6)
  h1 + h2
  dev.off()
  
  {
    ggplot() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.title = element_text(colour = "black", size = 10),
            legend.text = element_text(colour = "black", size = 10),
            axis.text = element_text(colour = "black", size = 10),
            axis.title = element_text(colour = "black", size = 10)) +
      geom_line(data = Global_node_mapping_df_2,
                aes(x = value, y = Type1_mean, group = rowIndex),
                color = "gray") +
      geom_point(data = Global_node_mapping_df_2,
                 aes(x = value, y = Type1_mean, color = variable)) +
      labs(x = "Pseudotime", y = "Mean interpolated expression of type1 gene set",
           color = "Trajectory") +
      scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                    "ptRef" = "#3A5FCD"))
  }
  
  
  Global_node_mapping_df_2$Type1_mean <- rowMeans(Global_node_mapping_df_2[,names(Type1)])
  Global_node_mapping_df_2$Type2_mean <- rowMeans(Global_node_mapping_df_2[,names(Type2)])
  
  pdf("./轨迹比较/CellAlign_4.positive_type2_mean_exp_Antler_Limb.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type2_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type2_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type2 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  pdf("./轨迹比较/CellAlign_4.positive_type1_mean_exp_Antler_Limb.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type1_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type1_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type1 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  
  metaNodesMap_Ref_exp_negative <- metaNodesMap_Ref_exp[,selected_sharedMarkers_cor_negative$Gene]
  metaNodesMap_Query_exp_negative <- metaNodesMap_Query_exp[,selected_sharedMarkers_cor_negative$Gene]
  cosine_df <- lines_cosine(metaNodesMap_Ref_exp_negative)
  Type3 <- cosine_df$Gene[cosine_df$y1 >= 0.9 & cosine_df$y2 < 0.9]
  Type3 <- unlist(apply(metaNodesMap_Ref_exp_negative[,Type3,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type3 <- Type3[order(Type3, decreasing = F)]
  Type4 <- cosine_df$Gene[cosine_df$y1 < 0.9 & cosine_df$y2 >= 0.9]
  Type4 <- unlist(apply(metaNodesMap_Ref_exp_negative[,Type4,
                                                      drop = F],
                        2, function(x){which.max(x)}))
  Type4 <- Type4[order(Type4, decreasing = F)]
  
  col_fun1 <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  col_anno_ref <- HeatmapAnnotation(
    Ref_Celltype = col_celltype,
    col = list(Ref_Celltype = RNA_Celltype_color)
  )
  col_anno_query <- HeatmapAnnotation(
    Query_Celltype = row_celltype,
    col = list(Query_Celltype = Limb_sub_color)
  )
  row_anno_gene <- rowAnnotation(
    Gene_type = c(rep("Type4", length(Type4)),
                  rep("Type3", length(Type3))),
    col = list(Gene_type = c("Type3" = "#836FFF",
                             "Type4" = "#D02090"))
  )
  h3 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Ref_exp_negative[,c(names(Type4), names(Type3))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun1(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof RefTraj",
                                top_annotation = col_anno_ref,
                                left_annotation = row_anno_gene)
  col_fun2 <- circlize::colorRamp2(c(0, 0.5, 1), c("#009ACD", "#FFFACD", "#551A8B"))
  h4 <- ComplexHeatmap::Heatmap(t(metaNodesMap_Query_exp_negative[,c(names(Type4), names(Type3))]),
                                cluster_rows = F, cluster_columns = F, col = col_fun2(seq(0, 1, length = 100)),
                                show_row_names = F, show_column_names = F,
                                name = "Interpolated expression\nof QueryTraj",
                                top_annotation = col_anno_query,
  )
  pdf("./轨迹比较/CellAlign_5.gene_heatmap_negative_Antler_Limb.pdf",
      width = 9, height = 5)
  h3 + h4
  dev.off()
  
  Global_node_mapping_df_2$Type3_mean <- rowMeans(Global_node_mapping_df_2[,names(Type3)])
  Global_node_mapping_df_2$Type4_mean <- rowMeans(Global_node_mapping_df_2[,names(Type4)])
  
  pdf("./轨迹比较/CellAlign_6.negative_type3_mean_exp_Antler_Limb.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type3_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type3_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type3 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  pdf("./轨迹比较/CellAlign_6.positive_type4_mean_exp_Antler_Limb.pdf",
      width = 6.2, height = 4.5)
  ggplot() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(colour = "black", size = 10),
          legend.text = element_text(colour = "black", size = 10),
          axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10)) +
    geom_line(data = Global_node_mapping_df_2,
              aes(x = value, y = Type4_mean, group = rowIndex),
              color = "gray") +
    geom_point(data = Global_node_mapping_df_2,
               aes(x = value, y = Type4_mean, color = variable)) +
    labs(x = "Pseudotime", y = "Mean interpolated expression of type4 gene set",
         color = "Trajectory") +
    scale_color_manual(values = c("ptQuery" = "#DA70D6",
                                  "ptRef" = "#3A5FCD"))
  dev.off()
  
  library(clusterProfiler)
  GeneList <- read.gmt("msigdb.v2023.2.Hs.symbols.gmt")
  GeneList <- GeneList[GeneList$gene %in% geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]],]
  Type1_enrich <- enricher(gene = names(Type1), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type1_enrich@result <- Type1_enrich@result[Type1_enrich@result$pvalue < 0.05,]
  Type1_enrich_result <- Type1_enrich@result
  Type1_enrich_result$ID <- unlist(lapply(Type1_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type1_enrich_result$Description <- unlist(lapply(Type1_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type1_enrich_result) <- Type1_enrich_result$ID
  Type1_enrich@result <- Type1_enrich_result
  openxlsx::write.xlsx(Type1_enrich@result,
                       "./轨迹比较/CellAlign_7.type1_genes_enrichment_Antler_Limb.xlsx")
  Type2_enrich <- enricher(gene = names(Type2), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type2_enrich@result <- Type2_enrich@result[Type2_enrich@result$pvalue < 0.05,]
  Type2_enrich_result <- Type2_enrich@result
  Type2_enrich_result$ID <- unlist(lapply(Type2_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type2_enrich_result$Description <- unlist(lapply(Type2_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type2_enrich_result) <- Type2_enrich_result$ID
  Type2_enrich@result <- Type2_enrich_result
  openxlsx::write.xlsx(Type2_enrich@result,
                       "./轨迹比较/CellAlign_7.type2_genes_enrichment_Antler_Limb.xlsx")
  Type3_enrich <- enricher(gene = names(Type3), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type3_enrich@result <- Type3_enrich@result[Type3_enrich@result$pvalue < 0.05,]
  Type3_enrich_result <- Type3_enrich@result
  Type3_enrich_result$ID <- unlist(lapply(Type3_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type3_enrich_result$Description <- unlist(lapply(Type3_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type3_enrich_result) <- Type3_enrich_result$ID
  Type3_enrich@result <- Type3_enrich_result
  openxlsx::write.xlsx(Type3_enrich@result,
                       "./轨迹比较/CellAlign_7.type3_genes_enrichment_Antler_Limb.xlsx")
  
  Type4_enrich <- enricher(gene = names(Type4), TERM2GENE = GeneList,
                           qvalueCutoff = 1, pvalueCutoff = 1,
                           universe = geneAnnotation@listData[["genes"]]@elementMetadata@listData[["symbol"]])
  Type4_enrich@result <- Type4_enrich@result[Type4_enrich@result$pvalue < 0.05,]
  Type4_enrich_result <- Type4_enrich@result
  Type4_enrich_result$ID <- unlist(lapply(Type4_enrich_result$ID,
                                          function(x){
                                            temp <- unlist(strsplit(x, "_", fixed = T))
                                            paste0(tolower(temp), collapse = " ")
                                          }))
  Type4_enrich_result$Description <- unlist(lapply(Type4_enrich_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_", fixed = T))
                                                     paste0(tolower(temp), collapse = " ")
                                                   }))
  rownames(Type4_enrich_result) <- Type4_enrich_result$ID
  Type4_enrich@result <- Type4_enrich_result
  openxlsx::write.xlsx(Type4_enrich@result,
                       "./轨迹比较/CellAlign_7.type4_genes_enrichment_Antler_Limb.xlsx")
  
  Gene_type_df <- data.frame(Type = c(rep("Type1", length(Type1)),
                                      rep("Type2", length(Type2)),
                                      rep("Type3", length(Type3)),
                                      rep("Type4", length(Type4))),
                             Gene = c(names(Type1), names(Type2),
                                      names(Type3), names(Type4)))
  Gene_type_df_list <- split.data.frame(Gene_type_df, f = list(Gene_type_df$Type))
  openxlsx::write.xlsx(Gene_type_df_list,
                       "./轨迹比较/CellAlign_8.type_genes_list_Antler_Limb.xlsx")
  
}

######### 废弃
#########
{
  #####
  # seurat 整合效果很差
  {
    
    Limb_sub <- NormalizeData(Limb_sub)
    Limb_sub <- FindVariableFeatures(Limb_sub)
    Limb_anchors <- FindTransferAnchors(reference = RNA_singlet_sub, query = Limb_sub,
                                        dims = 1:30, reference.reduction = "pca", k.filter = NA)
    Limb_predictions <- TransferData(anchorset = Limb_anchors,
                                     refdata = colnames(RNA_singlet_sub),
                                     dims = 1:30)
    Limb_sub$Limb_predictions <- Limb_predictions[colnames(Limb_sub), 1]
    Limb_sub$Limb_predictions_pseudotime <- pseudotime(RNA_monocle3)[Limb_sub$Limb_predictions]
    Limb_sub$Pseudotime <- pseudotime(Limb_monocle3)[colnames(Limb_sub)]
    pdf("./轨迹比较/Limb_Predicted_pseudotime.pdf")
    ggplot(data = Limb_sub@meta.data, aes(x = Limb_predictions_pseudotime,
                                          y = Pseudotime)) +
      geom_point() +
      geom_smooth() +
      # stat_cor() +
      labs(x = "Predicted Pseudotime", y = "Pseudotime")
    dev.off()
    
    Tumor_sub <- NormalizeData(Tumor_sub)
    Tumor_sub <- FindVariableFeatures(Tumor_sub)
    Tumor_anchors <- FindTransferAnchors(reference = RNA_singlet_sub, query = Tumor_sub,
                                         dims = 1:30, reference.reduction = "pca", k.filter = NA)
    Tumor_predictions <- TransferData(anchorset = Tumor_anchors,
                                      refdata = colnames(RNA_singlet_sub),
                                      dims = 1:30)
    Tumor_sub$Tumor_predictions <- Tumor_predictions[colnames(Tumor_sub), 1]
    Tumor_sub$Tumor_predictions_pseudotime <- pseudotime(RNA_monocle3)[Tumor_sub$Tumor_predictions]
    Tumor_sub$Pseudotime <- pseudotime(Tumor_monocle3)[colnames(Tumor_sub)]
    pdf("./轨迹比较/Tumor_Predicted_pseudotime.pdf")
    ggplot(data = Tumor_sub@meta.data, aes(x = Tumor_predictions_pseudotime,
                                           y = Pseudotime)) +
      geom_point() +
      geom_smooth() +
      # stat_cor() +
      labs(x = "Predicted Pseudotime", y = "Pseudotime")
    dev.off()
    
    Discs_sub <- NormalizeData(Discs_sub)
    Discs_sub <- FindVariableFeatures(Discs_sub)
    Discs_anchors <- FindTransferAnchors(reference = RNA_singlet_sub, query = Discs_sub,
                                         dims = 1:30, reference.reduction = "pca")
    Discs_predictions <- TransferData(anchorset = Discs_anchors,
                                      refdata = colnames(RNA_singlet_sub),
                                      dims = 1:30)
    
  }
  
  #####
  # 动态规划算法
  library(clusterProfiler)
  library(AUCell)
  C5 <- read.gmt("c5.all.v2023.2.Hs.symbols.gmt")
  C5 <- split.data.frame(C5, f = list(C5$term))
  C5 <- lapply(C5, function(x){
    x[,2]
  })
  C5 <- C5[unlist(lapply(C5, function(x){length(x)})) >= 5]
  C2 <- read.gmt("c2.all.v2023.2.Hs.symbols.gmt")
  C2 <- split.data.frame(C2, f = list(C2$term))
  C2 <- lapply(C2, function(x){
    x[,2]
  })
  C2 <- C2[unlist(lapply(C2, function(x){length(x)})) >= 5]
  pathway <- c(C5, C2)
  
}




