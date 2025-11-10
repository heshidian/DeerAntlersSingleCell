violin_linebox <- function(meta, group.by = "Day",
                           features = c("nFeature_RNA"),
                           group.color = group_color, linewidth = 2
) {
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
    geom_violin(data = meta, aes(x = meta[,group.by],
                                 y = meta[, features], fill = meta[, group.by]),
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
percent_bar <- function(meta, Ident, Group, fill_color = NULL, fill_label = "",
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
random_select_cells_v4 <- function(seurat_obj = Epi_tumor,
                                   group_by = "seurat_clusters",
                                   cluster_equally = TRUE,
                                   assay = "RNA",
                                   slot = "counts",
                                   cell_num = 10000,
                                   seed = 123) {
  if (!(assay %in% SeuratObject::Assays(seurat_obj))) {
    stop(paste0("The query assay ( ", assay, " ) was not existed!"))
  }
  # if (!(slot %in% names(seurat_obj@assays[[assay]]@layers))) {
  #   stop(paste0("The query slot ( ", slot, " ) was not existed in the assay ( ", assay, " )!"))
  # }
  groups <- unique(as.character(seurat_obj@meta.data[, group_by]))
  if (cluster_equally) {
    per_cluster_num <- ceiling(cell_num / length(groups))
    group_selected <- lapply(groups, function(x){
      cells <- colnames(seurat_obj)[as.character(seurat_obj@meta.data[, group_by]) == x]
      temp_num <- ifelse(per_cluster_num > length(cells),
                         length(cells), per_cluster_num)
      if (length(cells) <= 100) {
        message(paste0("Caution! the cell number of cluster ( ", x, " ) is ", length(cells), ", less than 100!"))
      }
      if (temp_num == length(cells)) {
        message(paste0("All cells of cluster ( ", x, " ) were selected."))
      }
      set.seed(seed)
      cluster_random <- sample(cells, temp_num, replace = FALSE)
      cluster_random
    })
  } else {
    cluster_percent <- table(as.character(seurat_obj@meta.data[, group_by]))
    cluster_percent <- as.data.frame(cluster_percent)
    cluster_percent$percent <- cluster_percent$Freq / sum(cluster_percent$Freq)
    cluster_percent$random_num <- ceiling(cluster_percent$percent * cell_num)
    group_selected <- lapply(seq_len(nrow(cluster_percent)), function(x){
      cells <- colnames(seurat_obj)[as.character(seurat_obj@meta.data[, group_by]) == as.character(cluster_percent$Var1[x])]
      temp_num <- ifelse(cluster_percent$random_num[x] > length(cells),
                         length(cells), cluster_percent$random_num[x])
      if (length(cells) <= 100) {
        message(paste0("Caution! the cell number of cluster ( ", x, " ) is ", length(cells), ", less than 100!"))
      }
      if (temp_num == length(cells)) {
        message(paste0("All cells of cluster ( ", x, " ) were selected."))
      }
      set.seed(seed)
      cluster_random <- sample(cells, temp_num, replace = FALSE)
      cluster_random
    })
  }
  
  selected_cells <- unlist(group_selected)
  result <- GetAssayData(seurat_obj,
                         assay = assay, slot = slot)
  result <- result[,selected_cells]
  result <- result[rowSums(result) > 0,]
  return(result)
}
getMarkers_custom <- function(
    seMarker = NULL,
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    n = NULL,
    returnGR = FALSE
){
  
  .validInput(input = seMarker, name = "seMarker", valid = c("SummarizedExperiment"))
  .validInput(input = cutOff, name = "cutOff", valid = c("character"))
  .validInput(input = n, name = "n", valid = c("integer", "null"))
  .validInput(input = returnGR, name = "returnGR", valid = c("boolean"))
  
  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }
  
  if(returnGR){
    
    if(metadata(seMarker)$Params$useMatrix != "PeakMatrix"){
      stop("Only markers can be returned as GRanges when PeakMatrix!")
    }
    
    rr <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start, rowData(seMarker)$end))
    
    grL <- lapply(seq_len(ncol(passMat)), function(x){
      idx <- which(passMat[, x])
      rrx <- rr[idx]
      rrx$AUC <- SummarizedExperiment::assays(seMarker[idx, ])[["AUC"]][, x]
      rrx$Log2FC <- SummarizedExperiment::assays(seMarker[idx, ])[["Log2FC"]][, x]
      rrx$FDR <- SummarizedExperiment::assays(seMarker[idx, ])[["FDR"]][, x]
      if("MeanDiff" %in% assayNames){
        rrx$MeanDiff <- SummarizedExperiment::assays(seMarker[idx, ])[["MeanDiff"]][, x]
      }
      rrx <- rrx[order(rrx$FDR),,drop=FALSE]
      if(!is.null(n)){
        if(n < nrow(rrx)){
          rrx <- rrx[seq_len(n), , drop = FALSE]
        }
      }
      rrx
    }) %>% SimpleList
    
    names(grL) <- colnames(seMarker)
    
    grL <- grL[gtools::mixedsort(names(grL))]
    
    return(grL)
    
  }else{
    
    markerList <- lapply(seq_len(ncol(passMat)), function(x){
      idx <- which(passMat[, x])
      rrx <- SummarizedExperiment::rowData(seMarker[idx,])
      rrx$AUC <- SummarizedExperiment::assays(seMarker[idx, ])[["AUC"]][, x]
      rrx$Log2FC <- SummarizedExperiment::assays(seMarker[idx, ])[["Log2FC"]][, x]
      rrx$FDR <- SummarizedExperiment::assays(seMarker[idx, ])[["FDR"]][, x]
      rrx$Pval <- SummarizedExperiment::assays(seMarker[idx, ])[["Pval"]][, x]
      if("MeanDiff" %in% assayNames){
        rrx$MeanDiff <- SummarizedExperiment::assays(seMarker[idx, ])[["MeanDiff"]][, x]
      }
      rrx <- rrx[order(rrx$Pval),,drop=FALSE]
      if(!is.null(n)){
        if(n < nrow(rrx)){
          rrx <- rrx[seq_len(n), , drop = FALSE]
        }
      }
      rrx
    }) %>% SimpleList
    
    names(markerList) <- colnames(seMarker)
    
    return(markerList)
    
  }
  
}

setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
.libPaths("/home/heshidian/R/x86_64-pc-linux-gnu-library/4.4")
library(Seurat)
library(monocle)
library(monocle3)
library(harmony)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(pheatmap)
library(scales)
library(egg)
library(org.Hs.eg.db)
library(CytoTRACE2)
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
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)

RNA_singlet <- readRDS("./RNA_1.RNA_singlet.rds")

RNA_group_color <- c("#E7211A", "#EFEA3C", "#72C8D5", "#6AB82D", "#18499E")
names(RNA_group_color) <- c("RM", "PC", "TZ", "CA", "MC")

RNA_Celltype_color <- c("#A4C9DD", "#2572A9", "#ADD487", "#399938",
                        "#F19695", "#D5231E", "#F5BB6F", "#EF7C1C",
                        "#C6B0D2", "#653B90")
names(RNA_Celltype_color) <- c("AnMCs", "Proliferative_progenitor cells", "Progenitor cells",
                               "Chondrocytes", "Hypertrophic chondrocytes", "Chondroclasts",
                               "Mural cells", "Endothelial cells", "Monocytes_Macrophages",
                               "Mast cells")

nFeature_RNA <- violin_linebox(RNA_singlet@meta.data,
                               features = "nFeature_RNA", linewidth = 1.5,
                               group.by = "Celltype_rename",
                               group.color = RNA_Celltype_color) +
  NoLegend()
pdf("RNA_QC_nFeature_RNA.pdf", width = 4.5, height = 4)
nFeature_RNA
dev.off()

nCount_RNA <- violin_linebox(RNA_singlet@meta.data,
                               features = "nCount_RNA", linewidth = 1.5,
                               group.by = "Celltype_rename",
                               group.color = RNA_Celltype_color) +
  NoLegend()
pdf("RNA_QC_nCount_RNA.pdf", width = 4.5, height = 4)
nCount_RNA
dev.off()

percent.mt <- violin_linebox(RNA_singlet@meta.data,
                             features = "percent.mt", linewidth = 1.5,
                             group.by = "Celltype_rename",
                             group.color = RNA_Celltype_color) +
  NoLegend()
pdf("RNA_QC_percent.mt.pdf", width = 4.5, height = 4)
percent.mt
dev.off()

RNA_singlet$Group <- factor(RNA_singlet$Group,
                            levels = c("RM", "PC", "TZ", "CA", "MC"))

nFeature_RNA <- violin_linebox(RNA_singlet@meta.data,
                               features = "nFeature_RNA", linewidth = 1.5,
                               group.by = "Group",
                               group.color = RNA_group_color) +
  NoLegend()
pdf("RNA_QC_nFeature_RNA_group.pdf", width = 3, height = 4)
nFeature_RNA
dev.off()

RNA_singlet$nCount_RNA_log <- log2(RNA_singlet$nCount_RNA)
nCount_RNA <- violin_linebox(RNA_singlet@meta.data,
                             features = "nCount_RNA", linewidth = 1.5,
                             group.by = "Group",
                             group.color = RNA_group_color) +
  NoLegend()
pdf("RNA_QC_nCount_RNA_group.pdf", width = 3, height = 4)
nCount_RNA
dev.off()

percent.mt <- violin_linebox(RNA_singlet@meta.data,
                             features = "percent.mt", linewidth = 1.5,
                             group.by = "Group",
                             group.color = RNA_group_color) +
  NoLegend()
pdf("RNA_QC_percent.mt_group.pdf", width = 3, height = 4)
percent.mt
dev.off()
RNA_num <- table(RNA_singlet$Celltype_rename)
RNA_num <- as.data.frame(RNA_num)
pdf("RNA_QC_Celltype_num.pdf", width = 4.8, height = 4.5)
ggplot() +
  geom_bar(data = RNA_num,
           aes(x = Var1, y = log10(Freq), fill = Var1),
           stat = "identity", width = 0.6) +
  geom_text(data = RNA_num,
            aes(x = Var1, y = log10(Freq), label = Freq),
            hjust = 0, vjust = 0, angle = 45,size = 3.3) +
  scale_fill_manual(values = RNA_Celltype_color) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = "") +
  NoLegend()
dev.off()

RNA_group_num <- table(RNA_singlet$Group)
RNA_group_num <- as.data.frame(RNA_group_num)
RNA_group_num$Var1 <- factor(RNA_group_num$Var1,
                        levels = c("RM", "PC", "TZ", "CA", "MC"))
pdf("RNA_QC_Group_num.pdf", width = 3, height = 4.5)
ggplot() +
  geom_bar(data = RNA_group_num,
           aes(x = Var1, y = log10(Freq), fill = Var1),
           stat = "identity", width = 0.6) +
  geom_text(data = RNA_group_num,
            aes(x = Var1, y = log10(Freq), label = Freq),
            hjust = 0, vjust = 0, angle = 45,size = 3.3) +
  scale_fill_manual(values = RNA_group_color) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = "") +
  NoLegend()
dev.off()







seurat_out <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/Stereo-seq/pythonProject/seurat_out_annotated.rds")
Celltype_color <- c("#A4C9DD", "#D5231E", "#399938", "#EF7C1C",
                    "#F19695", "#C6B0D2", "#F5BB6F", "#483D8B",
                    "#ADD487", "#2572A9", "#1E90FF")
names(Celltype_color) <- levels(seurat_out$Celltype)

seurat_out@meta.data$nCount_Spatial_log <- log2(seurat_out@meta.data$nCount_Spatial)
nCount_Spatial <- violin_linebox(seurat_out@meta.data,
                             features = "nCount_Spatial", linewidth = 1.5,
                             group.by = "Celltype",
                             group.color = Celltype_color) +
  NoLegend()
nCount_Spatial
pdf("Spatial_QC_nCount_Spatial.pdf", width = 5, height = 4.5)
nCount_Spatial
dev.off()

seurat_out@meta.data$nFeature_Spatial_log <- log2(seurat_out@meta.data$nFeature_Spatial)
nFeature_Spatial <- violin_linebox(seurat_out@meta.data,
                                 features = "nFeature_Spatial", linewidth = 1.5,
                                 group.by = "Celltype",
                                 group.color = Celltype_color) +
  NoLegend()
nFeature_Spatial

pdf("Spatial_QC_nFeature_Spatial.pdf", width = 5, height = 4.5)
nFeature_Spatial
dev.off()

spatial_num <- table(seurat_out$Celltype)
spatial_num <- as.data.frame(spatial_num)
pdf("Spatial_QC_Celltype_num.pdf", width = 5, height = 4.5)
ggplot() +
  geom_bar(data = spatial_num,
           aes(x = Var1, y = log10(Freq), fill = Var1),
           stat = "identity", width = 0.6) +
  geom_text(data = spatial_num,
            aes(x = Var1, y = log10(Freq), label = Freq),
            hjust = 0, vjust = 0, angle = 45,size = 3.3) +
  scale_fill_manual(values = Celltype_color) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 4.7)) +
  labs(x = "") +
   NoLegend()
dev.off()






library(ggsankey)
library(ArchR)
library(ggalt)
library(Seurat)
library(parallel)
library(clustree)
library(dplyr)
library(patchwork)
library(BSgenome.malu)
library(ggrastr)

rscripts <- list.files(path = "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/ArchR_R", pattern = ".R$", recursive = F, full.names = F)
for (i in rscripts) {
  source(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/ArchR_R/",i))
}
rscripts_fixed <- list.files(path = "/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/ArchR_R_bug_fixed", pattern = ".R$", recursive = F, full.names = F)
for (i in rscripts_fixed) {
  source(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/2.Rice_cell_atlas/Ref/ArchR_R_bug_fixed/",i))
}

geneAnnotation <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/ArchR_GeneAnnotation.rds")
genomeAnnotation <- readRDS("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/ArchR_GenomeAnnotation.rds")

addArchRThreads(threads = 30)

ATAC_singlet <- readRDS("./ATAC/ATAC_4.pass_filteration_singlet_annotated.rds")
table(ATAC_singlet$Celltype)

ATAC_num <- table(ATAC_singlet$Celltype)
ATAC_num <- as.data.frame(ATAC_num)
pdf("ATAC_QC_Celltype_num.pdf", width = 4.8, height = 4.5)
ggplot() +
  geom_bar(data = ATAC_num,
           aes(x = Var1, y = log10(Freq), fill = Var1),
           stat = "identity", width = 0.6) +
  geom_text(data = ATAC_num,
            aes(x = Var1, y = log10(Freq), label = Freq),
            hjust = 0, vjust = 0, angle = 45,size = 3.3) +
  scale_fill_manual(values = RNA_Celltype_color) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 5.7)) +
  labs(x = "") +
  NoLegend()
dev.off()


