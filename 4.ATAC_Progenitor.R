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

pathway_enrich_for_gene_cluster <- function(gmt,
                                            gene_cluster = data.frame(cluster, genes)) {
  universe <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/wx_create_gene_malu.txt",
                       sep = "\t", header = FALSE)
  universe <- unique(universe$V6)
  pathway <- lapply(gmt, function(x){
    read.gmt(x)
  })
  pathway <- data.table::rbindlist(pathway)
  pathway <- as.data.frame(pathway)
  gene_cluster_list <- split.data.frame(gene_cluster, f = list(gene_cluster[,1]))
  gene_num <- unlist(lapply(gene_cluster_list, function(x){nrow(x)}))
  gene_cluster_list <- gene_cluster_list[gene_num > 5]
  gene_enrich_list <- c()
  for (genes in 1:length(gene_cluster_list)) {
    genes <- gene_cluster_list[[genes]]
    temp <- enricher(gene = genes[,2], pvalueCutoff = 1,
                     qvalueCutoff = 1, TERM2GENE = pathway,
                     universe = universe)
    if (is.null(temp)) {
      gene_enrich_list <- c(gene_enrich_list,
                            list(temp_result))
    } else {
      temp_result <- temp@result
      temp_result <- temp_result[temp_result$pvalue < 0.05,]
      if (nrow(temp_result) == 0) {
        gene_enrich_list <- c(gene_enrich_list,
                              list(temp_result))
      } else {
        group <- unlist(lapply(temp_result$Description, function(x){
          unlist(strsplit(x, "_"))[1]
        }))
        temp_result <- data.frame(Group = group,
                                  temp_result)
        temp_result$Description <- unlist(lapply(temp_result$Description,
                                                 function(x){
                                                   temp <- unlist(strsplit(x, "_"))
                                                   temp <- temp[-1]
                                                   temp <- paste0(temp, collapse = " ")
                                                   temp <- tolower(temp)
                                                   Hmisc::capitalize(temp)
                                                 }))
        temp_result$ID <- temp_result$Description
        temp_result <- temp_result[!duplicated(temp_result$Description),]
        rownames(temp_result) <- temp_result$Description
        gene_enrich_list <- c(gene_enrich_list,
                              list(temp_result))
      }
    }
    
  }
  names(gene_enrich_list) <- names(gene_cluster_list)
  return(gene_enrich_list)
}

ArchR_DimPlot <- function(ArchR = ATAC_Singlet, reduction_dimname = c("tSNE_1", "tSNE_2"),
                          color = NULL, pt.size = 0.1, reduction = "TSNE_withBatchEffect", repel = FALSE,
                          color_by = "Sample", subset_cells = NULL, label = FALSE, label_size = 10) {
  embedding <- getEmbedding(ArchR, embedding = reduction)
  colnames(embedding) <- reduction_dimname
  ColData <- as.data.frame(getCellColData(ArchR))
  embedding <- as.data.frame(cbind(embedding,
                                   ColData))
  if (!is.null(subset_cells)) {
    embedding <- embedding[subset_cells,]
  }
  p <- ggplot() +
    geom_point(data = embedding,
               aes(x = embedding[,reduction_dimname[1]], y = embedding[,reduction_dimname[2]],
                   color = embedding[,color_by]), size = pt.size) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10, colour = "black")) +
    labs(x = reduction_dimname[1], y = reduction_dimname[2], color = color_by) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  if (!is.null(color)) {
    p <- p + scale_color_manual(values = color)
  }
  if (label == TRUE) {
    label_embedding <- aggregate.data.frame(embedding[,c(reduction_dimname[1],
                                                         reduction_dimname[2])],
                                            by = list(embedding[,color_by]),
                                            FUN = median)
    if (repel == TRUE) {
      p <- p + ggrepel::geom_text_repel(data = label_embedding,
                                        aes(x = label_embedding[,reduction_dimname[1]], y = label_embedding[,reduction_dimname[2]],
                                            label = Group.1), size = label_size)
    } else {
      p <- p + geom_text(data = label_embedding,
                         aes(x = label_embedding[,reduction_dimname[1]], y = label_embedding[,reduction_dimname[2]],
                             label = Group.1), size = label_size)
    }
    
  }
  p
}

peakAnnoEnrichment_custom <- function(
    Peak_list = peak_cluster_reorder_list,
    ArchRProj = ATAC_Chondrocyte,
    peakAnnotation = "Motif",
    matches = NULL,
    background = "all"
){
  
  tstart <- Sys.time()
  
  if(is.null(matches)){
    matches <- getMatches(ArchRProj, peakAnnotation)
  }
  
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  mcols(r1) <- NULL
  
  r2 <- getPeakSet(ArchRProj)
  pr2 <- paste(seqnames(r2),start(r2),end(r2),sep="_")
  mcols(r2) <- NULL
  
  rownames(matches) <- pr1
  
  if(tolower(background) %in% c("backgroundpeaks", "bgdpeaks", "background", "bgd")){
    method <- "bgd"
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(ArchRProj))
  }else{
    method <- "all"
  }
  
  enrichList <- lapply(1:length(Peak_list), function(x){
    idx <- match(Peak_list[[x]], rownames(matches))
    if(method == "bgd"){
      .computeEnrichment(matches, idx, c(idx, as.vector(bgdPeaks[idx,])))
    }else{
      .computeEnrichment(matches, idx, seq_len(nrow(matches)))
    }
  }) %>% SimpleList
  names(enrichList) <- names(Peak_list)
  
  assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
    d <- lapply(seq_along(enrichList), function(y){
      enrichList[[y]][colnames(matches),x,drop=FALSE]
    }) %>% Reduce("cbind",.)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
  names(assays) <- colnames(enrichList[[1]])
  assays <- rev(assays)
  out <- SummarizedExperiment::SummarizedExperiment(assays=assays)
  
  out
  
}

getTrajectory_custom <- function(
    ArchRProj = NULL,
    name = "Trajectory",
    useMatrix = "GeneScoreMatrix",
    groupEvery = 1,
    log2Norm = TRUE,
    scaleTo = 10000,
    smoothWindow = 11,
    threads = getArchRThreads(),
    break_n = 100
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = groupEvery, name = "groupEvery", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric", "null"))
  .validInput(input = smoothWindow, name = "smoothWindow", valid = c("integer"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  
  trajectory <- getCellColData(ArchRProj, name)
  trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]
  breaks <- seq(from = 0, to = max(trajectory[,1]), length.out = break_n + 1)
  if(!all(is.numeric(trajectory[,1]))){
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if(!all(trajectory[,1] >= 0 & trajectory[,1] <= 100)){
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }
  
  groupList <- lapply(seq_along(breaks), function(x){
    if(x == 1){
      NULL
    }else{
      rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] & trajectory[,1] <= breaks[x])]
    }
  })[-1]
  names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])
  
  if (any(unlist(lapply(groupList, function(x){length(x)})) == 0)) {
    message("There are no cells in some group lists")
    groupList <- groupList[unlist(lapply(groupList, function(x){length(x)})) > 0]
  }
  
  featureDF <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  matrixClass <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class")))
  
  message("Creating Trajectory Group Matrix..")
  groupMat <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF,
    groupList = groupList, 
    threads = threads, 
    verbose = FALSE, 
    useMatrix = useMatrix
  )
  
  #Scale
  if(!is.null(scaleTo)){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
    }else{
      groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
    }
  }
  
  if(log2Norm){
    if(any(groupMat < 0)){
      message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
    }else{
      groupMat <- log2(groupMat + 1)
    }
  }
  
  if(!is.null(smoothWindow)){
    
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) .centerRollMean(x, k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))
    
    #Create SE
    seTrajectory <- SummarizedExperiment(
      assays = SimpleList(
        smoothMat = as.matrix(smoothGroupMat), 
        mat = as.matrix(groupMat)
      ), 
      rowData = featureDF
    )
    if("name" %in% colnames(featureDF)){
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
    }else{
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
    }
    
  }else{
    
    colnames(groupMat) <- paste0(colnames(groupMat))
    
    #Create SE
    seTrajectory <- SummarizedExperiment(
      assays = SimpleList(
        mat = as.matrix(groupMat)
      ), 
      rowData = featureDF
    )
    if("name" %in% colnames(featureDF)){
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$name)
    }else{
      rownames(seTrajectory) <- paste0(featureDF$seqnames, ":", featureDF$start, "_", featureDF$end)
    }
    
  }
  
  metadata(seTrajectory)$Params <- list(
    useMatrix = useMatrix, 
    matrixClass = matrixClass,
    scaleTo = scaleTo, 
    log2Norm = log2Norm, 
    smoothWindow = smoothWindow, 
    date = Sys.Date()
  )
  
  return(list(groupList = groupList,
              seTrajectory = seTrajectory))
  
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

ArchR_FeaturePlot <- function(ArchR_obj, Exp, impute_weight,
                              gene, reduction = "UMAP_Harmony",
                              reduction_dim_names = c("UMAP_1", "UMAP_2"),
                              max_value = 0.9, min_value = 0.9) {
  embedding <- getEmbedding(ArchR_obj, embedding = reduction)
  colnames(embedding) <- reduction_dim_names
  embedding <- as.data.frame(cbind(embedding,
                                   Exp[rownames(embedding),]))
  feature <- gene
  max_value <- quantile(embedding[,feature], max_value)
  min_value <- quantile(embedding[,feature], 1-min_value)
  # min_value <- 0
  embedding[embedding[,feature] > max_value, feature] <- max_value
  embedding[embedding[,feature] < min_value, feature] <- min_value
  embedding <- embedding[order(embedding[,feature], decreasing = FALSE), ]
  p <- ggplot(data = embedding,
              aes(x = embedding[,1], y = embedding[,2], color = embedding[,feature])) +
    ggrastr::rasterise(geom_point(size = 0.1), dpi = 300) +
    theme_bw() +
    ggtitle(feature) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10, colour = "black"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
    labs(x = "UMAP_1", y = "UMAP_2", color = "Gene activity") +
    scale_color_gradientn(colours = paletteContinuous(set = "solarExtra", n = 100))
  return(p)
  
}

setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
.libPaths("/home/heshidian/R/x86_64-pc-linux-gnu-library/4.4")


library(Seurat)
library(monocle)
library(monocle3)
library(harmony)
library(tidyverse)
library(patchwork)
library(BSgenome.malu)
library(ArchR)
library(ggalt)
library(Seurat)
library(parallel)
library(clustree)
library(dplyr)
library(patchwork)
library(BSgenome.malu)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)



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


##############
setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/")

RNA_singlet <- readRDS("RNA_1.RNA_singlet.rds")
RNA_singlet$percent.mt
Idents(RNA_singlet) <- RNA_singlet$Celltype_rename
DimPlot(RNA_singlet)
RNA_Celltype_color <- c("#A4C9DD", "#2572A9", "#ADD487", "#399938",
                        "#F19695", "#D5231E", "#F5BB6F", "#EF7C1C",
                        "#C6B0D2", "#653B90")
names(RNA_Celltype_color) <- levels(RNA_singlet$Celltype_rename)
RNA_group_color <- c("#E7211A", "#EFEA3C", "#72C8D5", "#6AB82D", "#18499E")
names(RNA_group_color) <- c("RM", "PC", "TZ", "CA", "MC")

RNA_singlet_sub <- subset(RNA_singlet, subset = Celltype_rename %in% c("AnMCs",
                                                                       # "Proliferative_progenitor cells",
                                                                       "Progenitor cells",
                                                                       "Chondrocytes",
                                                                       "Hypertrophic chondrocytes"))
DimPlot(RNA_singlet_sub, cols = RNA_Celltype_color)

RNA_singlet_sub$Celltype_rename
table(RNA_singlet_sub$Celltype_rename)


pathToMacs2 <- "/home/heshidian/mambaforge/envs/common/bin/macs2"

# TF Motif
{
  library(motifmatchr)
  library(TFBSTools)
  TF_info <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/Homo_sapiens_2024_11_19_4_17_am/TF_Information.txt",
                      sep = "\t")
  TF_info <- TF_info[-which(TF_info$TF_Status == "N"),]
  Bind_site_length <- c()
  Average_icscore <- c()
  for (i in 1:nrow(TF_info)) {
    if (file.exists(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/Homo_sapiens_2024_11_19_4_17_am/pwms_all_motifs/",
                           TF_info$Motif_ID[i], ".txt"))) {
      temp <- universalmotif::read_cisbp(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/Homo_sapiens_2024_11_19_4_17_am/pwms_all_motifs/",
                                                TF_info$Motif_ID[i], ".txt"))
      Bind_site_length <- c(Bind_site_length,
                            ncol(temp@motif))
      Average_icscore <- c(Average_icscore,
                           temp@icscore / ncol(temp@motif))
    } else {
      Bind_site_length <- c(Bind_site_length,
                            NA)
      Average_icscore <- c(Average_icscore,
                           NA)
    }
    
  }
  TF_info$Bind_site_length <- Bind_site_length
  TF_info$Average_icscore <- Average_icscore
  TF_info <- TF_info[!is.na(TF_info$Average_icscore),]
  TF_info_list <- split.data.frame(TF_info,
                                   f = list(TF_info$TF_Name))
  for (i in 1:length(TF_info_list)) {
    TF_info_list[[i]] <- TF_info_list[[i]][which.max(TF_info_list[[i]]$Average_icscore),]
  }
  TF_info <- as.data.frame(rbindlist(TF_info_list))
  Motif_list <- c()
  for (i in 1:nrow(TF_info)) {
    temp <- universalmotif::read_cisbp(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/Homo_sapiens_2024_11_19_4_17_am/pwms_all_motifs/",
                                              TF_info$Motif_ID[i], ".txt"))
    temp@name <- TF_info$TF_Name[i]
    Motif_list <- c(Motif_list,
                    list(temp))
  }
  Motif_pwm <- universalmotif::convert_motifs(Motif_list, class = "TFBSTools-PWMatrix")
  Motif_pwm <- do.call(PWMatrixList, Motif_pwm)
  names(Motif_pwm) <- unlist(lapply(Motif_pwm, function(x){
    x@name
  }))
  
}

ATAC_singlet <- readRDS("./ATAC/ATAC_4.pass_filteration_singlet_annotated.rds")
cells <- rownames(ATAC_singlet)[ATAC_singlet$Celltype %in% c("Proliferative_progenitor cells",
                                                             "Progenitor cells")]
ATAC_Progenitor <- subsetArchRProject(ArchRProj = ATAC_singlet,
                                      cells = cells,
                                      outputDirectory = "ATAC_Progenitor",
                                      dropCells = TRUE,
                                      force = TRUE)
# ATAC_Progenitor <- readRDS("ATAC_Progenitor.rds")
# ATAC_Progenitor <- addBgdPeaks(ATAC_Progenitor, force = TRUE)
# ATAC_Progenitor <- addMotifAnnotations(ArchRProj = ATAC_Progenitor,
#                                        motifPWMs = Motif_pwm,
#                                        name = "Motif",
#                                        cutOff = 1e-05,
#                                        force = TRUE)
# ATAC_Progenitor <- addDeviationsMatrix(
#   ArchRProj = ATAC_Progenitor,
#   peakAnnotation = "Motif",
#   force = TRUE
# )

markersGS <- getMarkerFeatures(
  ArchRProj = ATAC_Progenitor, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersMotif <- getMarkerFeatures(
  ArchRProj = ATAC_Progenitor, 
  useMatrix = "MotifMatrix", 
  groupBy = "Celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersExp <- getMarkerFeatures(
  ArchRProj = ATAC_Progenitor, 
  useMatrix = "GeneIntegrationMatrix_Co_harmony_Annotation_Group", 
  groupBy = "Celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerListGS <- getMarkers(markersGS,
                           cutOff = "FDR <= 1 & Log2FC >= -Inf")
markerListGS <- as.data.frame(markerListGS@listData$`Proliferative_progenitor cells`)
rownames(markerListGS) <- markerListGS$name

markerListMotif <- getMarkers_custom(markersMotif,
                                     cutOff = "FDR <= 1 & AUC >= -Inf")
markerListMotif <- as.data.frame(markerListMotif@listData$`Proliferative_progenitor cells`)
rownames(markerListMotif) <- markerListMotif$name

markerListExp <- getMarkers(markersExp,
                            cutOff = "FDR <= 1 & Log2FC >= -Inf")
markerListExp <- as.data.frame(markerListExp@listData$`Proliferative_progenitor cells`)
rownames(markerListExp) <- markerListExp$name

GS_Exp <- data.frame(Gene = rownames(markerListGS),
                     GS_Log2FC = markerListGS$Log2FC,
                     Exp_Log2FC = markerListExp[rownames(markerListGS), "Log2FC"])
ggplot(data = GS_Exp,
       aes(x = GS_Log2FC, y = Exp_Log2FC)) +
  geom_point()


markerListGS <- markerListGS[markerListGS$Pval < 0.05,]
rownames(markerListGS) <- markerListGS$name
markerListExp <- markerListExp[markerListExp$Pval < 0.05,]
rownames(markerListExp) <- markerListExp$name
markerListMotif <- markerListMotif[markerListMotif$Pval < 0.05,]
markerListMotif <- markerListMotif[order(markerListMotif$AUC,
                                         decreasing = TRUE),]
rownames(markerListMotif) <- markerListMotif$name

GS_Exp_Motif <- data.frame(Gene = rownames(markerListMotif),
                           Motif_AUC = markerListMotif$AUC,
                           Motif_Diff = markerListMotif$MeanDiff,
                           Exp_Log2FC = markerListExp[rownames(markerListMotif), "Log2FC"],
                           Exp_Diff = markerListExp[rownames(markerListMotif), "MeanDiff"],
                           GS_Log2FC = markerListGS[rownames(markerListMotif), "Log2FC"],
                           GS_Diff = markerListGS[rownames(markerListMotif), "MeanDiff"])
rownames(GS_Exp_Motif) <- GS_Exp_Motif$Gene

genes <- getFeatures(ATAC_Progenitor, useMatrix = "GeneScoreMatrix")
GOBP_DNA_REPAIR <- clusterProfiler::read.gmt("./ATAC_Progenitor/GOBP_DNA_REPAIR.v2024.1.Hs.gmt")
GOBP_DNA_REPAIR <- GOBP_DNA_REPAIR[GOBP_DNA_REPAIR$gene %in% genes,]
GOBP_DNA_DAMAGE_RESPONSE_P53 <- clusterProfiler::read.gmt("./ATAC_Progenitor/GOBP_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR.v2024.1.Hs.gmt")
GOBP_DNA_DAMAGE_RESPONSE_P53 <- GOBP_DNA_DAMAGE_RESPONSE_P53[GOBP_DNA_DAMAGE_RESPONSE_P53$gene %in% genes,]
GOBP_DNA_DAMAGE_RESPONSE <- clusterProfiler::read.gmt("./ATAC_Progenitor/GOBP_DNA_DAMAGE_RESPONSE.v2024.1.Hs.gmt")
GOBP_DNA_DAMAGE_RESPONSE <- GOBP_DNA_DAMAGE_RESPONSE[GOBP_DNA_DAMAGE_RESPONSE$gene %in% genes,]
GOBP_APOPTOTIC_SIGNALING <- clusterProfiler::read.gmt("./ATAC_Progenitor/GOBP_APOPTOTIC_SIGNALING_PATHWAY.v2024.1.Hs.gmt")
GOBP_APOPTOTIC_SIGNALING <- GOBP_APOPTOTIC_SIGNALING[GOBP_APOPTOTIC_SIGNALING$gene %in% genes,]

ATAC_Progenitor <- addModuleScore(
  ArchRProj = ATAC_Progenitor,
  useMatrix = "GeneScoreMatrix",
  name = "Module",
  features = list("DNA_REPAIR" = GOBP_DNA_REPAIR$gene,
                  "DNA_DAMAGE_RESPONSE_P53" = GOBP_DNA_DAMAGE_RESPONSE_P53$gene,
                  "DNA_DAMAGE_RESPONSE" = GOBP_DNA_DAMAGE_RESPONSE$gene,
                  "APOPTOTIC_SIGNALING" = GOBP_APOPTOTIC_SIGNALING$gene),
)
ATAC_Progenitor_meta <- getCellColData(ATAC_Progenitor)
ATAC_Progenitor_meta <- as.data.frame(ATAC_Progenitor_meta)
ATAC_Progenitor_meta <- ATAC_Progenitor_meta[,c("Celltype",
                                                "Module.DNA_REPAIR",
                                                "Module.DNA_DAMAGE_RESPONSE_P53",
                                                "Module.DNA_DAMAGE_RESPONSE",
                                                "Module.APOPTOTIC_SIGNALING")]
ATAC_Progenitor_meta$Celltype <- factor(ATAC_Progenitor_meta$Celltype,
                                        levels = c("Proliferative_progenitor cells",
                                                   "Progenitor cells"))
{
  p1 <- ggplot(data = ATAC_Progenitor_meta,
         aes(x = Celltype, y = Module.DNA_DAMAGE_RESPONSE, fill = Celltype)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    coord_cartesian(ylim = c(min(ATAC_Progenitor_meta$Module.DNA_DAMAGE_RESPONSE),
                             max(ATAC_Progenitor_meta$Module.DNA_DAMAGE_RESPONSE)*1.08)) +
    ggpubr::stat_compare_means(comparisons = list(c("Proliferative_progenitor cells",
                                                    "Progenitor cells"))) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Module score", title = "DNA damage response") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, colour = "black", hjust = 0.5, face = "bold"),
          legend.position = "none")
  p1
  p2 <- ggplot(data = ATAC_Progenitor_meta,
         aes(x = Celltype, y = Module.DNA_REPAIR, fill = Celltype)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    coord_cartesian(ylim = c(min(ATAC_Progenitor_meta$Module.DNA_REPAIR),
                             max(ATAC_Progenitor_meta$Module.DNA_REPAIR)*1.08)) +
    ggpubr::stat_compare_means(comparisons = list(c("Proliferative_progenitor cells",
                                                    "Progenitor cells"))) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Module score", title = "DNA repair") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, colour = "black", hjust = 0.5, face = "bold"),
          legend.position = "none")
  p2
  p3 <- ggplot(data = ATAC_Progenitor_meta,
               aes(x = Celltype, y = Module.DNA_DAMAGE_RESPONSE_P53, fill = Celltype)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    coord_cartesian(ylim = c(min(ATAC_Progenitor_meta$Module.DNA_DAMAGE_RESPONSE_P53),
                             max(ATAC_Progenitor_meta$Module.DNA_DAMAGE_RESPONSE_P53)*1.08)) +
    ggpubr::stat_compare_means(comparisons = list(c("Proliferative_progenitor cells",
                                                    "Progenitor cells"))) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Module score", title = "DNA damage response\nsignal transduction by P53") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, colour = "black", hjust = 0.5, face = "bold"),
          legend.position = "none")
  p3
  p4 <- ggplot(data = ATAC_Progenitor_meta,
               aes(x = Celltype, y = Module.APOPTOTIC_SIGNALING, fill = Celltype)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    ggpubr::stat_compare_means(comparisons = list(c("Proliferative_progenitor cells",
                                                    "Progenitor cells"))) +
    scale_fill_manual(values = RNA_Celltype_color) +
    coord_cartesian(ylim = c(min(ATAC_Progenitor_meta$Module.APOPTOTIC_SIGNALING),
                             max(ATAC_Progenitor_meta$Module.APOPTOTIC_SIGNALING)*1.08)) +
    theme_bw() +
    labs(x = "", y = "Module score", title = "Apoptotic signaling") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, colour = "black", hjust = 0.5, face = "bold"),
          legend.position = "none") +
    coord_cartesian(ylim = c(0, 1.1))
  p4
  pdf("./ATAC_Progenitor/ATAC_pathway_score_V2.pdf",
      width = 4.5, height = 9)
  cowplot::plot_grid(p1, p2, p3, p4,
                     nrow = 2,
                     align = "hv", axis = "tblr")
  dev.off()
}

{
  library(AUCell)
  RNA_progenitor <- subset(RNA_singlet, subset = Celltype_rename %in% c("Proliferative_progenitor cells",
                                                                         "Progenitor cells"))
  RNA_AUCell <- AUCell_run(exprMat = RNA_progenitor@assays$RNA@counts,
             geneSets = list("DNA_REPAIR" = GOBP_DNA_REPAIR$gene,
                             "DNA_DAMAGE_RESPONSE_P53" = GOBP_DNA_DAMAGE_RESPONSE_P53$gene,
                             "DNA_DAMAGE_RESPONSE" = GOBP_DNA_DAMAGE_RESPONSE$gene,
                             "APOPTOTIC_SIGNALING" = GOBP_APOPTOTIC_SIGNALING$gene))
  RNA_AUCell <- RNA_AUCell@assays@data@listData[["AUC"]]
  RNA_AUCell <- as.data.frame(t(RNA_AUCell))
  RNA_AUCell$Celltype <- RNA_progenitor@meta.data[rownames(RNA_AUCell), "Celltype_rename"]
  RNA_AUCell$Celltype <- factor(as.character(RNA_AUCell$Celltype),
                                levels = rev(c("Progenitor cells", "Proliferative_progenitor cells")))
  p1 <- ggplot(data = RNA_AUCell,
               aes(x = Celltype, y = DNA_DAMAGE_RESPONSE, fill = Celltype)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    coord_cartesian(ylim = c(min(RNA_AUCell$DNA_DAMAGE_RESPONSE),
                             max(RNA_AUCell$DNA_DAMAGE_RESPONSE)*1.08)) +
    ggpubr::stat_compare_means(comparisons = list(c("Proliferative_progenitor cells",
                                                    "Progenitor cells"))) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "AUCell score", title = "DNA damage response") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, colour = "black", hjust = 0.5, face = "bold"),
          legend.position = "none")
  p1
  p2 <- ggplot(data = RNA_AUCell,
               aes(x = Celltype, y = DNA_REPAIR, fill = Celltype)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    coord_cartesian(ylim = c(min(RNA_AUCell$DNA_REPAIR),
                             max(RNA_AUCell$DNA_REPAIR)*1.08)) +
    ggpubr::stat_compare_means(comparisons = list(c("Proliferative_progenitor cells",
                                                    "Progenitor cells"))) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "AUCell score", title = "DNA repair") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, colour = "black", hjust = 0.5, face = "bold"),
          legend.position = "none")
  p2
  p3 <- ggplot(data = RNA_AUCell,
               aes(x = Celltype, y = DNA_DAMAGE_RESPONSE_P53, fill = Celltype)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    coord_cartesian(ylim = c(min(RNA_AUCell$DNA_DAMAGE_RESPONSE_P53),
                             max(RNA_AUCell$DNA_DAMAGE_RESPONSE_P53)*1.08)) +
    ggpubr::stat_compare_means(comparisons = list(c("Proliferative_progenitor cells",
                                                    "Progenitor cells"))) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "AUCell score", title = "DNA damage response\nsignal transduction by P53") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, colour = "black", hjust = 0.5, face = "bold"),
          legend.position = "none")
  p3
  p4 <- ggplot(data = RNA_AUCell,
               aes(x = Celltype, y = APOPTOTIC_SIGNALING, fill = Celltype)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    ggpubr::stat_compare_means(comparisons = list(c("Proliferative_progenitor cells",
                                                    "Progenitor cells"))) +
    scale_fill_manual(values = RNA_Celltype_color) +
    coord_cartesian(ylim = c(min(RNA_AUCell$APOPTOTIC_SIGNALING),
                             max(RNA_AUCell$APOPTOTIC_SIGNALING)*1.08)) +
    theme_bw() +
    labs(x = "", y = "AUCell score", title = "Apoptotic signaling") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, colour = "black", hjust = 0.5, face = "bold"),
          legend.position = "none") +
    coord_cartesian(ylim = c(0, 0.12))
  p4
  pdf("./ATAC_Progenitor/RNA_pathway_score_V2.pdf",
      width = 4.5, height = 9)
  cowplot::plot_grid(p1, p2, p3, p4,
                     nrow = 2,
                     align = "hv", axis = "tblr")
  dev.off()
}

GOBP_APOPTOTIC_SIGNALING <- clusterProfiler::read.gmt("./ATAC_Progenitor/GOBP_APOPTOTIC_SIGNALING_PATHWAY.v2024.1.Hs.gmt")
GOBP_APOPTOTIC_SIGNALING <- GOBP_APOPTOTIC_SIGNALING[GOBP_APOPTOTIC_SIGNALING$gene %in% GS_Exp_Motif$Gene,]
GOBP_APOPTOTIC_SIGNALING <- GS_Exp_Motif[GOBP_APOPTOTIC_SIGNALING$gene, ]

GOBP_DNA_DAMAGE_RESPONSE <- clusterProfiler::read.gmt("./ATAC_Progenitor/GOBP_DNA_DAMAGE_RESPONSE.v2024.1.Hs.gmt")
GOBP_DNA_DAMAGE_RESPONSE <- GOBP_DNA_DAMAGE_RESPONSE[GOBP_DNA_DAMAGE_RESPONSE$gene %in% GS_Exp_Motif$Gene,]
GOBP_DNA_DAMAGE_RESPONSE <- GS_Exp_Motif[GOBP_DNA_DAMAGE_RESPONSE$gene, ]

GOBP_DNA_DAMAGE_RESPONSE_P53 <- clusterProfiler::read.gmt("./ATAC_Progenitor/GOBP_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR.v2024.1.Hs.gmt")
GOBP_DNA_DAMAGE_RESPONSE_P53 <- GOBP_DNA_DAMAGE_RESPONSE_P53[GOBP_DNA_DAMAGE_RESPONSE_P53$gene %in% GS_Exp_Motif$Gene,]
GOBP_DNA_DAMAGE_RESPONSE_P53 <- GS_Exp_Motif[GOBP_DNA_DAMAGE_RESPONSE_P53$gene, ]

GOBP_DNA_REPAIR <- clusterProfiler::read.gmt("./ATAC_Progenitor/GOBP_DNA_REPAIR.v2024.1.Hs.gmt")
GOBP_DNA_REPAIR <- GOBP_DNA_REPAIR[GOBP_DNA_REPAIR$gene %in% GS_Exp_Motif$Gene,]
GOBP_DNA_REPAIR <- GS_Exp_Motif[GOBP_DNA_REPAIR$gene, ]


GS_Exp_Motif <- GS_Exp_Motif[!is.na(GS_Exp_Motif$Exp_Log2FC),]
GS_Exp_Motif <- GS_Exp_Motif[!is.na(GS_Exp_Motif$Motif_AUC),]
GS_Exp_Motif <- GS_Exp_Motif[abs(GS_Exp_Motif$Exp_Diff) > 0.05,]
GS_Exp_Motif$Type <- ""
GS_Exp_Motif$Type[GS_Exp_Motif$Exp_Log2FC >= 1 & GS_Exp_Motif$Motif_AUC >= 0.6] <- "Prolif"
GS_Exp_Motif$Type[GS_Exp_Motif$Exp_Log2FC <= -1 & GS_Exp_Motif$Motif_AUC <= 0.4] <- "Non-Prolif"
table(GS_Exp_Motif$Type)

Other <- GS_Exp_Motif[GS_Exp_Motif$Type == "",]
Prolif <- GS_Exp_Motif[GS_Exp_Motif$Type == "Prolif",]
Prolif_top5 <- Prolif[order(Prolif$Exp_Log2FC, decreasing = TRUE),]
Prolif_top5 <- Prolif_top5[1:5,]
Non_Prolif <- GS_Exp_Motif[GS_Exp_Motif$Type == "Non-Prolif",]
Non_Prolif_top5 <- Non_Prolif[order(Non_Prolif$Exp_Log2FC, decreasing = FALSE),]
Non_Prolif_top5 <- Non_Prolif_top5[1:5,]


library(ggfx)
pdf("./ATAC_Progenitor/1.positive_TFs.pdf",
    width = 4, height = 4)
ggplot() +
  with_outer_glow(
    geom_point(data = GS_Exp_Motif,
               aes(x = Exp_Log2FC, y = Motif_AUC),
               size = 3, color = "#D3D3D3", alpha = 0.4),
    colour = "#D3D3D3", sigma = 10, expand = 6, alpha = 0.6
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size = 12, color = "black",
                                  face = "bold", hjust = 0.5)) +
  geom_point(data = Prolif,
             aes(x = Exp_Log2FC, y = Motif_AUC),
             color = "black", size = 2.5, shape = 21,
             fill = "#2572A9", stroke = 1) +
  geom_point(data = Non_Prolif,
             aes(x = Exp_Log2FC, y = Motif_AUC),
             color = "black", size = 2.5, shape = 21,
             fill = "#ADD487", stroke = 1) +
  labs(x = "log2FC of TF expression",
       y = "AUC of TF motif deviation",
       title = "Proliferative PCs vs. PCs") +
  ggrepel::geom_text_repel(data = Non_Prolif_top5, aes(x = Exp_Log2FC, y = Motif_AUC, 
                                                       label = Gene),
                           size = 3.3, fontface = "italic",
                           box.padding = unit(1, "lines"),
                           point.padding = unit(0, "lines"), 
                           min.segment.length = 0,
                           segment.color = "black",
                           colour="#000000",
                           show.legend = FALSE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 10)) +
  ggrepel::geom_text_repel(data = Prolif_top5, aes(x = Exp_Log2FC, y = Motif_AUC, 
                                                   label = Gene),
                           size = 3.3, fontface = "italic",
                           box.padding = unit(1.8, "lines"),
                           point.padding = unit(0, "lines"), 
                           min.segment.length = 0,
                           segment.color = "black",
                           colour="#000000",
                           show.legend = FALSE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 10)) +
  geom_vline(xintercept = 1, colour = "696969", linetype = "dashed") +
  geom_vline(xintercept = -1, colour = "696969", linetype = "dashed") +
  geom_hline(yintercept = 0.6, colour = "696969", linetype = "dashed") +
  geom_hline(yintercept = 0.4, colour = "696969", linetype = "dashed")
dev.off()


ATAC_Progenitor <- addIterativeLSI(ArchRProj = ATAC_Progenitor,
                                   useMatrix = "TileMatrix",
                                   name = "IterativeLSI",
                                   force = TRUE)
ATAC_Progenitor <- addHarmony(
  ArchRProj = ATAC_Progenitor,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = T
)
ATAC_Progenitor <- addUMAP(
  ArchRProj = ATAC_Progenitor, 
  reducedDims = "IterativeLSI", 
  name = "UMAP_withBatchEffect", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)
ATAC_Progenitor <- addUMAP(
  ArchRProj = ATAC_Progenitor, 
  reducedDims = "Harmony", 
  name = "UMAP_Harmony", 
  nNeighbors = 30, 
  minDist = 0.2, 
  metric = "cosine",
  force = T,
  n_neighbors = 20,
  local_connectivity = 1.1
)

pdf("./ATAC_Progenitor/2.Celltypes_UMAP.pdf",
    width = 6.5, height = 4)
ArchR_DimPlot(ArchR = ATAC_Progenitor,
              reduction_dimname = c("UMAP_1", "UMAP_2"),
              reduction = "UMAP_Harmony",
              color_by = "Celltype", pt.size = 0.1,
              label = TRUE, label_size = 3,
              color = RNA_Celltype_color, repel = FALSE) +
  ggtitle("ATAC_UMAP_Harmony") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  labs(color = "Celltype")
dev.off()


ATAC_Progenitor <- addImputeWeights(ATAC_Progenitor)
impute_weight <- getImputeWeights(ATAC_Progenitor)
ATAC_Progenitor_Exp <- getMatrixFromProject(ATAC_Progenitor,
                                            useMatrix = "GeneIntegrationMatrix_Co_harmony_Annotation_Group")
ATAC_Progenitor_Exp@assays@data@listData[["GeneIntegrationMatrix_Co_harmony_Annotation_Group"]]@Dimnames[[1]] <- ATAC_Progenitor_Exp@elementMetadata@listData[["name"]]
ATAC_Progenitor_Exp <- ATAC_Progenitor_Exp@assays@data@listData[["GeneIntegrationMatrix_Co_harmony_Annotation_Group"]]
ATAC_Progenitor_Exp <- as.matrix(ATAC_Progenitor_Exp[c(Prolif$Gene,
                                                       Non_Prolif$Gene),])
ATAC_Progenitor_Exp_2 <- imputeMatrix(ATAC_Progenitor_Exp,
                                      imputeWeights = impute_weight,
                                      threads = 20)
ATAC_Progenitor_Exp_2 <- as.data.frame(t(ATAC_Progenitor_Exp_2))
ATAC_Progenitor_Exp <- as.data.frame(t(ATAC_Progenitor_Exp))

ATAC_Progenitor_Motif <- getMatrixFromProject(ATAC_Progenitor,
                                              useMatrix = "MotifMatrix")
ATAC_Progenitor_Motif <- ATAC_Progenitor_Motif@assays@data@listData[["z"]]
ATAC_Progenitor_Motif_2 <- imputeMatrix(ATAC_Progenitor_Motif,
                                        imputeWeights = impute_weight,
                                        threads = 20)
ATAC_Progenitor_Motif_2 <- as.data.frame(t(ATAC_Progenitor_Motif_2))
ATAC_Progenitor_Motif <- data.frame(t(ATAC_Progenitor_Motif))

pdf("./ATAC_Progenitor/3.TF_Motif_Exp.pdf",
    width = 9, height = 3.5)
for (i in c(Prolif$Gene,
            Non_Prolif$Gene)) {
  p1 <- ArchR_FeaturePlot(ArchR_obj = ATAC_Progenitor,
                          Exp = ATAC_Progenitor_Exp_2,
                          gene = i, reduction = "UMAP_Harmony",
                          reduction_dim_names = c("UMAP_1", "UMAP_2"),
                          max_value = 0.95, min_value = 1) +
    scale_color_gradientn(colors = c("#000080", "#008B8B", "#DCDCDC", "#FA8072", "red")) +
    labs(title = paste0(i, " | Expression"),
         color = "Expression")
  p2 <- ArchR_FeaturePlot(ArchR_obj = ATAC_Progenitor,
                          Exp = ATAC_Progenitor_Motif_2,
                          gene = i, reduction = "UMAP_Harmony",
                          reduction_dim_names = c("UMAP_1", "UMAP_2"),
                          max_value = 0.95, min_value = 1) +
    scale_color_gradientn(colors = c("#4B0082", "#7B68EE", "#DCDCDC", "#F08080", "#FF4500")) +
    labs(title = paste0(i, " | Motif deviation"),
         color = "Motif deviation")
  p3 <- cowplot::plot_grid(p1, p2)
  print(p3)
}
dev.off()


Non_Prolif <- Non_Prolif[order(Non_Prolif$Exp_Log2FC,
                               decreasing = FALSE),]
Prolif <- Prolif[order(Prolif$Exp_Log2FC,
                       decreasing = TRUE),]
openxlsx::write.xlsx(rbind(Prolif,
                           Non_Prolif),
                     "./ATAC_Progenitor/1.positiive_TFs.xlsx")



#Get Reduced Dims
rD <- getReducedDims(ATAC_Progenitor, reducedDims = "IterativeLSI", corCutOff = 0.75, dimsToUse = 1:30)
#Subsample
set.seed(123)
idx <- sample(seq_len(nrow(rD)), 500, replace = !nrow(rD) >= 500)
#KNN Matrix
knnObj <- .computeKNN(data = rD, query = rD[idx,], k = 100)
#Determin Overlap
keepKnn <- determineOverlapCpp(knnObj, floor(0.8 * 100))
#Keep Above Cutoff
knnObj <- knnObj[keepKnn==0,]
#Convert To Names List
knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
  rownames(rD)[knnObj[x, ]]
}) %>% SimpleList

#Features
geneSet <- .getFeatureDF(getArrowFiles(ATAC_Progenitor),
                         "GeneIntegrationMatrix_Co_harmony_Annotation_Group", threads = 300)
geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
geneDF <- mcols(geneStart)
geneDF$seqnames <- seqnames(geneStart)
#Group Matrix RNA expression
groupMatExp <- .getGroupMatrix(
  ArrowFiles = getArrowFiles(ATAC_Progenitor), 
  featureDF = geneDF, 
  groupList = knnObj, 
  useMatrix = "GeneIntegrationMatrix_Co_harmony_Annotation_Group",
  threads = 30,
  verbose = FALSE
)
rownames(groupMatExp) <- geneDF$name
groupMatExp <- log2(groupMatExp + 1)
groupMatExp <- t(groupMatExp)

Motif_peak_match <- getPeakAnnotation(ATAC_Progenitor, name = "Motif")
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

markersPeak <- getMarkerFeatures(
  ArchRProj = ATAC_Progenitor, 
  useMatrix = "PeakMatrix", 
  groupBy = "Celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerListPeak <- getMarkers(markersPeak,
                             cutOff = "FDR <= 1 & Log2FC >= -Inf")
markerListPeak <- as.data.frame(markerListPeak@listData$`Proliferative_progenitor cells`)
rownames(markerListPeak) <- paste0(markerListPeak$seqnames, "_",
                                   markerListPeak$start, "_",
                                   markerListPeak$end)
markerListPeak_Prolif <- markerListPeak[markerListPeak$Log2FC >= 1 & markerListPeak$Pval < 0.05,]
markerListPeak_NonProlif <- markerListPeak[markerListPeak$Log2FC <= -1 & markerListPeak$Pval < 0.05,]

TF_target_Prolif <- c()
for (i in unique(Prolif$Gene)) {
  temp_match <- Motif_peak_match[,i,drop = FALSE]
  temp_match$Gene <- peakset_df[rownames(temp_match), "nearestGene"]
  temp_match$Position <- peakset_df[rownames(temp_match), "peakType"]
  temp_match <- temp_match[temp_match[,1] == TRUE,]
  temp_match <- temp_match[temp_match$Position %in% c("Distal", "Promoter"),]
  # temp_match <- temp_match[temp_match$Position %in% c("Promoter"),]
  temp_match <- temp_match[rownames(temp_match) %in% rownames(markerListPeak_Prolif),]
  temp_match <- temp_match[temp_match$Gene %in% colnames(groupMatExp),]
  R <- c()
  P <- c()
  for (j in 1:nrow(temp_match)) {
    temp_test <- cor.test(groupMatExp[,temp_match$Gene[j]],
                          groupMatExp[,i])
    p <- temp_test$p.value
    r <- temp_test$estimate
    R <- c(R, r)
    P <- c(P, p)
  }
  temp_match$R <- R
  temp_match$P <- P
  temp_match <- data.frame(TF = i,
                           Target = temp_match[,2],
                           temp_match[,3:ncol(temp_match)])
  TF_target_Prolif <- as.data.frame(rbind(TF_target_Prolif,
                                          temp_match))
}
quantile(TF_target_Prolif$R, na.rm = TRUE)
TF_target_Prolif <- as.data.frame(na.omit(TF_target_Prolif))
TF_target_Prolif <- TF_target_Prolif[TF_target_Prolif$P < 0.05,]
TF_target_Prolif <- TF_target_Prolif[TF_target_Prolif$R > 0.3,]
{
  pathway_enrich_for_gene_cluster <- function(gmt,
                                              gene_cluster = data.frame(cluster, genes)) {
    universe <- read.csv("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/reference/wx_create_gene_malu.txt",
                         sep = "\t", header = FALSE)
    universe <- unique(universe$V6)
    pathway <- lapply(gmt, function(x){
      read.gmt(x)
    })
    pathway <- data.table::rbindlist(pathway)
    pathway <- as.data.frame(pathway)
    gene_cluster_list <- split.data.frame(gene_cluster, f = list(gene_cluster[,1]))
    gene_num <- unlist(lapply(gene_cluster_list, function(x){nrow(x)}))
    gene_cluster_list <- gene_cluster_list[gene_num > 5]
    gene_enrich_list <- c()
    for (genes in 1:length(gene_cluster_list)) {
      genes <- gene_cluster_list[[genes]]
      temp <- enricher(gene = genes[,2], pvalueCutoff = 1,
                       qvalueCutoff = 1, TERM2GENE = pathway,
                       universe = universe)
      if (is.null(temp)) {
        gene_enrich_list <- c(gene_enrich_list,
                              list(temp_result))
      } else {
        temp_result <- temp@result
        temp_result <- temp_result[temp_result$pvalue < 1,]
        if (nrow(temp_result) == 0) {
          gene_enrich_list <- c(gene_enrich_list,
                                list(temp_result))
        } else {
          group <- unlist(lapply(temp_result$Description, function(x){
            unlist(strsplit(x, "_"))[1]
          }))
          temp_result <- data.frame(Group = group,
                                    temp_result)
          temp_result$Description <- unlist(lapply(temp_result$Description,
                                                   function(x){
                                                     temp <- unlist(strsplit(x, "_"))
                                                     temp <- temp[-1]
                                                     temp <- paste0(temp, collapse = " ")
                                                     temp <- tolower(temp)
                                                     Hmisc::capitalize(temp)
                                                   }))
          temp_result$ID <- temp_result$Description
          temp_result <- temp_result[!duplicated(temp_result$Description),]
          rownames(temp_result) <- temp_result$Description
          gene_enrich_list <- c(gene_enrich_list,
                                list(temp_result))
        }
      }
      
    }
    names(gene_enrich_list) <- names(gene_cluster_list)
    return(gene_enrich_list)
  }
  Prolif_enrich <- pathway_enrich_for_gene_cluster(gmt = c("./Gene2Gene/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt",
                                          "./Gene2Gene/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt",
                                          # "./Gene2Gene/c5.go.bp.v2024.1.Hs.symbols.gmt",
                                          "./Gene2Gene/h.all.v2024.1.Hs.symbols.gmt"),
                                  gene_cluster = data.frame(cluster = "Prolif",
                                                            genes = unique(c(TF_target_Prolif$TF,
                                                                             TF_target_Prolif$Target))))
  Prolif_enrich <- Prolif_enrich[["Prolif"]]
  row_index <- grep(pattern = "damage", x = Prolif_enrich$Description, ignore.case = TRUE)
  temp <- Prolif_enrich[row_index,]
  Prolif_enrich <- Prolif_enrich[Prolif_enrich$pvalue < 0.05,]
  Prolif_enrich <- Prolif_enrich[order(Prolif_enrich$pvalue,
                                       decreasing = FALSE),]
  openxlsx::write.xlsx(Prolif_enrich, "./ATAC_Progenitor/Prolif_enrich.xlsx")
}

table(TF_target_Prolif$TF)
unique(TF_target_Prolif$TF)
unique(Prolif$Gene)
openxlsx::write.xlsx(TF_target_Prolif, "./ATAC_Progenitor/4.TF_target_Prolif.xlsx")

TF_target_Non_Prolif <- c()
for (i in unique(Non_Prolif$Gene)) {
  temp_match <- Motif_peak_match[,i,drop = FALSE]
  temp_match$Gene <- peakset_df[rownames(temp_match), "nearestGene"]
  temp_match$Position <- peakset_df[rownames(temp_match), "peakType"]
  temp_match <- temp_match[temp_match[,1] == TRUE,]
  temp_match <- temp_match[temp_match$Position %in% c("Distal", "Promoter"),]
  # temp_match <- temp_match[temp_match$Position %in% c("Promoter"),]
  temp_match <- temp_match[rownames(temp_match) %in% rownames(markerListPeak_NonProlif),]
  temp_match <- temp_match[temp_match$Gene %in% colnames(groupMatExp),]
  R <- c()
  P <- c()
  for (j in 1:nrow(temp_match)) {
    temp_test <- cor.test(groupMatExp[,temp_match$Gene[j]],
                          groupMatExp[,i])
    p <- temp_test$p.value
    r <- temp_test$estimate
    R <- c(R, r)
    P <- c(P, p)
  }
  temp_match$R <- R
  temp_match$P <- P
  temp_match <- data.frame(TF = i,
                           Target = temp_match[,2],
                           temp_match[,3:ncol(temp_match)])
  TF_target_Non_Prolif <- as.data.frame(rbind(TF_target_Non_Prolif,
                                              temp_match))
}
quantile(TF_target_Non_Prolif$R, na.rm = TRUE)
TF_target_Non_Prolif <- as.data.frame(na.omit(TF_target_Non_Prolif))
TF_target_Non_Prolif <- TF_target_Non_Prolif[TF_target_Non_Prolif$P < 0.05,]
TF_target_Non_Prolif <- TF_target_Non_Prolif[TF_target_Non_Prolif$R > 0.3,]
table(TF_target_Non_Prolif$TF)
unique(TF_target_Non_Prolif$TF)
unique(Non_Prolif$Gene)
openxlsx::write.xlsx(TF_target_Non_Prolif, "./ATAC_Progenitor/4.TF_target_Non_Prolif.xlsx")

{
  NonProlif_enrich <- pathway_enrich_for_gene_cluster(gmt = c("./Gene2Gene/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt",
                                                           "./Gene2Gene/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt",
                                                           # "./Gene2Gene/c5.go.bp.v2024.1.Hs.symbols.gmt",
                                                           "./Gene2Gene/h.all.v2024.1.Hs.symbols.gmt"),
                                                   gene_cluster = data.frame(cluster = "NonProlif",
                                                                             genes = unique(c(TF_target_Non_Prolif$TF,
                                                                                              TF_target_Non_Prolif$Target))))
  NonProlif_enrich <- NonProlif_enrich[["NonProlif"]]
  row_index <- grep(pattern = "damage", x = NonProlif_enrich$Description, ignore.case = TRUE)
  temp <- NonProlif_enrich[row_index,]
  NonProlif_enrich <- NonProlif_enrich[NonProlif_enrich$pvalue < 0.05,]
  NonProlif_enrich <- NonProlif_enrich[order(NonProlif_enrich$pvalue,
                                       decreasing = FALSE),]
  openxlsx::write.xlsx(NonProlif_enrich, "./ATAC_Progenitor/NonProlif_enrich.xlsx")
}

#######################
#######################
#######################
motifPositions <- getPositions(ATAC_Progenitor)
ATAC_Progenitor <- addGroupCoverages(ArchRProj = ATAC_Progenitor,
                                     groupBy = "Celltype",
                                     force = TRUE)
seFoot <- getFootprints(
  ArchRProj = ATAC_Progenitor, 
  positions = motifPositions[markerListMotif$name[750]], 
  groupBy = "Celltype"
)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = ATAC_Progenitor, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)




