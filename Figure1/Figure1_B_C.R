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
ArchR_FeaturePlot <- function(ArchR_obj, Exp, impute_weight,
                              gene, reduction = "UMAP_Harmony",
                              reduction_dim_names = c("UMAP_1", "UMAP_2")) {
  embedding <- getEmbedding(ArchR_obj, embedding = reduction)
  colnames(embedding) <- reduction_dim_names
  embedding <- as.data.frame(cbind(embedding,
                                   temp_mat[rownames(embedding),]))
  feature <- gene
  max_value <- quantile(embedding[,feature], 0.99)
  min_value <- 0
  embedding[embedding[,feature] > max_value, feature] <- max_value
  embedding[embedding[,feature] < min_value, feature] <- min_value
  embedding <- embedding[order(embedding[,feature], decreasing = FALSE), ]
  p <- ggplot(data = embedding,
              aes(x = embedding[,1], y = embedding[,2], color = embedding[,feature])) +
    rasterise(geom_point(size = 0.1), dpi = 300) +
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
DotPlot_2 <- function(Exp = ATAC_Singlet_GS_log,
                      group = as.character(ATAC_meta[colnames(ATAC_Singlet_GS_log),"Annotation"]),
                      genes = c("Ptprc", "B2m", "Col1a1", "Dcn", "Vim",
                                "Cntnap5b", "Ntm", "Ctss", "C1qb",
                                "Apoe", "Gja1", "Tph1", "Ddc", "Asmt"),
                      genes_order = c("Ptprc", "B2m", "Col1a1", "Dcn", "Vim",
                                      "Cntnap5b", "Ntm", "Ctss", "C1qb",
                                      "Apoe", "Gja1", "Tph1", "Ddc", "Asmt"),
                      group_order = c("Immunocyte", "VLMC", "Neuron",
                                      "Microglia", "Astrocyte", "Pinealocyte"),
                      fill_label = "Gene Activity", size_label = "Percent",
                      high_color = "red", low_color = "blue") {
  genes <- genes[genes %in% rownames(Exp)]
  genes_order <- genes_order[genes_order %in% genes]
  Exp <- Exp[genes,, drop = FALSE]
  Exp_binary <- Exp
  Exp_binary[Exp_binary > 0] <- 1
  group_num <- table(group)
  Percent <- aggregate.data.frame(t(Exp_binary), by = list(group), FUN = sum)
  Percent <- reshape::melt.data.frame(Percent, id.vars = "Group.1")
  Percent$Percent <- Percent$value / group_num[Percent$Group.1]
  Percent <- Percent[Percent$variable %in% genes,]
  Percent$Percent <- round(Percent$Percent * 100)
  Mean_exp <- aggregate.data.frame(t(Exp), by = list(group), FUN = mean)
  for (i in 2:ncol(Mean_exp)) {
    Mean_exp[,i] <- (Mean_exp[,i] - mean(Mean_exp[,i])) / sd(Mean_exp[,i])
  }
  Mean_exp <- reshape::melt.data.frame(Mean_exp, id.vars = "Group.1")
  Mean_exp <- Mean_exp[Mean_exp$variable %in% genes,]
  df <- data.frame(group = as.character(Mean_exp$Group.1),
                   gene = as.character(Mean_exp$variable),
                   percent = as.numeric(Percent$Percent),
                   exp = Mean_exp$value)
  df$exp[df$exp > 1.5] <- 1.5
  df$exp[df$exp < -1.5] <- -1.5
  df$gene <- factor(df$gene, levels = rev(genes_order))
  df$group <- factor(df$group, levels = group_order)
  p <- ggplot(data = df, aes(x = gene, y = group, color = exp, size = percent)) +
    geom_point() +
    scale_size_continuous(range = c(0.01, 8)) +
    scale_color_gradient2(low = low_color, mid = "white", high = high_color,
                          midpoint = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text = element_text(size = 10, color = "black"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black")) +
    labs(x = "", y = "", color = fill_label, size = size_label)
  p
}
Two_layer_sankey_plot <- function(layer1 = ATAC_singlet$Cluster_Annotation,
                                  layer2 = ATAC_singlet$predictedGroup_Co_harmony_Annotation,
                                  node_color = RNA_Celltype_color,
                                  layer1_name = "ATAC",
                                  laryer2_name = "RNA",
                                  plot_title = "Sankey Diagram for snATAC-seq to snRNA-seq Mapping"){
  df_sankey_1 <- data.frame(x = layer1_name,
                            node = layer1,
                            next_x = laryer2_name,
                            next_node = layer2)
  df_sankey_2 <- data.frame(x = laryer2_name,
                            node = layer2,
                            next_x = "End",
                            next_node = "End")
  df_sankey <- as.data.frame(rbind(df_sankey_1,
                                   df_sankey_2))
  df_sankey <- as_tibble(df_sankey)
  df_sankey <- df_sankey %>%
    mutate(x = as.factor(x)) %>%
    mutate(next_x = as.factor(next_x))
  
  Consistency <- sum(layer1 == layer2)
  Consistency <- Consistency / length(layer1)
  
  p <- ggplot(df_sankey, aes(x = x, next_x = next_x,
                             node = node, next_node = next_node,
                             fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = "black", width = 0.3) +
    geom_sankey_label(size = 3.3, color = "black", fill = NA, label.size = 0) +
    scale_y_continuous(expand = c(0.02, 0.02)) +
    theme_sankey(base_size = 10) +
    labs(title = paste0(plot_title,
                        "\n", "Consistency: ", round(Consistency*100, 2), "%")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 13, colour = "black", face = "bold"),
          axis.title = element_blank()) +
    scale_fill_manual(values = node_color)
  return(p)
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

setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/ATAC")

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

inputFiles <- list.files(path = "./fragments", pattern = "fragments.tsv.gz$",
                         recursive = T, full.names = T)
names(inputFiles) <- unlist(lapply(inputFiles, function(x){
  unlist(strsplit(x, ".", fixed = T))[3]
}))

arrowFiles <- createArrowFiles(inputFiles = inputFiles,
                               sampleNames = names(inputFiles),
                               outputNames = names(inputFiles),
                               minTSS = 0, minFrags = 0,
                               excludeChr = c("chrMT"),
                               addTileMat = F,
                               addGeneScoreMat = F,
                               subThreading = F,
                               geneAnnotation = geneAnnotation,
                               genomeAnnotation = genomeAnnotation,
                               force = TRUE)

# arrowFiles <- list.files(path = "./ATAC", pattern = ".arrow",
#                          recursive = T, full.names = T)

proj_all <- ArchRProject(
  ArrowFiles = arrowFiles, 
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  outputDirectory = "ArchR_ATAC",
  copyArrows = T #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
quantile(proj_all$TSSEnrichment)
quantile(proj_all$nFrags)

ChromSizes <- getChromSizes(proj_all)
names(ChromSizes) <- ChromSizes@seqnames
proj_all <- addTileMatrix(proj_all,
                          chromSizes = getChromSizes(proj_all),
                          blacklist = NULL,
                          excludeChr = c("chrMT"),
                          force = TRUE)
proj_all <- addGeneScoreMatrix(input = proj_all,
                               matrixName = "GeneScoreMatrix",
                               force = T,
                               excludeChr = c("chrMT"))
proj_all$NucleosomeSignal <- proj_all$nDiFrags / proj_all$nMonoFrags
quantile(proj_all$NucleosomeSignal)

samples <- unique(proj_all$Sample)

# TSS&Fragament QC Plots
{
  filterFrags <- 0
  filterTSS <- 0
  Metadata_all <- getCellColData(proj_all)
  Metadata_all <- as.data.frame(Metadata_all)
  for (i in 1:length(samples)) {
    print(i)
    QCDir <- paste0("./ArchR_ATAC/1.QCplot/")
    Metadata <- Metadata_all[Metadata_all$Sample == samples[i],]
    ggtitle <- sprintf("%s\n%s\n%s",
                       paste0(unlist(strsplit(samples[i],".arrow",fixed=T))[1], "\nnCells Pass Filter = ", nrow(Metadata)),
                       paste0("Median Frags = ", median(Metadata$nFrags)),
                       paste0("Median TSS Enrichment = ", median(Metadata$TSSEnrichment))
    )
    gg <- ggPoint(
      x = pmin(log10(Metadata$nFrags), 5) + rnorm(length(Metadata$nFrags), sd = 0.00001),
      y = Metadata$TSSEnrichment + rnorm(length(Metadata$nFrags), sd = 0.00001), 
      colorDensity = TRUE,
      xlim = c(2.5, 5),
      ylim = c(0, max(Metadata$TSSEnrichment) * 1.05),
      baseSize = 6,
      continuousSet = "sambaNight",
      xlabel = "Log 10 (Unique Fragments)",
      ylabel = "TSS Enrichment",
      title = ggtitle,
      rastr = TRUE) + 
      geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.25) +
      geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.25)
    pdf(file.path(QCDir,paste0(samples[i],"-TSS_by_Unique_Frags.pdf")),
        width=6,height=6,onefile=FALSE)
    .fixPlotSize(gg, plotWidth = 6, plotHeight = 6)
    dev.off()
    
    
    fragSummary <- .fastFragmentInfo(
      ArrowFile = getArrowFiles(proj_all)[samples[i]], 
      cellNames = rownames(Metadata_all)[Metadata_all$Sample == samples[i]], 
      nucLength = 147
    )
    Metadata <- fragSummary[[1]]
    plotDF <- data.frame(
      x = seq_along(fragSummary[[2]]), 
      percent = 100 * fragSummary[[2]]/sum(fragSummary[[2]])
    )
    gg <- ggplot(plotDF, aes(x = x, y = percent)) + theme_ArchR(baseSize = 7) + 
      geom_line(col = "dodgerblue4", size = 0.5) + 
      coord_cartesian(xlim = c(0,750), ylim = c(0, max(plotDF$percent) * 1.1), expand = FALSE) + 
      xlab("Size of Fragments (bp) \n") + 
      ylab("Fragments (%)") + 
      ggtitle(paste0(samples[i],"\nnFrags = ", round(sum(Metadata[,2])/10^6, 2)," M\nFragment Size Distribution"))
    pdf(file.path(QCDir,paste0(samples[i],"-Fragment_Size_Distribution.pdf")),width=5,height=5,onefile=FALSE)
    .fixPlotSize(gg, plotWidth = 6, plotHeight = 6)
    dev.off()
    
    plotDF <- plotTSSEnrichment(ArchRProj = proj_all,
                                returnDF = TRUE)
    plotDF <- data.frame(x=plotDF$x,v=plotDF$smoothValue,group=plotDF$group)
    plotDF <- plotDF[plotDF$group == samples[i],]
    p <- ggplot(plotDF, aes(x,v)) +
      geom_line(size = 1, color = "dodgerblue4") +
      theme_ArchR() +
      xlab("Distance From Center (bp)") +
      ylab("Normalized Insertion Profile") +
      # scale_color_manual(values=pal) +
      scale_y_continuous(limits = c(0, max(plotDF$v)*1.05), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(plotDF$x), max(plotDF$x)), expand = c(0,0))
    pdf(file.path(QCDir,paste0(samples[i],"-TSS_Enrichment.pdf")),
        width = 4.5, height = 4.5)
    .fixPlotSize(p, plotWidth = 6, plotHeight = 6)
    dev.off()
  }
}

proj_all$NucleosomeSignal <- proj_all$nDiFrags / proj_all$nMonoFrags
sum(is.nan(proj_all$NucleosomeSignal))
quantile(proj_all$NucleosomeSignal)

samples <- unique(proj_all$Sample)

proj_all$TSSGroup_1 <- "TSS Enrichment >= 1"
proj_all$TSSGroup_1[proj_all$TSSEnrichment < 1] <- "TSS Enrichment < 1"
proj_all$TSSGroup_2 <- "TSS Enrichment >= 2"
proj_all$TSSGroup_2[proj_all$TSSEnrichment < 2] <- "TSS Enrichment < 2"
proj_all$TSSGroup_3 <- "TSS Enrichment >= 3"
proj_all$TSSGroup_3[proj_all$TSSEnrichment < 3] <- "TSS Enrichment < 3"
proj_all$TSSGroup_4 <- "TSS Enrichment >= 4"
proj_all$TSSGroup_4[proj_all$TSSEnrichment < 4] <- "TSS Enrichment < 4"
proj_all$TSSGroup_5 <- "TSS Enrichment >= 5"
proj_all$TSSGroup_5[proj_all$TSSEnrichment < 5] <- "TSS Enrichment < 5"

proj_all$Sample_TSSGroup_1 <- paste0(proj_all$Sample,
                                     " | ", proj_all$TSSGroup_1)
temp <- table(proj_all$Sample_TSSGroup_1)
temp2 <- paste0(names(temp), "\n(n = ", as.numeric(temp), ")")
names(temp2) <- names(temp)
proj_all$Sample_TSSGroup_1 <- temp2[proj_all$Sample_TSSGroup_1]
proj_all$Sample_TSSGroup_2 <- paste0(proj_all$Sample,
                                     " | ", proj_all$TSSGroup_2)
temp <- table(proj_all$Sample_TSSGroup_2)
temp2 <- paste0(names(temp), "\n(n = ", as.numeric(temp), ")")
names(temp2) <- names(temp)
proj_all$Sample_TSSGroup_2 <- temp2[proj_all$Sample_TSSGroup_2]
proj_all$Sample_TSSGroup_3 <- paste0(proj_all$Sample,
                                     " | ", proj_all$TSSGroup_3)
temp <- table(proj_all$Sample_TSSGroup_3)
temp2 <- paste0(names(temp), "\n(n = ", as.numeric(temp), ")")
names(temp2) <- names(temp)
proj_all$Sample_TSSGroup_3 <- temp2[proj_all$Sample_TSSGroup_3]
proj_all$Sample_TSSGroup_4 <- paste0(proj_all$Sample,
                                     " | ", proj_all$TSSGroup_4)
temp <- table(proj_all$Sample_TSSGroup_4)
temp2 <- paste0(names(temp), "\n(n = ", as.numeric(temp), ")")
names(temp2) <- names(temp)
proj_all$Sample_TSSGroup_4 <- temp2[proj_all$Sample_TSSGroup_4]
proj_all$Sample_TSSGroup_5 <- paste0(proj_all$Sample,
                                     " | ", proj_all$TSSGroup_5)
temp <- table(proj_all$Sample_TSSGroup_5)
temp2 <- paste0(names(temp), "\n(n = ", as.numeric(temp), ")")
names(temp2) <- names(temp)
proj_all$Sample_TSSGroup_5 <- temp2[proj_all$Sample_TSSGroup_5]

QCDir <- "./ArchR_ATAC/1.QCplot/TSS_enrich/"
for (g in c("Sample_TSSGroup_1", "Sample_TSSGroup_2", "Sample_TSSGroup_3",
            "Sample_TSSGroup_4", "Sample_TSSGroup_5")) {
  plotDF <- plotTSSEnrichment(ArchRProj = proj_all,
                              groupBy = g,
                              returnDF = TRUE)
  plotDF <- data.frame(x = plotDF$x, v = plotDF$smoothValue,
                       group = plotDF$group)
  
  p <- ggplot(plotDF, aes(x, v, color = group)) +
    geom_line(size = 1) +
    theme_ArchR() +
    xlab("Distance From Center (bp)") +
    ylab("Normalized Insertion Profile") +
    # scale_color_manual(values=pal) +
    scale_y_continuous(limits = c(0, max(plotDF$v)*1.05), expand = c(0,0)) +
    scale_x_continuous(limits = c(min(plotDF$x), max(plotDF$x)), expand = c(0,0)) +
    facet_wrap(".~group", ncol = 2)
  pdf(file.path(QCDir,paste0(g,"-TSS_Enrichment.pdf")),
      width = 12, height = 45)
  .fixPlotSize(p, plotWidth = 12, plotHeight = 40)
  dev.off()
}


proj_all$PassQC <- 0
quantile(proj_all$nFrags)
quantile(proj_all$TSSEnrichment)
quantile(proj_all$NucleosomeSignal)
proj_all$PassQC[proj_all$nFrags > 1000 & proj_all$TSSEnrichment > 1 & proj_all$NucleosomeSignal < 2] <- 1
table(proj_all$PassQC) # 去除18,497个细胞
Cells <- getCellNames(proj_all)
Cells <- Cells[proj_all$PassQC == 1]
saveRDS(proj_all, "ATAC_1.proj_all.rds")

# TSS&Fragament QC Plots
{
  filterFrags <- 1000
  filterTSS <- 1
  Metadata_all <- getCellColData(proj_all)
  Metadata_all <- as.data.frame(Metadata_all)
  Metadata_all <- Metadata_all[Cells,]
  for (i in 1:length(samples)) {
    print(i)
    QCDir <- paste0("./ArchR_ATAC/1.QCplot/")
    Metadata <- Metadata_all[Metadata_all$Sample == samples[i],]
    ggtitle <- sprintf("%s\n%s\n%s",
                       paste0(unlist(strsplit(samples[i],".arrow",fixed=T))[1], "\nnCells Pass Filter = ", nrow(Metadata)),
                       paste0("Median Frags = ", median(Metadata$nFrags)),
                       paste0("Median TSS Enrichment = ", median(Metadata$TSSEnrichment))
    )
    gg <- ggPoint(
      x = pmin(log10(Metadata$nFrags), 5) + rnorm(length(Metadata$nFrags), sd = 0.00001),
      y = Metadata$TSSEnrichment + rnorm(length(Metadata$nFrags), sd = 0.00001), 
      colorDensity = TRUE,
      xlim = c(2.5, 5),
      ylim = c(0, max(Metadata$TSSEnrichment) * 1.05),
      baseSize = 6,
      continuousSet = "sambaNight",
      xlabel = "Log 10 (Unique Fragments)",
      ylabel = "TSS Enrichment",
      title = ggtitle,
      rastr = TRUE) + 
      geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.25) +
      geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.25)
    pdf(file.path(QCDir,paste0(samples[i],"-TSS_by_Unique_Frags(filtered).pdf")),
        width=6,height=6,onefile=FALSE)
    .fixPlotSize(gg, plotWidth = 6, plotHeight = 6)
    dev.off()
    
    
    fragSummary <- .fastFragmentInfo(
      ArrowFile = getArrowFiles(proj_all)[samples[i]], 
      cellNames = rownames(Metadata_all)[Metadata_all$Sample == samples[i]], 
      nucLength = 147
    )
    Metadata <- fragSummary[[1]]
    plotDF <- data.frame(
      x = seq_along(fragSummary[[2]]), 
      percent = 100 * fragSummary[[2]]/sum(fragSummary[[2]])
    )
    gg <- ggplot(plotDF, aes(x = x, y = percent)) + theme_ArchR(baseSize = 7) + 
      geom_line(col = "dodgerblue4", size = 0.5) + 
      coord_cartesian(xlim = c(0,750), ylim = c(0, max(plotDF$percent) * 1.1), expand = FALSE) + 
      xlab("Size of Fragments (bp) \n") + 
      ylab("Fragments (%)") + 
      ggtitle(paste0(samples[i],"\nnFrags = ", round(sum(Metadata[,2])/10^6, 2)," M\nFragment Size Distribution"))
    pdf(file.path(QCDir,paste0(samples[i],"-Fragment_Size_Distribution(filtered).pdf")),width=5,height=5,onefile=FALSE)
    .fixPlotSize(gg, plotWidth = 6, plotHeight = 6)
    dev.off()
    
    plotDF <- plotTSSEnrichment(ArchRProj = proj_all,
                                returnDF = TRUE)
    plotDF <- data.frame(x=plotDF$x,v=plotDF$smoothValue,group=plotDF$group)
    plotDF <- plotDF[plotDF$group == samples[i],]
    p <- ggplot(plotDF, aes(x,v)) +
      geom_line(size = 1, color = "dodgerblue4") +
      theme_ArchR() +
      xlab("Distance From Center (bp)") +
      ylab("Normalized Insertion Profile") +
      # scale_color_manual(values=pal) +
      scale_y_continuous(limits = c(0, max(plotDF$v)*1.05), expand = c(0,0)) +
      scale_x_continuous(limits = c(min(plotDF$x), max(plotDF$x)), expand = c(0,0))
    pdf(file.path(QCDir,paste0(samples[i],"-TSS_Enrichment(filtered).pdf")),
        width = 4.5, height = 4.5)
    .fixPlotSize(p, plotWidth = 6, plotHeight = 6)
    dev.off()
  }
}

proj_atac_filter <- subsetArchRProject(ArchRProj = proj_all,
                                       cells = Cells,
                                       outputDirectory = "ATAC_Filter",
                                       dropCells = TRUE,
                                       force = TRUE)

samples <- unique(proj_atac_filter$Sample)
samples_pro <- list()
for (s in samples) {
  sample_cells <- proj_atac_filter$cellNames[c(proj_atac_filter$Sample == s)]
  temp_pro_sample <- proj_atac_filter[sample_cells,]
  temp_list <- list(temp_pro_sample)
  names(temp_list) <- s
  samples_pro <- c(samples_pro, temp_list)
}
for (s in samples) {
  temp_pro_sample <- samples_pro[[s]]
  temp_pro_sample <- addDoubletScores(
    input = temp_pro_sample,
    nTrials = 20,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
    LSIMethod = 1,
    dimsToUse = 1:30
  )
  samples_pro[[s]] <- temp_pro_sample
}
for (s in samples) {
  temp_pro_sample <- samples_pro[[s]]
  df <- getCellColData(temp_pro_sample)
  df <- as.data.frame(df)
  df$Cells <- rownames(df)
  df <- df[order(df$DoubletEnrichment, decreasing = TRUE), ]
  n <- nrow(df)
  cutEnrich <- 1
  cutScore <- -Inf
  filterRatio <- 1
  if(!is.null(cutEnrich)){
    df <- df[which(df$DoubletEnrichment >= cutEnrich), ]
  }
  if(!is.null(cutScore)){
    df <- df[which(df$DoubletScore >= cutScore), ]
  }
  temp_doublet <- head(rownames(df), filterRatio * n * (n / 100000))
  temp_pro_sample$DoubletIdent <- ifelse(temp_pro_sample$cellNames %in% temp_doublet,
                                         "Doublet", "Singlet")
  samples_pro[[s]] <- temp_pro_sample
}

for (s in samples) {
  temp_pro_sample <- samples_pro[[s]]
  temp_pro_sample <- addIterativeLSI(ArchRProj = temp_pro_sample,
                                     useMatrix = "TileMatrix",
                                     name = "IterativeLSI")
  temp_pro_sample <- addUMAP(
    ArchRProj = temp_pro_sample, 
    reducedDims = "IterativeLSI", 
    name = "UMAP_IterativeLSI", 
    nNeighbors = 40, 
    minDist = 0.4, 
    metric = "cosine"
  )
  samples_pro[[s]] <- temp_pro_sample
}

for (s in samples) {
  temp_pro_sample <- samples_pro[[s]]
  df <- getEmbedding(temp_pro_sample, embedding = "UMAP_IterativeLSI")
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$DoubletIdent <- as.data.frame(getCellColData(temp_pro_sample))[rownames(df),"DoubletIdent"]
  df <- df[order(df$DoubletIdent, decreasing = T),]
  
  p <- ggplot(data = df, aes(x = UMAP_1, y = UMAP_2, color = DoubletIdent)) +
    geom_point(size = 0.4) +
    theme_classic() +
    ggtitle(paste0(s, " | ", as.numeric(table(df$DoubletIdent)["Doublet"]), " doublets")) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    scale_color_manual(values = c("red", "black"))
  pdf(paste0("./ArchR_ATAC/1.QCplot/Doublet/1.ATAC_", s, ".pdf"), width = 5.5, height = 4.5)
  print(p)
  dev.off()
}
samples_pro # ATAC

doublets_cells <- c()
for (s in samples) {
  temp_pro_sample <- samples_pro[[s]]
  df <- getCellColData(temp_pro_sample)
  df <- as.data.frame(df)
  doublets_cells <- c(doublets_cells,
                      rownames(df)[df$DoubletIdent == "Doublet"])
}

proj_atac_filter$DoubletIdent <- ifelse(proj_atac_filter$cellNames %in% doublets_cells,
                                        "Doublet", "Singlet")
table(proj_atac_filter$DoubletIdent)
# 10,574个doublet
singlet <- proj_atac_filter$cellNames[proj_atac_filter$DoubletIdent == "Singlet"]
saveRDS(proj_atac_filter, "ATAC_2.pass_filteration.rds")

ATAC_singlet <- subsetArchRProject(ArchRProj = proj_atac_filter,
                                   cells = singlet,
                                   outputDirectory = "ATAC_Singlet",
                                   dropCells = TRUE,
                                   force = TRUE)
saveRDS(ATAC_singlet, "ATAC_3.pass_filteration_singlet.rds")

RNA_singlet <- readRDS("sce.0.2.integrate_umap_last.3000.ham.RDS")
DimPlot(RNA_singlet, reduction = "umap", label = T)
RNA_singlet <- RenameIdents(RNA_singlet,
                            "Antler_Mural cells" = "Mural cells",
                            "Antler_proliferative_progenitor cells" = "Proliferative_progenitor cells",
                            "Antler_progenitor cells" = "Progenitor cells",
                            "Antler_Endothelial cells" = "Endothelial cells",
                            "Antler_Chondrocytes1" = "Chondrocytes",
                            "Antler_AnMCs" = "AnMCs",
                            "Antler_Monocytes_Macrophages" = "Monocytes_Macrophages",
                            "Antler_T cells" = "Mast cells",
                            "Antler_Chondroclasts" = "Chondroclasts",
                            "Antler_Hypertrophic chondrocytes" = "Hypertrophic chondrocytes"
)
RNA_singlet$Celltype_rename <- factor(as.character(Idents(RNA_singlet)),
                                      levels = c("AnMCs", "Proliferative_progenitor cells",
                                                 "Progenitor cells", "Chondrocytes", "Hypertrophic chondrocytes",
                                                 "Chondroclasts", "Mural cells", "Endothelial cells",
                                                 "Monocytes_Macrophages", "Mast cells"))
Idents(RNA_singlet) <- RNA_singlet$Celltype_rename
DimPlot(RNA_singlet)
RNA_Celltype_color <- c("#A4C9DD", "#2572A9", "#ADD487", "#399938",
                        "#F19695", "#D5231E", "#F5BB6F", "#EF7C1C",
                        "#C6B0D2", "#653B90")
names(RNA_Celltype_color) <- levels(RNA_singlet$Celltype_rename)

DimPlot(RNA_singlet, cols = RNA_Celltype_color)
RNA_singlet$Group <- unlist(lapply(RNA_singlet$sample, function(x){
  unlist(strsplit(x, ".", fixed = T))[2]
}))
table(RNA_singlet$Group)
RNA_group_color <- c("#E7211A", "#EFEA3C", "#72C8D5", "#6AB82D", "#18499E")
names(RNA_group_color) <- c("RM", "PC", "TZ", "CA", "MC")

ATAC_singlet <- readRDS("ATAC_3.pass_filteration_singlet.rds")
ATAC_singlet$Group <- unlist(lapply(ATAC_singlet$Sample,
                                    function(x){
                                      temp <- unlist(strsplit(x, "", fixed = T))
                                      paste0(temp[1:2], collapse = "")
                                    }))

ATAC_singlet <- addIterativeLSI(ArchRProj = ATAC_singlet,
                                useMatrix = "TileMatrix",
                                name = "IterativeLSI",
                                force = TRUE)
ATAC_singlet <- addHarmony(
  ArchRProj = ATAC_singlet,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = T
)
ATAC_singlet <- addClusters(
  input = ATAC_singlet,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Clusters_0.5",
  resolution = 0.5,
  force = TRUE,
  dimsToUse = 1:30,
  maxClusters = 50
)
ATAC_singlet <- addUMAP(
  ArchRProj = ATAC_singlet, 
  reducedDims = "IterativeLSI", 
  name = "UMAP_withBatchEffect", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)
ATAC_singlet <- addUMAP(
  ArchRProj = ATAC_singlet, 
  reducedDims = "Harmony", 
  name = "UMAP_Harmony", 
  nNeighbors = 30, 
  minDist = 0.2, 
  metric = "cosine",
  force = T,
  n_neighbors = 20,
  local_connectivity = 1.1
)

ATAC_singlet <- addTSNE(
  ArchRProj = ATAC_singlet, 
  reducedDims = "Harmony", 
  name = "TSNE_Harmony", 
  dimsToUse = 1:30
)
ATAC_singlet <- addTSNE(
  ArchRProj = ATAC_singlet, 
  reducedDims = "IterativeLSI", 
  name = "TSNE_withBatchEffect", 
  dimsToUse = 1:30
)

cluster_num <- length(unique(ATAC_singlet$Harmony_Clusters_1.2))
ATAC_singlet$Harmony_Clusters_1.2 <- factor(ATAC_singlet$Harmony_Clusters_1.2,
                                            levels = paste0("C", 1:cluster_num))

set.seed(123)
ATAC_Cluster_color <- randomcoloR::distinctColorPalette(k = cluster_num)
names(ATAC_Cluster_color) <- paste0("C", 1:cluster_num)
pdf("./ATAC_Singlet/ATAC_1.Cluster_TSNE.pdf",
    height = 4.2, width = 5.7)
ArchR_DimPlot(ArchR = ATAC_singlet, reduction_dimname = c("TSNE_1", "TSNE_2"),
              reduction = "TSNE_Harmony", color_by = "Harmony_Clusters_1.2", pt.size = 0.1,
              label = TRUE, label_size = 3, color = ATAC_Cluster_color, repel = FALSE) +
  ggtitle("ATAC_TSNE_Harmony") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  labs(color = "Cluster")
dev.off()

pdf("./ATAC_Singlet/ATAC_1.Cluster_UMAP.pdf",
    height = 4.2, width = 5.7)
ArchR_DimPlot(ArchR = ATAC_singlet, reduction_dimname = c("UMAP_1", "UMAP_2"),
              reduction = "UMAP_Harmony", color_by = "Harmony_Clusters_1.2", pt.size = 0.1,
              label = TRUE, label_size = 3, color = ATAC_Cluster_color, repel = FALSE) +
  ggtitle("ATAC_UMAP_Harmony") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  labs(color = "Cluster")
dev.off()

pdf("./ATAC_Singlet/ATAC_2.Group_UMAP.pdf",
    height = 4.2, width = 5.2)
ArchR_DimPlot(ArchR = ATAC_singlet, reduction_dimname = c("UMAP_1", "UMAP_2"),
              reduction = "UMAP_Harmony", color_by = "Group", pt.size = 0.1,
              label = FALSE, label_size = 4, color = RNA_group_color, repel = FALSE) +
  ggtitle("ATAC_UMAP_Harmony") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  labs(color = "Group")
dev.off()

pdf("./ATAC_Singlet/ATAC_2.Group_TSNE.pdf",
    height = 4.2, width = 5.2)
ArchR_DimPlot(ArchR = ATAC_singlet, reduction_dimname = c("TSNE_1", "TSNE_2"),
              reduction = "TSNE_Harmony", color_by = "Group", pt.size = 0.1,
              label = FALSE, label_size = 4, color = RNA_group_color, repel = FALSE) +
  ggtitle("ATAC_TSNE_Harmony") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  labs(color = "Group")
dev.off()