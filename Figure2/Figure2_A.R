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
source("./monocle3_util.R")
#### 参数设置
dim_use <- 30
groups <- "Celltype_rename"
batch <- "orig.ident"
reduction_method <- "UMAP"

#### 运行 Monocle3_Step1 构建轨迹
data <- GetAssayData(RNA_singlet_sub, assay = 'RNA', slot = 'counts')
cell_metadata <- RNA_singlet_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
monocle3_cds <- new_cell_data_set(data,
                                  cell_metadata = cell_metadata,
                                  gene_metadata = gene_annotation)
monocle3_cds <- preprocess_cds(monocle3_cds, num_dim = dim_use)
monocle3_cds <- align_cds(monocle3_cds, alignment_group = batch,
                          alignment_k = dim_use)
monocle3_cds <- reduce_dimension(monocle3_cds,
                                 reduction_method = "UMAP",
                                 preprocess_method = "Aligned",
                                 umap.min_dist = 0.1,
                                 umap.n_neighbors = 30)
cds.embed <- monocle3_cds@int_colData$reducedDims$UMAP
colnames(cds.embed) <- c("UMAP_1", "UMAP_2")
cds.embed <- as.data.frame(cds.embed)
cds.embed$Celltype <- RNA_singlet_sub@meta.data[rownames(cds.embed), "Celltype_rename"]
ggplot(data = cds.embed, aes(x = UMAP_1, y = UMAP_2, colour = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = RNA_Celltype_color)
monocle3_cds <- cluster_cells(monocle3_cds,
                              reduction_method = "UMAP",
                              k = dim_use)
int.embed <- Embeddings(RNA_singlet_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
monocle3_cds@int_colData$reducedDims$UMAP <- int.embed

monocle3_cds <- learn_graph(monocle3_cds,
                            use_partition = F,
                            close_loop = T,
                            learn_graph_control = list(nn.k = 15, # 15
                                                       minimal_branch_len = 16)) # 16
#### 可视化
# 可视化轨迹
umap <- Embeddings(RNA_singlet, reduction = "umap")
umap <- as.data.frame(umap)
p <- plot_cells_2(cds = monocle3_cds,
                  reduction_method = reduction_method, # 可以选择TSNE， PCA
                  color_cells_by = groups, # 展示细胞的分组，这里设置的是Sample
                  show_trajectory_graph = TRUE, trajectory_graph_color = "grey28", trajectory_graph_segment_size = 1.5, # 展示轨迹
                  group_cells_color = RNA_Celltype_color, # 设置细胞分组的颜色，默认是NULL，可以自己设置
                  label_cell_groups = FALSE, label_groups_by_cluster = TRUE, # 设置是否展示细胞分组
                  group_label_size = 6, # 设置分组label的大小
                  label_principal_points = FALSE, label_leaves = FALSE, graph_label_size = 0, # 设置轨迹上的分至点以及叶节点
                  cell_size = 0.8, # 设置细胞大小
                  background_cell = umap,
                  background_cell_size = 0.8
)
pdf(paste0("All_cells_Monocle3_Step1_", reduction_method, "_", groups, ".pdf"),
    width = 8, height = 5)
print(p)
dev.off()

monocle3_cds <- Monocle3_Step2(CDS = monocle3_cds, # cds 数据对象，上一步得到的结果
                               reduction_method = reduction_method, # 降维方法
                               root_cells = NULL # 手动设置根细胞，默认为NULL，即用鼠标圈选。如果不想用鼠标圈，可以直接输入细胞的barcode，比如 TTTCACAGTAGTCGTT-2，可以输入多个细胞作为起点
)
#### 可视化拟时序
p <- plot_cells_2(cds = monocle3_cds,
                  reduction_method = reduction_method, # 可以选择TSNE， PCA
                  color_cells_by = "pseudotime", # 展示细胞的分组，这里设置的是Sample
                  show_trajectory_graph = TRUE, trajectory_graph_color = "grey28", trajectory_graph_segment_size = 1.5, # 展示轨迹
                  group_cells_color = NULL, # 设置细胞分组的颜色，默认是NULL，可以自己设置
                  label_cell_groups = FALSE, label_groups_by_cluster = TRUE, # 设置是否展示细胞分组
                  group_label_size = 6, # 设置分组label的大小
                  label_principal_points = FALSE, label_leaves = FALSE, graph_label_size = 0, # 设置轨迹上的分至点以及叶节点
                  cell_size = 0.8, # 设置细胞大小
                  background_cell = umap,
                  background_cell_size = 0.8
)
pdf(paste0("All_cells_Monocle3_Step2_",
           reduction_method, "_pseudotime", ".pdf"),
    width = 6.8, height = 5)
print(p)
dev.off()

### 
library(clusterProfiler)
library(AUCell)
monocle3_pseudotime <- pseudotime(monocle3_cds)
{
  gene_files <- list.files("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/09.geneset/blood vessel")
  GeneList <- c()
  for (i in gene_files) {
    temp <- read.gmt(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/09.geneset/blood vessel/",i))
    GeneList <- as.data.frame(rbind(GeneList,
                                    temp))
  }
  GeneList <- split.data.frame(GeneList, f = list(GeneList$term))
  GeneList <- lapply(GeneList, function(x){
    x[,2]
  })
  length(GeneList)
  cells_AUC <- AUCell_run(RNA_singlet_sub@assays$RNA@counts, GeneList)
  cells_AUC_df <- cells_AUC@assays@data@listData[["AUC"]]
  cells_AUC_df <- data.frame(t(cells_AUC_df),
                             check.rows = F, check.names = F)
  monocle3_pseudotime <- monocle3_pseudotime[rownames(cells_AUC_df)]
  monocle3_pseudotime_cor <- c()
  for (i in 1:ncol(cells_AUC_df)) {
    print(i)
    cor_result <- WGCNA::corAndPvalue(monocle3_pseudotime,
                                      cells_AUC_df[,i])
    monocle3_pseudotime_cor <- as.data.frame(rbind(monocle3_pseudotime_cor,
                                                   data.frame(Pathway = colnames(cells_AUC_df)[i],
                                                              Cor = cor_result$cor,
                                                              P = cor_result$p)))
  }
  
  
  monocle3_pseudotime_cor$Pathway <- unlist(lapply(monocle3_pseudotime_cor$Pathway,
                                                   function(x){
                                                     x <- tolower(x)
                                                   }))
  
  umap_sub <- Embeddings(RNA_singlet_sub, reduction = "umap")
  umap_sub <- as.data.frame(umap_sub)
  df <- data.frame("Celltype" = RNA_singlet_sub@meta.data[names(monocle3_pseudotime), "Celltype_rename"],
                   "Pseudotime" = monocle3_pseudotime,
                   cells_AUC_df
  )
  df <- as.data.frame(cbind(umap_sub,
                            df))
  max(df$Pseudotime)
  intervals <- cut(df$Pseudotime, breaks = seq(from = 0,
                                               to = 28, by = 2),
                   include.lowest = TRUE)
  df$Pseudotime_bin <- intervals
  for (i in colnames(cells_AUC_df)) {
    pathway <- unlist(strsplit(i, split = "_", fixed = T))
    pathway <- pathway[-1]
    pathway <- unlist(lapply(pathway, function(x){
      tolower(x)
    }))
    pathway <- paste0(pathway, collapse = " ")
    pathway <- Hmisc::capitalize(pathway)
    pdf(paste0("./轨迹打分/blood vessel/All_cells_Monocle3_", pathway, ".pdf"),
        width = 4.5, height = 4)
    p <- jitter_linebox(meta = df, group.by = "Pseudotime_bin",
                        features = i,
                        group.color = viridis::viridis(n = length(levels(df$Pseudotime_bin)),
                                                       option = "C"),
                        linewidth = 2, jitter_width = 0.2, stat_cor_size = 4) +
      theme(legend.position = "none",
            plot.title = element_text(size = 10)) +
      labs(x = "Pseudotime bin", y = "AUCell score") +
      ggtitle(pathway)
    print(p)
    dev.off()
  }
  
}
{
  gene_files <- list.files("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/09.geneset/hallmark")
  GeneList <- c()
  for (i in gene_files) {
    temp <- read.gmt(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/09.geneset/hallmark/",i))
    GeneList <- as.data.frame(rbind(GeneList,
                                    temp))
  }
  GeneList <- split.data.frame(GeneList, f = list(GeneList$term))
  GeneList <- lapply(GeneList, function(x){
    x[,2]
  })
  length(GeneList)
  cells_AUC <- AUCell_run(RNA_singlet_sub@assays$RNA@counts, GeneList)
  cells_AUC_df <- cells_AUC@assays@data@listData[["AUC"]]
  cells_AUC_df <- data.frame(t(cells_AUC_df),
                             check.rows = F, check.names = F)
  monocle3_pseudotime <- monocle3_pseudotime[rownames(cells_AUC_df)]
  monocle3_pseudotime_cor <- c()
  for (i in 1:ncol(cells_AUC_df)) {
    print(i)
    cor_result <- WGCNA::corAndPvalue(monocle3_pseudotime,
                                      cells_AUC_df[,i])
    monocle3_pseudotime_cor <- as.data.frame(rbind(monocle3_pseudotime_cor,
                                                   data.frame(Pathway = colnames(cells_AUC_df)[i],
                                                              Cor = cor_result$cor,
                                                              P = cor_result$p)))
  }
  
  
  monocle3_pseudotime_cor$Pathway <- unlist(lapply(monocle3_pseudotime_cor$Pathway,
                                                   function(x){
                                                     x <- tolower(x)
                                                   }))
  
  umap_sub <- Embeddings(RNA_singlet_sub, reduction = "umap")
  umap_sub <- as.data.frame(umap_sub)
  df <- data.frame("Celltype" = RNA_singlet_sub@meta.data[names(monocle3_pseudotime), "Celltype_rename"],
                   "Pseudotime" = monocle3_pseudotime,
                   cells_AUC_df
  )
  df <- as.data.frame(cbind(umap_sub,
                            df))
  max(df$Pseudotime)
  intervals <- cut(df$Pseudotime, breaks = seq(from = 0,
                                               to = 28, by = 2),
                   include.lowest = TRUE)
  df$Pseudotime_bin <- intervals
  for (i in colnames(cells_AUC_df)) {
    pathway <- unlist(strsplit(i, split = "_", fixed = T))
    pathway <- pathway[-1]
    pathway <- unlist(lapply(pathway, function(x){
      tolower(x)
    }))
    pathway <- paste0(pathway, collapse = " ")
    pathway <- Hmisc::capitalize(pathway)
    pdf(paste0("./轨迹打分/hallmark/All_cells_Monocle3_", pathway, ".pdf"),
        width = 4.5, height = 4)
    p <- jitter_linebox(meta = df, group.by = "Pseudotime_bin",
                        features = i,
                        group.color = viridis::viridis(n = length(levels(df$Pseudotime_bin)),
                                                       option = "C"),
                        linewidth = 2, jitter_width = 0.2, stat_cor_size = 4) +
      theme(legend.position = "none",
            plot.title = element_text(size = 10)) +
      labs(x = "Pseudotime bin", y = "AUCell score") +
      ggtitle(pathway)
    print(p)
    dev.off()
  }
  
}
{
  gene_files <- list.files("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/09.geneset/HYPOXIA")
  GeneList <- c()
  for (i in gene_files) {
    temp <- read.gmt(paste0("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous/09.geneset/HYPOXIA/",i))
    GeneList <- as.data.frame(rbind(GeneList,
                                    temp))
  }
  GeneList <- split.data.frame(GeneList, f = list(GeneList$term))
  GeneList <- lapply(GeneList, function(x){
    x[,2]
  })
  length(GeneList)
  cells_AUC <- AUCell_run(RNA_singlet_sub@assays$RNA@counts, GeneList)
  cells_AUC_df <- cells_AUC@assays@data@listData[["AUC"]]
  cells_AUC_df <- data.frame(t(cells_AUC_df),
                             check.rows = F, check.names = F)
  monocle3_pseudotime <- monocle3_pseudotime[rownames(cells_AUC_df)]
  monocle3_pseudotime_cor <- c()
  for (i in 1:ncol(cells_AUC_df)) {
    print(i)
    cor_result <- WGCNA::corAndPvalue(monocle3_pseudotime,
                                      cells_AUC_df[,i])
    monocle3_pseudotime_cor <- as.data.frame(rbind(monocle3_pseudotime_cor,
                                                   data.frame(Pathway = colnames(cells_AUC_df)[i],
                                                              Cor = cor_result$cor,
                                                              P = cor_result$p)))
  }
  
  
  monocle3_pseudotime_cor$Pathway <- unlist(lapply(monocle3_pseudotime_cor$Pathway,
                                                   function(x){
                                                     x <- tolower(x)
                                                   }))
  
  umap_sub <- Embeddings(RNA_singlet_sub, reduction = "umap")
  umap_sub <- as.data.frame(umap_sub)
  df <- data.frame("Celltype" = RNA_singlet_sub@meta.data[names(monocle3_pseudotime), "Celltype_rename"],
                   "Pseudotime" = monocle3_pseudotime,
                   cells_AUC_df
  )
  df <- as.data.frame(cbind(umap_sub,
                            df))
  max(df$Pseudotime)
  intervals <- cut(df$Pseudotime, breaks = seq(from = 0,
                                               to = 28, by = 2),
                   include.lowest = TRUE)
  df$Pseudotime_bin <- intervals
  for (i in colnames(cells_AUC_df)) {
    pathway <- unlist(strsplit(i, split = "_", fixed = T))
    pathway <- pathway[-1]
    pathway <- unlist(lapply(pathway, function(x){
      tolower(x)
    }))
    pathway <- paste0(pathway, collapse = " ")
    pathway <- Hmisc::capitalize(pathway)
    pdf(paste0("./轨迹打分/HYPOXIA/All_cells_Monocle3_", pathway, ".pdf"),
        width = 4.5, height = 4)
    p <- jitter_linebox(meta = df, group.by = "Pseudotime_bin",
                        features = i,
                        group.color = viridis::viridis(n = length(levels(df$Pseudotime_bin)),
                                                       option = "C"),
                        linewidth = 2, jitter_width = 0.2, stat_cor_size = 4) +
      theme(legend.position = "none",
            plot.title = element_text(size = 10)) +
      labs(x = "Pseudotime bin", y = "AUCell score") +
      ggtitle(pathway)
    print(p)
    dev.off()
  }
  
}