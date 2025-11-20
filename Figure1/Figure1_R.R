
all_gene_sets <- clusterProfiler::read.gmt("./打分补充/h.all.v2025.1.Hs.symbols.gmt")
all_gene_sets$term <- as.character(all_gene_sets$term)
all_gene_sets$term <- unlist(lapply(all_gene_sets$term,
                                    function(x){
                                      temp <- unlist(strsplit(x, "_", fixed = TRUE))
                                      temp <- temp[-1]
                                      temp <- tolower(temp)
                                      temp[1] <- Hmisc::capitalize(temp[1])
                                      paste0(temp, collapse = " ")
                                    }))
all_gene_sets <- all_gene_sets[all_gene_sets$term %in% c("Apoptosis",
                                                         "P53 pathway",
                                                         "Myc targets v1",
                                                         "Myc targets v2",
                                                         "Dna repair"),]
all_gene_sets <- split.data.frame(all_gene_sets,
                                  f = list(all_gene_sets$term))
all_gene_sets <- lapply(all_gene_sets, function(x){
  x[,2]
})
genes <- getFeatures(Prolif_PC, useMatrix = "GeneScoreMatrix")
all_gene_sets <- lapply(all_gene_sets, function(x){
  x[x %in% genes]
})
Prolif_PC <- addModuleScore(Prolif_PC,
                            useMatrix = "GeneScoreMatrix",
                            name = "Module",
                            features = all_gene_sets)
Prolif_PC_meta <- getCellColData(Prolif_PC)
Prolif_PC_meta <- as.data.frame(Prolif_PC_meta)

Prolif_PC_score <- data.frame(Pathway = c("Apoptosis"),
                              Score = Prolif_PC_meta$Module.Apoptosis)
Prolif_PC_score <- as.data.frame(rbind(Prolif_PC_score,
                                       data.frame(Pathway = c("P53 pathway"),
                                                  Score = Prolif_PC_meta$Module.P53.pathway)))
Prolif_PC_score <- as.data.frame(rbind(Prolif_PC_score,
                                       data.frame(Pathway = c("Myc targets v1"),
                                                  Score = Prolif_PC_meta$Module.Myc.targets.v1)))
Prolif_PC_score <- as.data.frame(rbind(Prolif_PC_score,
                                       data.frame(Pathway = c("Myc targets v2"),
                                                  Score = Prolif_PC_meta$Module.Myc.targets.v2)))
Prolif_PC_score <- as.data.frame(rbind(Prolif_PC_score,
                                       data.frame(Pathway = c("Dna repair"),
                                                  Score = Prolif_PC_meta$Module.Dna.repair)))
Prolif_PC_score$Pathway <- factor(Prolif_PC_score$Pathway,
                                  levels = c("Apoptosis", "P53 pathway",
                                             "Myc targets v1", "Myc targets v2",
                                             "Dna repair"))
ggplot(data = Prolif_PC_score,
       aes(x = Pathway, y = Score)) +
  geom_violin()
Prolif_PC <- addImputeWeights(Prolif_PC)
Prolif_PC_imputeWeights <- getImputeWeights(Prolif_PC)

Prolif_PC_score_MAGIC <- c()
for (i in c("Apoptosis", "P53 pathway",
            "Myc targets v1", "Myc targets v2",
            "Dna repair")) {
  colorMat <- matrix(Prolif_PC_score$Score[Prolif_PC_score$Pathway == i], nrow = 1)
  colnames(colorMat) <- rownames(Prolif_PC_meta)
  colorMat <- imputeMatrix(mat = colorMat,
                           imputeWeights = Prolif_PC_imputeWeights)
  Prolif_PC_score_MAGIC <- as.data.frame(rbind(Prolif_PC_score_MAGIC,
                                               data.frame(Pathway = i,
                                                          Score = colorMat)))
}
Prolif_PC_score_MAGIC$Pathway <- factor(Prolif_PC_score_MAGIC$Pathway,
                                        levels = c("Apoptosis", "P53 pathway",
                                                   "Myc targets v1", "Myc targets v2",
                                                   "Dna repair"))


g <- violin_linebox(
  meta = Prolif_PC_score_MAGIC,
  group.by = Pathway,   # 不加引号！
  features = Score,
  group.color = c("#EF7C1C", "#F5BB6F", "#C6B0D2", "#A4C9DD", "#F19695"),
  linewidth = 2,
  compare_list = list(
    c("Apoptosis", "Myc targets v1"),
    c("Apoptosis", "Myc targets v2"),
    c("Apoptosis", "Dna repair"),
    c("P53 pathway", "Myc targets v1"),
    c("P53 pathway", "Myc targets v2"),
    c("P53 pathway", "Dna repair"))
)
g <- g + ggtitle("Proliferative PCs") + NoLegend()
pdf("./打分补充/Proliferative PCs.pdf", width = 3, height = 5)
g
dev.off()




PC <- addModuleScore(PC,
                     useMatrix = "GeneScoreMatrix",
                     name = "Module",
                     features = all_gene_sets)
PC_meta <- getCellColData(PC)
PC_meta <- as.data.frame(PC_meta)

PC_score <- data.frame(Pathway = c("Apoptosis"),
                       Score = PC_meta$Module.Apoptosis)
PC_score <- as.data.frame(rbind(PC_score,
                                data.frame(Pathway = c("P53 pathway"),
                                           Score = PC_meta$Module.P53.pathway)))
PC_score <- as.data.frame(rbind(PC_score,
                                data.frame(Pathway = c("Myc targets v1"),
                                           Score = PC_meta$Module.Myc.targets.v1)))
PC_score <- as.data.frame(rbind(PC_score,
                                data.frame(Pathway = c("Myc targets v2"),
                                           Score = PC_meta$Module.Myc.targets.v2)))
PC_score <- as.data.frame(rbind(PC_score,
                                data.frame(Pathway = c("Dna repair"),
                                           Score = PC_meta$Module.Dna.repair)))


PC_score$Pathway <- factor(PC_score$Pathway,
                           levels = c("Apoptosis", "P53 pathway",
                                      "Myc targets v1", "Myc targets v2",
                                      "Dna repair"))

PC <- addImputeWeights(PC)
PC_imputeWeights <- getImputeWeights(PC)

PC_score_MAGIC <- c()
for (i in c("Apoptosis", "P53 pathway",
            "Myc targets v1", "Myc targets v2",
            "Dna repair")) {
  colorMat <- matrix(PC_score$Score[PC_score$Pathway == i], nrow = 1)
  colnames(colorMat) <- rownames(PC_meta)
  colorMat <- imputeMatrix(mat = colorMat,
                           imputeWeights = PC_imputeWeights)
  PC_score_MAGIC <- as.data.frame(rbind(PC_score_MAGIC,
                                        data.frame(Pathway = i,
                                                   Score = colorMat)))
}
PC_score_MAGIC$Pathway <- factor(PC_score_MAGIC$Pathway,
                                 levels = c("Apoptosis", "P53 pathway",
                                            "Myc targets v1", "Myc targets v2",
                                            "Dna repair"))


g <- violin_linebox(
  meta = PC_score_MAGIC,
  group.by = Pathway,   # 不加引号！
  features = Score,
  group.color = c("#EF7C1C", "#F5BB6F", "#C6B0D2", "#A4C9DD", "#F19695"),
  linewidth = 2,
  compare_list = list(
    c("Apoptosis", "Myc targets v1"),
    c("Apoptosis", "Myc targets v2"),
    c("Apoptosis", "Dna repair"),
    c("P53 pathway", "Myc targets v1"),
    c("P53 pathway", "Myc targets v2"),
    c("P53 pathway", "Dna repair"))
)
g <- g + ggtitle("PCs") + NoLegend()
pdf("./打分补充/PCs.pdf", width = 3, height = 5)
g
dev.off()





#############


all_gene_sets <- clusterProfiler::read.gmt("../../打分补充/h.all.v2025.1.Hs.symbols.gmt")
all_gene_sets$term <- as.character(all_gene_sets$term)
all_gene_sets$term <- unlist(lapply(all_gene_sets$term,
                                    function(x){
                                      temp <- unlist(strsplit(x, "_", fixed = TRUE))
                                      temp <- temp[-1]
                                      temp <- tolower(temp)
                                      temp[1] <- Hmisc::capitalize(temp[1])
                                      paste0(temp, collapse = " ")
                                    }))
all_gene_sets <- all_gene_sets[all_gene_sets$term %in% c("Apoptosis",
                                                         "P53 pathway",
                                                         "Myc targets v1",
                                                         "Myc targets v2",
                                                         "Dna repair"),]
all_gene_sets <- split.data.frame(all_gene_sets,
                                  f = list(all_gene_sets$term))
all_gene_sets <- lapply(all_gene_sets, function(x){
  x[,2]
})
genes <- getFeatures(Project, useMatrix = "GeneScoreMatrix")
all_gene_sets <- lapply(all_gene_sets, function(x){
  x[x %in% genes]
})
Project <- addModuleScore(Project,
                          useMatrix = "GeneScoreMatrix",
                          name = "Module",
                          features = all_gene_sets)
Project_meta <- getCellColData(Project)
Project_meta <- as.data.frame(Project_meta)

Project_score <- data.frame(Pathway = c("Apoptosis"),
                            Score = Project_meta$Module.Apoptosis)
Project_score <- as.data.frame(rbind(Project_score,
                                     data.frame(Pathway = c("P53 pathway"),
                                                Score = Project_meta$Module.P53.pathway)))
Project_score <- as.data.frame(rbind(Project_score,
                                     data.frame(Pathway = c("Myc targets v1"),
                                                Score = Project_meta$Module.Myc.targets.v1)))
Project_score <- as.data.frame(rbind(Project_score,
                                     data.frame(Pathway = c("Myc targets v2"),
                                                Score = Project_meta$Module.Myc.targets.v2)))
Project_score <- as.data.frame(rbind(Project_score,
                                     data.frame(Pathway = c("Dna repair"),
                                                Score = Project_meta$Module.Dna.repair)))
Project_score$Pathway <- factor(Project_score$Pathway,
                                levels = c("Apoptosis", "P53 pathway",
                                           "Myc targets v1", "Myc targets v2",
                                           "Dna repair"))
ggplot(data = Project_score,
       aes(x = Pathway, y = Score)) +
  geom_violin()
Project <- addImputeWeights(Project)
Project_imputeWeights <- getImputeWeights(Project)

Project_score_MAGIC <- c()
for (i in c("Apoptosis", "P53 pathway",
            "Myc targets v1", "Myc targets v2",
            "Dna repair")) {
  colorMat <- matrix(Project_score$Score[Project_score$Pathway == i], nrow = 1)
  colnames(colorMat) <- rownames(Project_meta)
  colorMat <- imputeMatrix(mat = colorMat,
                           imputeWeights = Project_imputeWeights)
  Project_score_MAGIC <- as.data.frame(rbind(Project_score_MAGIC,
                                             data.frame(Pathway = i,
                                                        Score = colorMat)))
}
Project_score_MAGIC$Pathway <- factor(Project_score_MAGIC$Pathway,
                                      levels = c("Apoptosis", "P53 pathway",
                                                 "Myc targets v1", "Myc targets v2",
                                                 "Dna repair"))


g <- violin_linebox(
  meta = Project_score_MAGIC,
  group.by = Pathway,   # 不加引号！
  features = Score,
  group.color = c("#EF7C1C", "#F5BB6F", "#C6B0D2", "#A4C9DD", "#F19695"),
  linewidth = 2,
  compare_list = list(
    c("Apoptosis", "Myc targets v1"),
    c("Apoptosis", "Myc targets v2"),
    c("Apoptosis", "Dna repair"),
    c("P53 pathway", "Myc targets v1"),
    c("P53 pathway", "Myc targets v2"),
    c("P53 pathway", "Dna repair"))
)
g <- g + ggtitle("Osteosarcoma cells") + NoLegend()
pdf("./通路打分.pdf", width = 3, height = 5)
g
dev.off()