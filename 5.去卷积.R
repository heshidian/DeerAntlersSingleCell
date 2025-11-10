library(BayesPrism)
library(Seurat)
library(MuSiC)
library(ggplot2)

##### BayesPrism
{
  setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
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
  
  bulk <- read.csv("antler-bulkrna.matrix.gene", sep = "\t",
                   header = TRUE, row.names = 1)
  sum(rownames(bulk) %in% rownames(RNA_singlet))
  bulk <- t(bulk)
  
  RNA_singlet_exp <- RNA_singlet@assays$RNA@counts
  RNA_singlet_exp <- t(RNA_singlet_exp)
  
  common_genes <- intersect(colnames(RNA_singlet_exp),
                            colnames(bulk))
  RNA_singlet_exp <- RNA_singlet_exp[,common_genes]
  bulk <- bulk[,common_genes]
  
  cell.type.labels <- as.character(RNA_singlet$Celltype_rename)
  names(cell.type.labels) <- colnames(RNA_singlet)
  cell.type.labels <- cell.type.labels[rownames(RNA_singlet_exp)]
  
  plot.cor.phi(input = RNA_singlet_exp,
               input.labels = cell.type.labels,
               title="cell type correlation",
               #specify pdf.prefix if need to output to pdf
               #pdf.prefix="gbm.cor.cs",
               cexRow=0.7, cexCol=0.7,
               margins=c(2,2))
  
  sc.stat <- plot.scRNA.outlier(
    input=RNA_singlet_exp, #make sure the colnames are gene symbol or ENSMEBL ID
    cell.type.labels=cell.type.labels,
    species="hs", #currently only human(hs) and mouse(mm) annotations are supported
    return.raw=TRUE #return the data used for plotting.
    #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
  )
  
  bk.stat <- plot.bulk.outlier(
    bulk.input=bulk,#make sure the colnames are gene symbol or ENSMEBL ID
    sc.input=RNA_singlet_exp, #make sure the colnames are gene symbol or ENSMEBL ID
    cell.type.labels=cell.type.labels,
    species="hs", #currently only human(hs) and mouse(mm) annotations are supported
    return.raw=TRUE
    #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
  )
  
  sc.dat.filtered <- cleanup.genes (input=RNA_singlet_exp,
                                    input.type="count.matrix",
                                    species="hs",
                                    gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1") ,
                                    exp.cells=1)
  plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                   bulk.input = bulk
                   #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
  )
  
  sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,
                                         gene.type = "protein_coding")
  
  myPrism <- new.prism(
    reference=sc.dat.filtered.pc,
    mixture=bulk,
    input.type="count.matrix",
    cell.type.labels = cell.type.labels,
    cell.state.labels = cell.type.labels,
    key=NULL,
    outlier.cut=0.01,
    outlier.fraction=0.1,
  )
  bp.res <- run.prism(prism = myPrism, n.cores=50)
  theta <- get.fraction (bp=bp.res,
                         which.theta="final",
                         state.or.type="type")
  theta <- data.frame(t(theta),
                      check.rows = FALSE, check.names = FALSE)
  saveRDS(theta, "BayesPrism_result.rds")
}

##### MuSiC
{
  music_prop <- function (bulk.mtx, sc.sce, markers = NULL, clusters, samples, 
                          select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE, 
                          iter.max = 1000, nu = 1e-04, eps = 0.01, centered = FALSE, 
                          normalize = FALSE, ...) {
    bulk.gene = rownames(bulk.mtx)[rowMeans(bulk.mtx) != 0]
    bulk.mtx = bulk.mtx[bulk.gene, ]
    if (is.null(markers)) {
      sc.markers = bulk.gene
    } else {
      sc.markers = intersect(bulk.gene, unlist(markers))
    }
    sc.basis = music_basis(sc.sce, non.zero = TRUE, markers = sc.markers, 
                           clusters = clusters, samples = samples, select.ct = select.ct, 
                           cell_size = cell_size, ct.cov = ct.cov, verbose = verbose)
    message("music_basis Done!")
    cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
    if (is.null(markers)) {
      if (length(cm.gene) < 0.2 * min(length(bulk.gene), nrow(sc.sce))) 
        stop("Too few common genes!")
    } else {
      if (length(cm.gene) < 0.2 * length(unlist(markers))) 
        stop("Too few common genes!")
    }
    if (verbose) {
      message(paste("Used", length(cm.gene), "common genes..."))
    }
    m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx))
    m.bulk = match(cm.gene, bulk.gene)
    D1 = sc.basis$Disgn.mtx[m.sc, ]
    M.S = colMeans(sc.basis$S, na.rm = T)
    if (!is.null(cell_size)) {
      if (!is.data.frame(cell_size)) {
        stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
      } else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
        stop("Cell type names in cell_size must match clusters")
      } else if (any(is.na(as.numeric(cell_size[, 2])))) {
        stop("Cell sizes should all be numeric")
      }
      my_ms_names <- names(M.S)
      cell_size <- cell_size[my_ms_names %in% cell_size[, 1], 
      ]
      M.S <- cell_size[match(my_ms_names, cell_size[, 1]), 
      ]
      M.S <- M.S[, 2]
      names(M.S) <- my_ms_names
    }
    Yjg = relative.ab(bulk.mtx[m.bulk, ])
    N.bulk = ncol(bulk.mtx)
    if (ct.cov) {
      Sigma.ct = sc.basis$Sigma.ct[, m.sc]
      Est.prop.allgene = NULL
      Est.prop.weighted = NULL
      Weight.gene = NULL
      r.squared.full = NULL
      Var.prop = NULL
      for (i in 1:N.bulk) {
        if (sum(Yjg[, i] == 0) > 0) {
          D1.temp = D1[Yjg[, i] != 0, ]
          Yjg.temp = Yjg[Yjg[, i] != 0, i]
          Sigma.ct.temp = Sigma.ct[, Yjg[, i] != 0]
          if (verbose) 
            message(paste(colnames(Yjg)[i], "has common genes", 
                          sum(Yjg[, i] != 0), "..."))
        } else {
          D1.temp = D1
          Yjg.temp = Yjg[, i]
          Sigma.ct.temp = Sigma.ct
          if (verbose) 
            message(paste(colnames(Yjg)[i], "has common genes", 
                          sum(Yjg[, i] != 0), "..."))
        }
        lm.D1.weighted = music.iter.ct(Yjg.temp, D1.temp, 
                                       M.S, Sigma.ct.temp, iter.max = iter.max, nu = nu, 
                                       eps = eps, centered = centered, normalize = normalize)
        Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
        Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
        weight.gene.temp = rep(NA, nrow(Yjg))
        weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
        Weight.gene = cbind(Weight.gene, weight.gene.temp)
        r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
        Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
      }
    } else {
      Sigma = sc.basis$Sigma[m.sc, ]
      valid.ct = (colSums(is.na(Sigma)) == 0) & (colSums(is.na(D1)) == 
                                                   0) & (!is.na(M.S))
      if (sum(valid.ct) <= 1) {
        stop("Not enough valid cell type!")
      }
      if (verbose) {
        message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
      }
      D1 = D1[, valid.ct]
      M.S = M.S[valid.ct]
      Sigma = Sigma[, valid.ct]
      Est.prop.allgene = NULL
      Est.prop.weighted = NULL
      Weight.gene = NULL
      r.squared.full = NULL
      Var.prop = NULL
      for (i in 1:N.bulk) {
        if (sum(Yjg[, i] == 0) > 0) {
          D1.temp = D1[Yjg[, i] != 0, ]
          Yjg.temp = Yjg[Yjg[, i] != 0, i]
          Sigma.temp = Sigma[Yjg[, i] != 0, ]
          if (verbose) 
            message(paste(colnames(Yjg)[i], "has common genes", 
                          sum(Yjg[, i] != 0), "..."))
        } else {
          D1.temp = D1
          Yjg.temp = Yjg[, i]
          Sigma.temp = Sigma
          if (verbose) 
            message(paste(colnames(Yjg)[i], "has common genes", 
                          sum(Yjg[, i] != 0), "..."))
        }
        lm.D1.weighted = music.iter(Yjg.temp, D1.temp, M.S, 
                                    Sigma.temp, iter.max = iter.max, nu = nu, eps = eps, 
                                    centered = centered, normalize = normalize)
        Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
        Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
        weight.gene.temp = rep(NA, nrow(Yjg))
        weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
        Weight.gene = cbind(Weight.gene, weight.gene.temp)
        r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
        Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
      }
    }
    colnames(Est.prop.weighted) = colnames(D1)
    rownames(Est.prop.weighted) = colnames(Yjg)
    colnames(Est.prop.allgene) = colnames(D1)
    rownames(Est.prop.allgene) = colnames(Yjg)
    names(r.squared.full) = colnames(Yjg)
    colnames(Weight.gene) = colnames(Yjg)
    rownames(Weight.gene) = cm.gene
    colnames(Var.prop) = colnames(D1)
    rownames(Var.prop) = colnames(Yjg)
    return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene, 
                Weight.gene = Weight.gene, r.squared.full = r.squared.full, 
                Var.prop = Var.prop))
  }
  music_basis <- function (x, non.zero = TRUE, markers = NULL, clusters, samples, 
                           select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE) {
    if (!is.null(select.ct)) {
      x = x[, x@colData[, clusters] %in% select.ct]
    }
    if (non.zero) {
      x <- x[rowSums(counts(x)) > 0, ]
    }
    clusters <- as.character(colData(x)[, clusters])
    samples <- as.character(colData(x)[, samples])
    M.theta <- sapply(unique(clusters), function(ct) {
      my.rowMeans(sapply(unique(samples), function(sid) {
        y = counts(x)[, clusters %in% ct & samples %in% sid]
        if (is.null(dim(y))) {
          return(y/sum(y))
        } else {
          return(rowSums(y)/sum(y))
        }
      }), na.rm = TRUE)
    })
    if (verbose) {
      message("Creating Relative Abudance Matrix...")
    }
    if (ct.cov) {
      nGenes = nrow(x)
      n.ct = length(unique(clusters))
      nSubs = length(unique(samples))
      Theta <- sapply(unique(clusters), function(ct) {
        sapply(unique(samples), function(sid) {
          y = counts(x)[, clusters %in% ct & samples %in% 
                          sid]
          if (is.null(dim(y))) {
            return(y/sum(y))
          } else {
            return(rowSums(y)/sum(y))
          }
        })
      })
      if (!is.null(select.ct)) {
        m.ct = match(select.ct, colnames(Theta))
        Theta = Theta[, m.ct]
      }
      Sigma.ct = sapply(1:nGenes, function(g) {
        sigma.temp = Theta[nGenes * (0:(nSubs - 1)) + g, 
        ]
        Cov.temp = cov(sigma.temp)
        Cov.temp1 = cov(sigma.temp[rowSums(is.na(Theta[nGenes * 
                                                         (0:(nSubs - 1)) + 1, ])) == 0, ])
        Cov.temp[which(colSums(is.na(sigma.temp)) > 0), ] = Cov.temp1[which(colSums(is.na(sigma.temp)) > 
                                                                              0), ]
        Cov.temp[, which(colSums(is.na(sigma.temp)) > 0)] = Cov.temp1[, 
                                                                      which(colSums(is.na(sigma.temp)) > 0)]
        return(Cov.temp)
      })
      colnames(Sigma.ct) = rownames(x)
      if (!is.null(markers)) {
        ids <- intersect(unlist(markers), rownames(x))
        m.ids = match(ids, rownames(x))
        Sigma.ct <- Sigma.ct[, m.ids]
      }
      if (verbose) {
        message("Creating Covariance Matrix...")
      }
    } else {
      Sigma <- sapply(unique(clusters), function(ct) {
        apply(sapply(unique(samples), function(sid) {
          y = counts(x)[, clusters %in% ct & samples %in% 
                          sid]
          if (is.null(dim(y))) {
            return(y/sum(y))
          } else {
            return(rowSums(y)/sum(y))
          }
        }), 1, stats::var, na.rm = TRUE)
      })
      if (!is.null(select.ct)) {
        m.ct = match(select.ct, colnames(Sigma))
        Sigma = Sigma[, m.ct]
      }
      if (!is.null(markers)) {
        ids <- intersect(unlist(markers), rownames(x))
        m.ids = match(ids, rownames(x))
        Sigma <- Sigma[m.ids, ]
      }
      if (verbose) {
        message("Creating Variance Matrix...")
      }
    }
    S <- sapply(unique(clusters), function(ct) {
      my.rowMeans(sapply(unique(samples), function(sid) {
        y = counts(x)[, clusters %in% ct & samples %in% sid]
        if (is.null(dim(y))) {
          return(sum(y))
        } else {
          return(sum(y)/ncol(y))
        }
      }), na.rm = TRUE)
    })
    if (verbose) {
      message("Creating Library Size Matrix...")
    }
    S[S == 0] = NA
    M.S = colMeans(S, na.rm = TRUE)
    if (!is.null(cell_size)) {
      if (!is.data.frame(cell_size)) {
        stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
      } else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
        stop("Cell type names in cell_size must match clusters")
      } else if (any(is.na(as.numeric(cell_size[, 2])))) {
        stop("Cell sizes should all be numeric")
      }
      my_ms_names <- names(M.S)
      cell_size <- cell_size[my_ms_names %in% cell_size[, 1], 
      ]
      M.S <- cell_size[match(my_ms_names, cell_size[, 1]), 
      ]
      M.S <- M.S[, 2]
      names(M.S) <- my_ms_names
    }
    D <- t(t(M.theta) * M.S)
    if (!is.null(select.ct)) {
      m.ct = match(select.ct, colnames(D))
      D = D[, m.ct]
      S = S[, m.ct]
      M.S = M.S[m.ct]
      M.theta = M.theta[, m.ct]
    }
    if (!is.null(markers)) {
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      D <- D[m.ids, ]
      M.theta <- M.theta[m.ids, ]
    }
    if (ct.cov) {
      return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, 
                  Sigma.ct = Sigma.ct))
    } else {
      return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, 
                  Sigma = Sigma))
    }
  }
  music.iter <- function (Y, D, S, Sigma, iter.max = 1000, nu = 1e-04, eps = 0.01, 
                          centered = FALSE, normalize = FALSE) {
    if (length(S) != ncol(D)) {
      common.cell.type = intersect(colnames(D), names(S))
      if (length(common.cell.type) <= 1) {
        stop("Not enough cell types!")
      }
      D = D[, match(common.cell.type, colnames(D))]
      S = S[match(common.cell.type, names(S))]
    }
    if (ncol(Sigma) != ncol(D)) {
      common.cell.type = intersect(colnames(D), colnames(Sigma))
      if (length(common.cell.type) <= 1) {
        stop("Not enough cell type!")
      }
      D = D[, match(common.cell.type, colnames(D))]
      Sigma = Sigma[, match(common.cell.type, colnames(Sigma))]
      S = S[match(common.cell.type, names(S))]
    }
    k = ncol(D)
    common.gene = intersect(names(Y), rownames(D))
    if (length(common.gene) < 0.1 * min(length(Y), nrow(D))) {
      stop("Not enough common genes!")
    }
    Y = Y[match(common.gene, names(Y))]
    D = D[match(common.gene, rownames(D)), ]
    Sigma = Sigma[match(common.gene, rownames(Sigma)), ]
    X = D
    if (centered) {
      X = X - mean(X)
      Y = Y - mean(Y)
    }
    if (normalize) {
      X = X/stats::sd(as.vector(X))
      S = S * stats::sd(as.vector(X))
      Y = Y/stats::sd(Y)
    } else {
      Y = Y * 100
    }
    lm.D = music.basic(Y, X, S, Sigma, iter.max = iter.max, nu = nu, 
                       eps = eps)
    return(lm.D)
  }
  music.basic <- function (Y, X, S, Sigma, iter.max, nu, eps) {
    k = ncol(X)
    lm.D = nnls(X, Y)
    r = resid(lm.D)
    weight.gene = 1/(nu + r^2 + colSums((lm.D$x * S)^2 * t(Sigma)))
    Y.weight = Y * sqrt(weight.gene)
    D.weight = sweep(X, 1, sqrt(weight.gene), "*")
    lm.D.weight = nnls(D.weight, Y.weight)
    p.weight = lm.D.weight$x/sum(lm.D.weight$x)
    p.weight.iter = p.weight
    r = resid(lm.D.weight)
    for (iter in 1:iter.max) {
      weight.gene = 1/(nu + r^2 + colSums((lm.D.weight$x * 
                                             S)^2 * t(Sigma)))
      Y.weight = Y * sqrt(weight.gene)
      D.weight = X * as.matrix(sqrt(weight.gene))[, rep(1, 
                                                        k)]
      lm.D.weight = nnls(D.weight, Y.weight)
      p.weight.new = lm.D.weight$x/sum(lm.D.weight$x)
      r.new = resid(lm.D.weight)
      if (sum(abs(p.weight.new - p.weight)) < eps) {
        p.weight = p.weight.new
        r = r.new
        R.squared = 1 - stats::var(Y - X %*% as.matrix(lm.D.weight$x))/stats::var(Y)
        fitted = X %*% as.matrix(lm.D.weight$x)
        var.p = diag(solve(t(D.weight) %*% D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
        return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x, 
                    fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D), 
                    p.weight = p.weight, q.weight = lm.D.weight$x, 
                    fit.weight = fitted, resid.weight = Y - X %*% 
                      as.matrix(lm.D.weight$x), weight.gene = weight.gene, 
                    converge = paste0("Converge at ", iter), rsd = r, 
                    R.squared = R.squared, var.p = var.p))
      }
      p.weight = p.weight.new
      r = r.new
    }
    fitted = X %*% as.matrix(lm.D.weight$x)
    R.squared = 1 - stats::var(Y - X %*% as.matrix(lm.D.weight$x))/stats::var(Y)
    var.p = diag(solve(t(D.weight) %*% D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
    return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x, 
                fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D), p.weight = p.weight, 
                q.weight = lm.D.weight$x, fit.weight = fitted, resid.weight = Y - 
                  X %*% as.matrix(lm.D.weight$x), weight.gene = weight.gene, 
                converge = "Reach Maxiter", rsd = r, R.squared = R.squared, 
                var.p = var.p))
  }
  
  
}
{
  setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
  RNA_singlet <- readRDS("RNA_1.RNA_singlet.rds")
  bulk_data <- read.csv("antler-bulkrna.matrix.gene", sep = "\t",
                   header = TRUE, row.names = 1)
  common_genes <- intersect(rownames(bulk_data),
                            rownames(RNA_singlet))
  bulk_data <- bulk_data[common_genes,]
  RNA_singlet_2 <- subset(RNA_singlet, features = common_genes)
  
  # 假设 bulk_data 是一个矩阵，行为基因，列为样本
  bulk_sample_ids <- colnames(bulk_data)
  
  # 创建 phenoData（样本信息）
  bulk_pheno_data <- data.frame(sample_id = bulk_sample_ids)
  rownames(bulk_pheno_data) <- bulk_sample_ids
  bulk_pheno_data <- AnnotatedDataFrame(bulk_pheno_data)
  
  # 创建 ExpressionSet 对象
  bulk_eset <- ExpressionSet(assayData = as.matrix(bulk_data), phenoData = bulk_pheno_data)
  
  bulk.mtx <- Biobase::exprs(bulk_eset)
  
  sc_sce <- as.SingleCellExperiment(RNA_singlet_2)
  
  library(SingleCellExperiment)
  # 确保 SingleCellExperiment 包含细胞类型和样本信息
  colData(sc_sce)$cell_type <- as.character(RNA_singlet_2$Celltype_rename)
  colData(sc_sce)$sample_id <- as.character(RNA_singlet_2$orig.ident)
  
  results <- music_prop(
    bulk.mtx = bulk.mtx,       # bulk RNA-seq ExpressionSet
    sc.sce = sc_sce,           # 单细胞 RNA-seq ExpressionSet
    clusters = "cell_type",      # 细胞类型标签列
    samples = "sample_id"        # 样本标签列
  )
  saveRDS(results, "MuSiC_result.rds")
  
  temp <- results$Est.prop.weighted
  temp <- cor(t(temp))
  pheatmap::pheatmap(as.matrix(temp),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE)
}

#### BisqueRNA
library(BisqueRNA)
{
  setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
  RNA_singlet <- readRDS("RNA_1.RNA_singlet.rds")
  bulk_data <- read.csv("antler-bulkrna.matrix.gene", sep = "\t",
                        header = TRUE, row.names = 1)
  common_genes <- intersect(rownames(bulk_data),
                            rownames(RNA_singlet))
  bulk_data <- bulk_data[common_genes,]
  RNA_singlet_2 <- subset(RNA_singlet, features = common_genes)
  
  # 提取表达矩阵、细胞类型和样本信息
  expr_matrix <- GetAssayData(RNA_singlet_2, slot = "counts")  # 或者 "data" slot
  cell_types <- RNA_singlet_2$Celltype_rename
  sample_ids <- RNA_singlet_2$orig.ident  # 假设 "sample_id" 是样本信息列
  # 创建 phenoData（样本信息）
  pheno_data <- data.frame(cell_type = cell_types,
                           SubjectName = sample_ids)
  rownames(pheno_data) <- colnames(expr_matrix)
  pheno_data <- AnnotatedDataFrame(pheno_data)
  # 创建 ExpressionSet 对象
  sc_eset <- ExpressionSet(assayData = as.matrix(expr_matrix), phenoData = pheno_data)
  
  # 假设 bulk_data 是一个矩阵，行为基因，列为样本
  bulk_sample_ids <- colnames(bulk_data)
  # 创建 phenoData（样本信息）
  bulk_pheno_data <- data.frame(sample_id = bulk_sample_ids)
  rownames(bulk_pheno_data) <- bulk_sample_ids
  bulk_pheno_data <- AnnotatedDataFrame(bulk_pheno_data)
  # 创建 ExpressionSet 对象
  bulk_eset <- ExpressionSet(assayData = as.matrix(bulk_data),
                             phenoData = bulk_pheno_data)
  
  
  result <- ReferenceBasedDecomposition(bulk_eset,
                                        sc_eset,
                                        cell.types = "cell_type",
                                        use.overlap = FALSE)
  saveRDS(result[["bulk.props"]],
          "BisqueRNA_result.rds")
  temp <- result[["bulk.props"]]
  temp <- as.matrix(cor(temp))
  pheatmap::pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE)
}

##### Scaden
{
  setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
  RNA_singlet <- readRDS("RNA_1.RNA_singlet.rds")
  bulk_data <- read.csv("antler-bulkrna.matrix.gene", sep = "\t",
                        header = TRUE, row.names = 1)
  common_genes <- intersect(rownames(bulk_data),
                            rownames(RNA_singlet))
  bulk_data <- bulk_data[common_genes,]
  RNA_singlet_2 <- subset(RNA_singlet, features = common_genes)
  
  bulk_data <- bulk_data[rowSums(bulk_data) > 0,]
  data.table::fwrite(bulk_data,
                     file = "./scaden/bulk_data.txt",
                     sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  sc_exp <- RNA_singlet_2@assays$RNA@counts
  sc_exp <- t(sc_exp)
  sc_exp <- sc_exp[,rownames(bulk_data)]
  sc_exp <- as.matrix(sc_exp)
  data.table::fwrite(sc_exp,
                     file = "./scaden/sce_data.txt",
                     sep = "\t", quote = FALSE,
                     row.names = TRUE, col.names = TRUE)
  data.table::fwrite(data.frame(Celltype = as.character(RNA_singlet_2$Celltype_rename)),
                     file = "./scaden/sce_meta.txt",
                     sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE)
  
  temp <- read.csv("./scaden/scaden_predictions.txt",
                   sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE)
  temp <- as.matrix(cor(t(temp)))
  pheatmap::pheatmap(temp, cluster_rows = FALSE, cluster_cols = FALSE)
}


##### 结果可视化_V1
{
  RNA_Celltype_color <- c("#A4C9DD", "#2572A9", "#ADD487", "#399938",
                          "#F19695", "#D5231E", "#F5BB6F", "#EF7C1C",
                          "#C6B0D2", "#653B90")
  names(RNA_Celltype_color) <- c("AnMCs", "Proliferative_progenitor cells", "Progenitor cells",
                                 "Chondrocytes", "Hypertrophic chondrocytes", "Chondroclasts",
                                 "Mural cells", "Endothelial cells", "Monocytes_Macrophages",
                                 "Mast cells")
  Single_Percent <- as.data.frame.array(table(RNA_singlet$Celltype_rename,
                                              RNA_singlet$orig.ident))
  colnames(Single_Percent) <- unlist(lapply(colnames(Single_Percent), function(x){
    unlist(strsplit(x, ".", fixed = TRUE))[2]
  }))
  Single_Percent_2 <- data.frame(Celltype = rownames(Single_Percent),
                                 Single_Percent,
                                 check.rows = FALSE, check.names = FALSE)
  Single_Percent_2 <- reshape::melt.data.frame(Single_Percent_2,
                                               id.vars = "Celltype")
  pdf("./反卷积结果/1.单细胞_细胞类型比例barplot.pdf",
      width = 6.8, height = 4)
  ggplot(data = Single_Percent_2,
         aes(x = variable, y = value, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Percentage", fill = "Cell type") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  
  Single_Percent <- as.data.frame.array(table(RNA_singlet$Celltype_rename,
                                              RNA_singlet$orig.ident))
  colnames(Single_Percent) <- unlist(lapply(colnames(Single_Percent), function(x){
    unlist(strsplit(x, ".", fixed = TRUE))[2]
  }))
  colnames(Single_Percent) <- paste0(colnames(Single_Percent),
                                     "_sc")
  for (i in 1:ncol(Single_Percent)) {
    Single_Percent[,i] <- Single_Percent[,i] / sum(Single_Percent[,i])
  }
  pdf("./反卷积结果/1.单细胞样本_细胞类型比例heatmap.pdf",
      width = 5, height = 5)
  pheatmap::pheatmap(as.matrix(cor(Single_Percent)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
}

# BayesPrism
{
  BayesPrism <- readRDS("BayesPrism_result.rds")
  BayesPrism_2 <- data.frame(Celltype = rownames(BayesPrism),
                           BayesPrism,
                           check.rows = FALSE, check.names = FALSE)
  BayesPrism_2 <- reshape::melt.data.frame(BayesPrism_2,
                                           id.vars = "Celltype")
  pdf("./反卷积结果/2.BayesPrism_细胞类型比例barplot.pdf",
      width = 6.8, height = 4)
  ggplot(data = BayesPrism_2,
         aes(x = variable, y = value, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Percentage", fill = "Cell type") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  colnames(BayesPrism) <- paste0(colnames(BayesPrism), "_bk")
  pdf("./反卷积结果/2.BayesPrism_细胞类型比例heatmap.pdf",
      width = 6, height = 6)
  pheatmap::pheatmap(as.matrix(cor(BayesPrism)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
  
  pdf("./反卷积结果/2.BayesPrism_细胞类型比例heatmap.pdf",
      width = 6, height = 6)
  pheatmap::pheatmap(as.matrix(cor(BayesPrism)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
  
  celltypes <- rownames(BayesPrism)
  merged <- as.data.frame(cbind(BayesPrism,
                                Single_Percent[celltypes,]))
  samples <- c("CA1_bk", "CA2_bk", "CA3_bk", "CA1_sc", "CA2_sc", "CA3_sc",
               "MC1_bk", "MC2_bk", "MC3_bk", "MC1_sc", "MC2_sc", "MC3_sc",
               "PC1_bk", "PC2_bk", "PC3_bk", "PC1_sc", "PC2_sc",
               "RM1_bk", "RM2_bk", "RM3_bk", "RM1_sc", "RM2_sc",
               "TZ1_bk", "TZ2_bk", "TZ3_bk", "TZ1_sc", "TZ2_sc", "TZ3_sc")
  merged <- merged[,samples]
  pdf("./反卷积结果/2.BayesPrism_bulk_细胞类型比例heatmap.pdf",
      width = 10, height = 10)
  pheatmap::pheatmap(as.matrix(cor(merged)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
}

# MuSiC
{
  MuSiC <- readRDS("MuSiC_result.rds")
  MuSiC <- MuSiC[["Est.prop.weighted"]]
  MuSiC <- data.frame(t(MuSiC),
                      check.rows = FALSE, check.names = FALSE)
  MuSiC_2 <- data.frame(Celltype = rownames(MuSiC),
                        MuSiC,
                        check.rows = FALSE, check.names = FALSE)
  MuSiC_2 <- reshape::melt.data.frame(MuSiC_2,
                                      id.vars = "Celltype")
  pdf("./反卷积结果/3.MuSiC_细胞类型比例barplot.pdf",
      width = 6.8, height = 4)
  ggplot(data = MuSiC_2,
         aes(x = variable, y = value, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Percentage", fill = "Cell type") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  colnames(MuSiC) <- paste0(colnames(MuSiC), "_bk")
  pdf("./反卷积结果/3.MuSiC_细胞类型比例heatmap.pdf",
      width = 6, height = 6)
  pheatmap::pheatmap(as.matrix(cor(MuSiC)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
  
  celltypes <- rownames(MuSiC)
  merged <- as.data.frame(cbind(MuSiC,
                                Single_Percent[celltypes,]))
  samples <- c("CA1_bk", "CA2_bk", "CA3_bk", "CA1_sc", "CA2_sc", "CA3_sc",
               "MC1_bk", "MC2_bk", "MC3_bk", "MC1_sc", "MC2_sc", "MC3_sc",
               "PC1_bk", "PC2_bk", "PC3_bk", "PC1_sc", "PC2_sc",
               "RM1_bk", "RM2_bk", "RM3_bk", "RM1_sc", "RM2_sc",
               "TZ1_bk", "TZ2_bk", "TZ3_bk", "TZ1_sc", "TZ2_sc", "TZ3_sc")
  merged <- merged[,samples]
  pdf("./反卷积结果/3.MuSiC_bulk_细胞类型比例heatmap.pdf",
      width = 10, height = 10)
  pheatmap::pheatmap(as.matrix(cor(merged)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
}

# BisqueRNA
{
  BisqueRNA <- readRDS("BisqueRNA_result.rds")
  BisqueRNA_2 <- data.frame(Celltype = rownames(BisqueRNA),
                            BisqueRNA,
                        check.rows = FALSE, check.names = FALSE)
  BisqueRNA_2 <- reshape::melt.data.frame(BisqueRNA_2,
                                      id.vars = "Celltype")
  pdf("./反卷积结果/4.BisqueRNA_细胞类型比例barplot.pdf",
      width = 6.8, height = 4)
  ggplot(data = BisqueRNA_2,
         aes(x = variable, y = value, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Percentage", fill = "Cell type") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  colnames(BisqueRNA) <- paste0(colnames(BisqueRNA), "_bk")
  pdf("./反卷积结果/4.BisqueRNA_细胞类型比例heatmap.pdf",
      width = 6, height = 6)
  pheatmap::pheatmap(as.matrix(cor(BisqueRNA)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
  
  celltypes <- rownames(BisqueRNA)
  merged <- as.data.frame(cbind(BisqueRNA,
                                Single_Percent[celltypes,]))
  samples <- c("CA1_bk", "CA2_bk", "CA3_bk", "CA1_sc", "CA2_sc", "CA3_sc",
               "MC1_bk", "MC2_bk", "MC3_bk", "MC1_sc", "MC2_sc", "MC3_sc",
               "PC1_bk", "PC2_bk", "PC3_bk", "PC1_sc", "PC2_sc",
               "RM1_bk", "RM2_bk", "RM3_bk", "RM1_sc", "RM2_sc",
               "TZ1_bk", "TZ2_bk", "TZ3_bk", "TZ1_sc", "TZ2_sc", "TZ3_sc")
  merged <- merged[,samples]
  pdf("./反卷积结果/4.BisqueRNA_bulk_细胞类型比例heatmap.pdf",
      width = 10, height = 10)
  pheatmap::pheatmap(as.matrix(cor(merged)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
}

# Scaden
{
  Scaden <- read.csv("./scaden/scaden_predictions.txt",
                     sep = "\t", header = TRUE, row.names = 1,
                     check.names = FALSE)
  Scaden <- data.frame(t(Scaden),
                       check.rows = FALSE, check.names = FALSE)
  Scaden_2 <- data.frame(Celltype = rownames(Scaden),
                         Scaden,
                         check.rows = FALSE, check.names = FALSE)
  Scaden_2 <- reshape::melt.data.frame(Scaden_2,
                                       id.vars = "Celltype")
  pdf("./反卷积结果/5.Scaden_细胞类型比例barplot.pdf",
      width = 6.8, height = 4)
  ggplot(data = Scaden_2,
         aes(x = variable, y = value, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Percentage", fill = "Cell type") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  colnames(Scaden) <- paste0(colnames(Scaden), "_bk")
  pdf("./反卷积结果/5.Scaden_细胞类型比例heatmap.pdf",
      width = 6, height = 6)
  pheatmap::pheatmap(as.matrix(cor(Scaden)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
  
  celltypes <- rownames(Scaden)
  merged <- as.data.frame(cbind(Scaden,
                                Single_Percent[celltypes,]))
  samples <- c("CA1_bk", "CA2_bk", "CA3_bk", "CA1_sc", "CA2_sc", "CA3_sc",
               "MC1_bk", "MC2_bk", "MC3_bk", "MC1_sc", "MC2_sc", "MC3_sc",
               "PC1_bk", "PC2_bk", "PC3_bk", "PC1_sc", "PC2_sc",
               "RM1_bk", "RM2_bk", "RM3_bk", "RM1_sc", "RM2_sc",
               "TZ1_bk", "TZ2_bk", "TZ3_bk", "TZ1_sc", "TZ2_sc", "TZ3_sc")
  merged <- merged[,samples]
  pdf("./反卷积结果/5.Scaden_bulk_细胞类型比例heatmap.pdf",
      width = 10, height = 10)
  pheatmap::pheatmap(as.matrix(cor(merged)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
}


#### BisqueRNA_V2
library(BisqueRNA)
{
  setwd("/media/heshidian/RAID5_42TB/0.MyFiles/1.BGI/6.Cartialgenous")
  RNA_singlet <- readRDS("RNA_1.RNA_singlet.rds")
  bulk_data <- read.csv("antler-bulkrna.matrix.gene", sep = "\t",
                        header = TRUE, row.names = 1)
  # RM
  {
    RM_bulk_data <- bulk_data[,c("RM1", "RM2", "RM3")]
    RM_bulk_data <- RM_bulk_data[rowSums(RM_bulk_data) > 0,]
    RM_sc <- subset(RNA_singlet, subset = sample == "1.RM")
    common_genes <- intersect(rownames(RM_bulk_data),
                              rownames(RM_sc))
    RM_bulk_data <- RM_bulk_data[common_genes,]
    RM_sc <- subset(RM_sc, features = common_genes)

    expr_matrix <- GetAssayData(RM_sc, slot = "counts") 
    cell_types <- RM_sc$Celltype_rename
    sample_ids <- RM_sc$orig.ident  
    pheno_data <- data.frame(cell_type = cell_types,
                             SubjectName = sample_ids)
    rownames(pheno_data) <- colnames(expr_matrix)
    pheno_data <- AnnotatedDataFrame(pheno_data)
    sc_eset <- ExpressionSet(assayData = as.matrix(expr_matrix),
                             phenoData = pheno_data)
    
    bulk_sample_ids <- colnames(RM_bulk_data)
    bulk_pheno_data <- data.frame(sample_id = bulk_sample_ids)
    rownames(bulk_pheno_data) <- bulk_sample_ids
    bulk_pheno_data <- AnnotatedDataFrame(bulk_pheno_data)
    # 创建 ExpressionSet 对象
    bulk_eset <- ExpressionSet(assayData = as.matrix(RM_bulk_data),
                               phenoData = bulk_pheno_data)
    
    RM_result <- ReferenceBasedDecomposition(bulk_eset,
                                          sc_eset,
                                          cell.types = "cell_type",
                                          use.overlap = FALSE)
    RM_result <- RM_result[["bulk.props"]]
    
  }
  # PC
  {
    PC_bulk_data <- bulk_data[,c("PC1", "PC2", "PC3")]
    PC_bulk_data <- PC_bulk_data[rowSums(PC_bulk_data) > 0,]
    PC_sc <- subset(RNA_singlet, subset = sample == "2.PC")
    common_genes <- intersect(rownames(PC_bulk_data),
                              rownames(PC_sc))
    PC_bulk_data <- PC_bulk_data[common_genes,]
    PC_sc <- subset(PC_sc, features = common_genes)
    
    expr_matrix <- GetAssayData(PC_sc, slot = "counts") 
    cell_types <- PC_sc$Celltype_rename
    sample_ids <- PC_sc$orig.ident  
    pheno_data <- data.frame(cell_type = cell_types,
                             SubjectName = sample_ids)
    rownames(pheno_data) <- colnames(expr_matrix)
    pheno_data <- AnnotatedDataFrame(pheno_data)
    sc_eset <- ExpressionSet(assayData = as.matrix(expr_matrix),
                             phenoData = pheno_data)
    
    bulk_sample_ids <- colnames(PC_bulk_data)
    bulk_pheno_data <- data.frame(sample_id = bulk_sample_ids)
    rownames(bulk_pheno_data) <- bulk_sample_ids
    bulk_pheno_data <- AnnotatedDataFrame(bulk_pheno_data)
    # 创建 ExpressionSet 对象
    bulk_eset <- ExpressionSet(assayData = as.matrix(PC_bulk_data),
                               phenoData = bulk_pheno_data)
    
    PC_result <- ReferenceBasedDecomposition(bulk_eset,
                                             sc_eset,
                                             cell.types = "cell_type",
                                             use.overlap = FALSE)
    PC_result <- PC_result[["bulk.props"]]
    
  }
  # TZ
  {
    TZ_bulk_data <- bulk_data[,c("TZ1", "TZ2", "TZ3")]
    TZ_bulk_data <- TZ_bulk_data[rowSums(TZ_bulk_data) > 0,]
    TZ_sc <- subset(RNA_singlet, subset = sample == "3.TZ")
    common_genes <- intersect(rownames(TZ_bulk_data),
                              rownames(TZ_sc))
    TZ_bulk_data <- TZ_bulk_data[common_genes,]
    TZ_sc <- subset(TZ_sc, features = common_genes)
    
    expr_matrix <- GetAssayData(TZ_sc, slot = "counts") 
    cell_types <- TZ_sc$Celltype_rename
    sample_ids <- TZ_sc$orig.ident  
    pheno_data <- data.frame(cell_type = cell_types,
                             SubjectName = sample_ids)
    rownames(pheno_data) <- colnames(expr_matrix)
    pheno_data <- AnnotatedDataFrame(pheno_data)
    sc_eset <- ExpressionSet(assayData = as.matrix(expr_matrix),
                             phenoData = pheno_data)
    
    bulk_sample_ids <- colnames(TZ_bulk_data)
    bulk_pheno_data <- data.frame(sample_id = bulk_sample_ids)
    rownames(bulk_pheno_data) <- bulk_sample_ids
    bulk_pheno_data <- AnnotatedDataFrame(bulk_pheno_data)
    # 创建 ExpressionSet 对象
    bulk_eset <- ExpressionSet(assayData = as.matrix(TZ_bulk_data),
                               phenoData = bulk_pheno_data)
    
    TZ_result <- ReferenceBasedDecomposition(bulk_eset,
                                             sc_eset,
                                             cell.types = "cell_type",
                                             use.overlap = FALSE)
    TZ_result <- TZ_result[["bulk.props"]]
    
  }
  # CA
  {
    CA_bulk_data <- bulk_data[,c("CA1", "CA2", "CA3")]
    CA_bulk_data <- CA_bulk_data[rowSums(CA_bulk_data) > 0,]
    CA_sc <- subset(RNA_singlet, subset = sample == "4.CA")
    common_genes <- intersect(rownames(CA_bulk_data),
                              rownames(CA_sc))
    CA_bulk_data <- CA_bulk_data[common_genes,]
    CA_sc <- subset(CA_sc, features = common_genes)
    
    expr_matrix <- GetAssayData(CA_sc, slot = "counts") 
    cell_types <- CA_sc$Celltype_rename
    sample_ids <- CA_sc$orig.ident  
    pheno_data <- data.frame(cell_type = cell_types,
                             SubjectName = sample_ids)
    rownames(pheno_data) <- colnames(expr_matrix)
    pheno_data <- AnnotatedDataFrame(pheno_data)
    sc_eset <- ExpressionSet(assayData = as.matrix(expr_matrix),
                             phenoData = pheno_data)
    
    bulk_sample_ids <- colnames(CA_bulk_data)
    bulk_pheno_data <- data.frame(sample_id = bulk_sample_ids)
    rownames(bulk_pheno_data) <- bulk_sample_ids
    bulk_pheno_data <- AnnotatedDataFrame(bulk_pheno_data)
    # 创建 ExpressionSet 对象
    bulk_eset <- ExpressionSet(assayData = as.matrix(CA_bulk_data),
                               phenoData = bulk_pheno_data)
    
    CA_result <- ReferenceBasedDecomposition(bulk_eset,
                                             sc_eset,
                                             cell.types = "cell_type",
                                             use.overlap = FALSE)
    CA_result <- CA_result[["bulk.props"]]
    
  }
  # MC
  {
    MC_bulk_data <- bulk_data[,c("MC1", "MC2", "MC3")]
    MC_bulk_data <- MC_bulk_data[rowSums(MC_bulk_data) > 0,]
    MC_sc <- subset(RNA_singlet, subset = sample == "5.MC")
    common_genes <- intersect(rownames(MC_bulk_data),
                              rownames(MC_sc))
    MC_bulk_data <- MC_bulk_data[common_genes,]
    MC_sc <- subset(MC_sc, features = common_genes)
    
    expr_matrix <- GetAssayData(MC_sc, slot = "counts") 
    cell_types <- MC_sc$Celltype_rename
    sample_ids <- MC_sc$orig.ident  
    pheno_data <- data.frame(cell_type = cell_types,
                             SubjectName = sample_ids)
    rownames(pheno_data) <- colnames(expr_matrix)
    pheno_data <- AnnotatedDataFrame(pheno_data)
    sc_eset <- ExpressionSet(assayData = as.matrix(expr_matrix),
                             phenoData = pheno_data)
    
    bulk_sample_ids <- colnames(MC_bulk_data)
    bulk_pheno_data <- data.frame(sample_id = bulk_sample_ids)
    rownames(bulk_pheno_data) <- bulk_sample_ids
    bulk_pheno_data <- AnnotatedDataFrame(bulk_pheno_data)
    # 创建 ExpressionSet 对象
    bulk_eset <- ExpressionSet(assayData = as.matrix(MC_bulk_data),
                               phenoData = bulk_pheno_data)
    
    MC_result <- ReferenceBasedDecomposition(bulk_eset,
                                             sc_eset,
                                             cell.types = "cell_type",
                                             use.overlap = FALSE)
    MC_result <- MC_result[["bulk.props"]]
    
  }
  
  View(RM_result)
  celltypes <- as.character(unique(RNA_singlet$Celltype_rename))
  
  RM_result <- as.data.frame(RM_result)[celltypes,]
  rownames(RM_result) <- celltypes
  PC_result <- as.data.frame(PC_result)[celltypes,]
  rownames(PC_result) <- celltypes
  TZ_result <- as.data.frame(TZ_result)[celltypes,]
  rownames(TZ_result) <- celltypes
  CA_result <- as.data.frame(CA_result)[celltypes,]
  rownames(CA_result) <- celltypes
  MC_result <- as.data.frame(MC_result)[celltypes,]
  rownames(MC_result) <- celltypes
  temp <- as.data.frame(cbind(RM_result,
                      PC_result,
                      TZ_result,
                      CA_result,
                      MC_result))
  temp <- as.matrix(temp)
  temp[is.na(temp)] <- 0
  
  BisqueRNA_V2 <- data.frame(Celltype = rownames(temp),
                            temp,
                            check.rows = FALSE, check.names = FALSE)
  BisqueRNA_V2 <- reshape::melt.data.frame(BisqueRNA_V2,
                                          id.vars = "Celltype")
  pdf("./反卷积结果/V2_BisqueRNA_细胞类型比例barplot.pdf",
      width = 6.8, height = 4)
  ggplot(data = BisqueRNA_V2,
         aes(x = variable, y = value, fill = Celltype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = RNA_Celltype_color) +
    theme_bw() +
    labs(x = "", y = "Percentage", fill = "Cell type") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10, color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  dev.off()
  
  colnames(temp) <- paste0(colnames(temp), "_bk")
  pdf("./反卷积结果/V2_BisqueRNA_细胞类型比例heatmap.pdf",
      width = 6, height = 6)
  pheatmap::pheatmap(as.matrix(cor(temp)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
  
  Single_Percent <- as.data.frame.array(table(RNA_singlet$Celltype_rename,
                                              RNA_singlet$orig.ident))
  colnames(Single_Percent) <- unlist(lapply(colnames(Single_Percent), function(x){
    unlist(strsplit(x, ".", fixed = TRUE))[2]
  }))
  colnames(Single_Percent) <- paste0(colnames(Single_Percent),
                                     "_sc")
  for (i in 1:ncol(Single_Percent)) {
    Single_Percent[,i] <- Single_Percent[,i] / sum(Single_Percent[,i])
  }
  celltypes <- rownames(temp)
  merged <- as.data.frame(cbind(temp,
                                Single_Percent[celltypes,]))
  samples <- c("CA1_bk", "CA2_bk", "CA3_bk", "CA1_sc", "CA2_sc", "CA3_sc",
               "MC1_bk", "MC2_bk", "MC3_bk", "MC1_sc", "MC2_sc", "MC3_sc",
               "PC1_bk", "PC2_bk", "PC3_bk", "PC1_sc", "PC2_sc",
               "RM1_bk", "RM2_bk", "RM3_bk", "RM1_sc", "RM2_sc",
               "TZ1_bk", "TZ2_bk", "TZ3_bk", "TZ1_sc", "TZ2_sc", "TZ3_sc")
  merged <- merged[,samples]
  pdf("./反卷积结果/V2_BisqueRNA_bulk_细胞类型比例heatmap.pdf",
      width = 10, height = 10)
  pheatmap::pheatmap(as.matrix(cor(merged)),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     display_numbers = TRUE,
                     number_color = "black",
                     fontsize_number = 6,
                     cellwidth = 20, cellheight = 20)
  dev.off()
  
  
}

