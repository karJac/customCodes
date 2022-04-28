### Custom functions by KJ ###

convertHumanGeneList <- function(x) {
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}



############################################################################


notUniqueGenes <- function(object) {
  rownames(object)[!(isUnique(rownames(object)))]
}



############################################################################


sumNotUniqueGenes <- function(macExp) {
  genList <- notUniqueGenes(macExp)
  genList <- unique(genList)
  n <- 0
  myList <- list()
  for (i in genList) {
    n <- n + 1
    myList[[n]] <- colSums(macExp[rownames(macExp) == i, ])
    names(myList)[n] <- i
  }
  macExp <- macExp[!(rownames(macExp) %in% genList), ]
  tmp <- do.call(rbind, myList)
  rownames(tmp) <- names(myList)
  macExp <- rbind(macExp, tmp)
  return(macExp)
}



############################################################################


find.pK <- function(sweep.stats) # suppres creation of graphics
{
  "%ni%" <- Negate("%in%")
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    bc.mvn <- as.data.frame(matrix(0L,
      nrow = length(unique(sweep.stats$pK)),
      ncol = 5
    ))
    colnames(bc.mvn) <- c(
      "ParamID", "pK", "MeanBC", "VarBC",
      "BCmetric"
    )
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"]) / (sd(sweep.stats[
        ind,
        "BCreal"
      ])^2)
    }
    par(mar = rep(1, 4))
    return(bc.mvn)
  }
  if ("AUC" %in% colnames(sweep.stats) == TRUE) {
    bc.mvn <- as.data.frame(matrix(0L,
      nrow = length(unique(sweep.stats$pK)),
      ncol = 6
    ))
    colnames(bc.mvn) <- c(
      "ParamID", "pK", "MeanAUC", "MeanBC",
      "VarBC", "BCmetric"
    )
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"]) / (sd(sweep.stats[
        ind,
        "BCreal"
      ])^2)
    }
    par(mar = rep(1, 4))
    x <- plot(
      x = bc.mvn$ParamID, y = bc.mvn$MeanAUC, pch = 18,
      col = "black", cex = 0.75, xlab = NA, ylab = NA
    )
    x <- lines(
      x = bc.mvn$ParamID, y = bc.mvn$MeanAUC, col = "black",
      lty = 2
    )
    par(new = TRUE)
    return(bc.mvn)
  }
}



############################################################################


library(viridis)
umapp <- function(macExp, i) {
  theme_feature <- theme_classic() +
    theme(
      axis.text = element_text(24), axis.ticks = element_blank(),
      axis.title = element_text(size = 20),
      legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 20),
      plot.title = element_text(face = "italic", size = 60),
      legend.key.size = unit(1, "cm")
    )

  plot <- ggplot(macExp, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
    geom_jitter(size = 0.5, alpha = 0.5, aes_string(color = as.name(i))) + # czemu okregi a nie koła?
    labs(title = i, color = "log2(RPKM+1)") +
    scale_colour_gradientn(colours = rev(viridis(50))) + # rev(inferno(50))
    # breaks=c(min,max))+
    # , labels=c("min", "max"))+
    xlab("wnnUMAP_1") +
    ylab("wnnUMAP_2") +
    theme_feature
  return(plot)
}



############################################################################


testmac <- function(macExp, j) {
  for (i in j) {
    if (identical(colnames(macExp)[which(grepl(j, colnames(macExp)))], character(0))) {
      print(paste(i, " ", "is not in genes pool"))
    } else {
      return(colnames(macExp)[which(grepl(j, colnames(macExp)))])
    }
  }
}



############################################################################

testmacrev <- function(macExp, j) {
  for (i in j) {
    if (identical(rownames(macExp)[which(grepl(j, rownames(macExp)))], character(0))) {
      print(paste(i, " ", "is not in genes pool"))
    } else {
      return(rownames(macExp)[which(grepl(j, rownames(macExp)))])
    }
  }
}



############################################################################

testseu <- function(seuobj, j) {
  geneList <- list()
  n <- 0
  for (i in j) {
    n <- n + 1
    if (identical(rownames(seuobj)[which(grepl(i, rownames(seuobj)))], character(0))) {
      print(paste(i, " ", "is not in genes pool"))
    } else {
      print(rownames(seuobj)[which(grepl(i, rownames(seuobj)))])
      geneList[[n]] <- rownames(seuobj)[which(grepl(i, rownames(seuobj)))]
    }
  }
  return(unlist(geneList))
}



############################################################################

head2 <- function(i) {
  i[1:5, 1:5]
}



############################################################################

perc99 <- function(macExp) {
  for (i in colnames(macExp[, -c(1:4)])) {
    q99 <- quantile(macExp[which(macExp[, i] > 0), i], probs = 0.99)
    macExp[(which(macExp[, i] > q99)), i] <- q99
  }
  return(macExp)
}



############################################################################

greyClusters <- function(seuObj, idents, redu = "wnn.umap") {
  cellsInCluster <- WhichCells(seuObj, idents = idents)
  DimPlot(seuObj,
    reduction = redu,
    label = F, group.by = "ident",
    cells.highlight = list(cellsInCluster), cols.highlight = c("darkred"),
    cols = "grey", pt.size = 0.3
  ) + ggtitle(idents) + theme(legend.position = "none")
}




############################################################################

percTab <- function(alldata) {
  table(Idents(alldata), alldata$labels) %>%
    prop.table(margin = 2) %>%
    round(digits = 4) * 100
}




############################################################################

mybarplot <- function(height, x = "Count", color = "p.adjust", showCategory = 8,
                      font.size = 12, title = "", label_format = 30, ...) { # barplot working on enrichKEGG and enrichPathway result objects (default stopped to working xd)
  object <- height
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
  }
  df <- fortify(object@result[1:showCategory, ], by = x)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  if (colorBy %in% colnames(df)) {
    p <- ggplot(df, aes_string(
      x = x, y = "Description",
      fill = colorBy
    )) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue", name = color, guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(
      x = x, y = "Description",
      fill = "Description"
    )) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  p + geom_col() + ggtitle(title) + xlab(NULL) + ylab(NULL)
}




############################################################################


RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames
  } else {
    "Unequal gene sets: nrow(RNA) != nrow(newnames)"
  }
  obj@assays$RNA <- RNA
  return(obj)
}




############################################################################

h2m <- function(i) { # human-mouse
  lapply(i, function(j) {
    tmp <- unlist(strsplit(j, ""))
    tmp <- c(tmp[[1]], tolower(tmp[-1]))
    return(paste(tmp, collapse = ""))
  })
}
# as.character(expression(CD62L, CCR7, CXCR4))


#############################################################################







############################################################################


fp <- function(object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
                 c("lightgrey", "#ff0000", "#00ff00")
               } else {
                 c("#FDE725FF", "#440154FF")
               }, pt.size=NULL, pt.s = 1, order = FALSE, min.cutoff = NA, max.cutoff = NA,
               reduction = NULL, split.by = NULL, keep.scale = "feature",
               shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5,
               label = FALSE, label.size = 7, label.color = "midnightblue", label.bg.color="white", bg.r = 0.1, repel = TRUE, max.iter=1,autosize=FALSE,
               alpha=0.7, ncol = NULL, coord.fixed = FALSE, by.col = TRUE, sort.cell = NULL,
               interactive = FALSE, combine = TRUE, raster = NULL, split.ncol=NULL, split.nrow=NULL) {
  library(patchwork)
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ",
      "parameter instead for equivalent functionality.",
      call. = FALSE, immediate. = TRUE
    )
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(IFeaturePlot(
      object = object, feature = features[1],
      dims = dims, reduction = reduction, slot = slot
    ))
  }
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c(
    "feature",
    "all"
  ))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  no.right <- theme(
    axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_text(
      face = "bold.italic",
      size = 14, margin = margin(r = 7)
    )
  )
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
      `0` = {
        warning("No colors provided, using default colors",
          call. = FALSE, immediate. = TRUE
        )
        default.colors
      },
      `1` = {
        warning("Only one color provided, assuming specified is double-negative and augmenting with default colors",
          call. = FALSE, immediate. = TRUE
        )
        c(cols, default.colors[2:3])
      },
      `2` = {
        warning("Only two colors provided, assuming specified are for features and agumenting with '",
          default.colors[1], "' for double-negatives",
          call. = FALSE, immediate. = TRUE
        )
        c(default.colors[1], cols)
      },
      `3` = cols,
      {
        warning("More than three colors provided, using only first three",
          call. = FALSE, immediate. = TRUE
        )
        cols[1:3]
      }
    )
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(
    dims, "ident",
    features
  ), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ", paste(features,
      collapse = ", "
    ), " in slot ", slot, call. = FALSE)
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[
      ,
      feature
    ]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[
      ,
      feature
    ]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(
    features, min.cutoff,
    max.cutoff
  ), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index -
        3], data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index -
        3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      } else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(EXPR = split.by,
      ident = Idents(object = object)[cells,
        drop = TRUE
      ],
      object[[split.by, drop = TRUE]][cells,
        drop = TRUE
      ]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend,
    yes = 4, no = length(x = features) * length(x = levels(x = data$split))
  ))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[
    ,
    dims[1]
  ])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[
    ,
    dims[2]
  ])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(
      two.colors = cols[2:3], col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(
      color.matrix[, 1], color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, ,
      drop = FALSE
    ]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[
        ,
        features
      ]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ",
          paste(no.expression, collapse = ", "),
          call. = FALSE
        )
      }
      data.plot <- cbind(
        data.plot[, c(dims, "ident")],
        BlendExpression(data = data.plot[, features[1:2]])
      )
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[
          ,
          feature
        ])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(
        dims, "ident", feature,
        shape.by
      )]
      plot <- mySingleDimPlot(
        data = data.single, dims = dims, alpha=alpha,
        col.by = feature, order = order, pt.size = pt.size,
        cols = cols.use, shape.by = shape.by, label = FALSE,
        raster = raster, my.pt.size = pt.s, autosize = autosize,
      ) + scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) + theme_cowplot() +
        CenterTitle()
      if (label) {
        plot <- fplabelClustersWithOutilnedText(plot = plot, id = "ident", repel = repel, size = label.size, color = label.color, bg.color=label.bg.color, bg.r=bg.r,max.iter=max.iter)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(
          fill = NA,
          colour = "black"
        ))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(
            sec.axis = dup_axis(name = ident),
            limits = ylims
          ) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning(
              "All cells have the same value (",
              unique.feature.exp, ") of ", feature, "."
            )
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            } else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(
            colors = cols.grad,
            guide = "colorbar"
          ))
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" &&
        !blend) {
        max.feature.value <- max(data.single[, feature])
        min.feature.value <- min(data.single[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(
          colors = rev(viridis(50)),
          limits = c(min.feature.value, max.feature.value)
        ))
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(blend.legend + scale_y_continuous(
          sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
            1, yes = levels(x = data$split)[ii], no = "")),
          expand = c(0, 0)
        ) + labs(
          x = features[1], y = features[2],
          title = if (ii == 1) {
            paste("Color threshold:", blend.threshold)
          } else {
            NULL
          }
        ) + no.right), after = 4 * ii - 1
      ))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol, 
                 no = length(x = features))
  legend <- if (blend) {
    "none"
  } else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() +
          ggtitle("") + scale_y_continuous(
            sec.axis = dup_axis(name = ""),
            limits = ylims
          ) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) +
        1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
          scale_y_continuous(
            sec.axis = dup_axis(name = features[[idx]]),
            limits = ylims
          ) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) ==
        1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
           if (is.null(split.ncol)){
                nrow <- split.by %iff% length(x = levels(x = data$split))
           }
           else {
               nrow <- split.ncol
               ncol <- split.nrow
           }
      }
      plots <- plots[c(do.call(what = rbind, args = split(
        x = 1:length(x = plots),
        f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
      )))]
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
        length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" &&
      !blend) {
      max.feature.value <- max(data.plot[, features])
      min.feature.value <- min(data.plot[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(
        colors = rev(viridis(50)),
        limits = c(min.feature.value, max.feature.value)
      ))
    }
  }
  print("lol")
  return(plots)
}









############################################################################

featurePlotMultiple <- function(object.list, features.to.plot,
                                feature.labels = NULL,
                                filename,
                                assay = NULL,
                                ncol = 2,
                                minCutOff = NA,
                                maxCutOff = NA,
                                colorPalette = brewer.pal(11, "RdYlBu")[11:1]) {
  if (is.null(feature.labels)) {
    feature.labels <- features.to.plot
  }
  if (!is.null(assay)) {
    for (i in 1:length(object.list)) {
      DefaultAssay(object.list[[i]]) <- assay
    }
  }

  p.list <- lapply(object.list, function(object.to.plot) {
    lapply(seq_along(features.to.plot), function(i) {
      if (features.to.plot[i] %in% rownames(object.to.plot)) {
        FeaturePlot(object.to.plot,
          features = features.to.plot[i],
          pt.size = 0.9, min.cutoff = minCutOff, max.cutoff = maxCutOff
        ) +
          scale_colour_gradientn(colours = colorPalette) +
          ggtitle(feature.labels[i])
      } else {
        ggplot() +
          theme_void()
      }
    })
  })

  png(file = filename, width = 3840 / 2, height = 2400 / 2)
  print(CombinePlots(plots = unlist(p.list, recursive = F), ncol = ncol))
  dev.off()
}


splitDataFrame <- function(df, listToSplitBy) {
  n <- 0
  marklist <- list()
  for (i in listToSplitBy) {
    n <- n + 1
    marklist[[n]] <- markers[df$cluster == i, ]
  }
  return(marklist)
}









############################################################################


fplabelClustersWithOutilnedText <- function(plot, id, clusters = NULL, labels = NULL, split.by = NULL, repel = TRUE, box = FALSE, geom = "GeomPoint", position = "median", color = "white", bg.color = "midnightblue", bg.r=0.1,
                                            max.iter=1,
                                          ...) {
     xynames <- unlist(
          x = Seurat:::GetXYAesthetics(plot = plot, geom = geom),
          use.names = TRUE
     )
     if (!id %in% colnames(x = plot$data)) {
          stop("Cannot find variable ", id, " in plotting data")
     }
     if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
          warning("Cannot find splitting variable ", id, " in plotting data")
          split.by <- NULL
     }
     data <- plot$data[, c(xynames, id, split.by)]
     possible.clusters <- as.character(x = na.omit(object = unique(x = data[
          ,
          id
     ])))
     groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[
          ,
          id
     ])))
     if (any(!groups %in% possible.clusters)) {
          stop("The following clusters were not found: ", paste(groups[!groups %in%
                                                                            possible.clusters], collapse = ","))
     }
     pb <- ggplot_build(plot = plot)
     if (geom == "GeomSpatial") {
          xrange.save <- layer_scales(plot = plot)$x$range$range
          yrange.save <- layer_scales(plot = plot)$y$range$range
          data[, xynames["y"]] <- max(data[, xynames["y"]]) - data[
               ,
               xynames["y"]
          ] + min(data[, xynames["y"]])
          if (!pb$plot$plot_env$crop) {
               y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) -
                    pb$layout$panel_params[[1]]$y.range
               data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
          }
     }
     data <- cbind(data, color = pb$data[[1]][[1]])
     labels.loc <- lapply(X = groups, FUN = function(group) {
          data.use <- data[data[, id] == group, , drop = FALSE]
          data.medians <- if (!is.null(x = split.by)) {
               do.call(what = "rbind", args = lapply(X = unique(x = data.use[
                    ,
                    split.by
               ]), FUN = function(split) {
                    medians <- apply(
                         X = data.use[data.use[, split.by] ==
                                           split, xynames, drop = FALSE], MARGIN = 2,
                         FUN = median, na.rm = TRUE
                    )
                    medians <- as.data.frame(x = t(x = medians))
                    medians[, split.by] <- split
                    return(medians)
               }))
          } else {
               as.data.frame(x = t(x = apply(X = data.use[, xynames,
                                                          drop = FALSE
               ], MARGIN = 2, FUN = median, na.rm = TRUE)))
          }
          data.medians[, id] <- group
          data.medians$color <- data.use$color[1]
          return(data.medians)
     })
     if (position == "nearest") {
          labels.loc <- lapply(X = labels.loc, FUN = function(x) {
               group.data <- data[as.character(x = data[, id]) ==
                                       as.character(x[3]), ]
               nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(
                    1,
                    2
               )]), k = 1)$nn.idx
               x[1:2] <- group.data[nearest.point, 1:2]
               return(x)
          })
     }
     labels.loc <- do.call(what = "rbind", args = labels.loc)
     labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[
          ,
          id
     ]))
     labels <- labels %||% groups
     if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
          stop(
               "Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (",
               length(x = labels.loc), ")."
          )
     }
     names(x = labels) <- groups
     for (group in groups) {
          labels.loc[labels.loc[, id] == group, id] <- labels[group]
     }
     if (box) {
          geom.use <- ifelse(test = repel, yes = geom_label_repel,
                             no = geom_label
          )
          plot <- plot + geom.use(
               data = labels.loc, mapping = aes_string(
                    x = xynames["x"],
                    y = xynames["y"], label = id, fill = id
               ), show.legend = FALSE,
               ...
          ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[
               ,
               id
          ])])
     } else {
          geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
          plot <- plot + geom.use(
               data = labels.loc, mapping = aes_string(
                    x = xynames["x"],
                    y = xynames["y"], label = id
               ), max.iter=max.iter,
               color = color, bg.color = bg.color, bg.r = bg.r, show.legend = FALSE,
               ...
          )
     }
     if (geom == "GeomSpatial") {
          plot <- suppressMessages(expr = plot + coord_fixed(
               xlim = xrange.save,
               ylim = yrange.save
          ))
     }
     return(plot)
}









############################################################################


dplabelClustersWithOutilnedText <- function(plot, id, clusters = NULL, labels = NULL, split.by = NULL, repel = TRUE, box = FALSE, geom = "GeomPoint", position = "median", color = "white", bg.r = 0.1, max.iter=1,
                                            ...) {
     xynames <- unlist(
          x = GetXYAesthetics(plot = plot, geom = geom),
          use.names = TRUE
     )
     if (!id %in% colnames(x = plot$data)) {
          stop("Cannot find variable ", id, " in plotting data")
     }
     if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
          warning("Cannot find splitting variable ", id, " in plotting data")
          split.by <- NULL
     }
     data <- plot$data[, c(xynames, id, split.by)]
     possible.clusters <- as.character(x = na.omit(object = unique(x = data[
          ,
          id
     ])))
     groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[
          ,
          id
     ])))
     if (any(!groups %in% possible.clusters)) {
          stop("The following clusters were not found: ", paste(groups[!groups %in%
                                                                            possible.clusters], collapse = ","))
     }
     pb <- ggplot_build(plot = plot)
     if (geom == "GeomSpatial") {
          xrange.save <- layer_scales(plot = plot)$x$range$range
          yrange.save <- layer_scales(plot = plot)$y$range$range
          data[, xynames["y"]] <- max(data[, xynames["y"]]) - data[
               ,
               xynames["y"]
          ] + min(data[, xynames["y"]])
          if (!pb$plot$plot_env$crop) {
               y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) -
                    pb$layout$panel_params[[1]]$y.range
               data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
          }
     }
     data <- cbind(data, color = pb$data[[1]][[1]])
     labels.loc <- lapply(X = groups, FUN = function(group) {
          data.use <- data[data[, id] == group, , drop = FALSE]
          data.medians <- if (!is.null(x = split.by)) {
               do.call(what = "rbind", args = lapply(X = unique(x = data.use[
                    ,
                    split.by
               ]), FUN = function(split) {
                    medians <- apply(
                         X = data.use[data.use[, split.by] ==
                                           split, xynames, drop = FALSE], MARGIN = 2,
                         FUN = median, na.rm = TRUE
                    )
                    medians <- as.data.frame(x = t(x = medians))
                    medians[, split.by] <- split
                    return(medians)
               }))
          } else {
               as.data.frame(x = t(x = apply(X = data.use[, xynames,
                                                          drop = FALSE
               ], MARGIN = 2, FUN = median, na.rm = TRUE)))
          }
          data.medians[, id] <- group
          data.medians$color <- data.use$color[1]
          return(data.medians)
     })
     if (position == "nearest") {
          labels.loc <- lapply(X = labels.loc, FUN = function(x) {
               group.data <- data[as.character(x = data[, id]) ==
                                       as.character(x[3]), ]
               nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(
                    1,
                    2
               )]), k = 1)$nn.idx
               x[1:2] <- group.data[nearest.point, 1:2]
               return(x)
          })
     }
     labels.loc <- do.call(what = "rbind", args = labels.loc)
     labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[
          ,
          id
     ]))
     labels <- labels %||% groups
     if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
          stop(
               "Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (",
               length(x = labels.loc), ")."
          )
     }
     names(x = labels) <- groups
     for (group in groups) {
          labels.loc[labels.loc[, id] == group, id] <- labels[group]
     }
     if (box) {
          geom.use <- ifelse(test = repel, yes = geom_label_repel,
                             no = geom_label
          )
          plot <- plot + geom.use(
               data = labels.loc, mapping = aes_string(
                    x = xynames["x"],
                    y = xynames["y"], label = id, fill = id
               ), show.legend = FALSE,
               ...
          ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[
               ,
               id,
          ])])
     } else {
          geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
          #labels.loc <- labels.loc[match(0:(dim(labels.loc)[[1]]-1),labels.loc$ident),]
          plot <- plot + geom.use(
               data = labels.loc, mapping = aes_string(
                    x = xynames["x"],
                    y = xynames["y"], label = id
               ), max.iter=max.iter,
               color = color, bg.color = labels.loc$color[order(labels.loc[,id,])], bg.r = bg.r, show.legend = FALSE,
               ...
          )
     }
     if (geom == "GeomSpatial") {
          plot <- suppressMessages(expr = plot + coord_fixed(
               xlim = xrange.save,
               ylim = yrange.save
          ))
     }
     return(plot)
}










############################################################################


dp <- function(object, dims = c(1, 2), cells = NULL, cols = NULL, pt.size = 1, reduction = NULL, group.by = NULL, split.by = NULL, shape.by = NULL, order = NULL, shuffle = FALSE, seed = 1, label = TRUE, label.size = 7, label.color = "black", bg.r=0.1, label.box = FALSE, max.iter=1, repel = TRUE, alpha = 0.5, cells.highlight = NULL, cols.highlight = "#DE2D26", sizes.highlight = 1, autosize = FALSE, na.value = "grey50", ncol = NULL, combine = TRUE, raster = NULL, raster.dpi = c(512, 512)) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[["ident"]] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  data <- cbind(data, object[[group.by]][cells, , drop = FALSE])
  group.by <- colnames(x = data)[3:ncol(x = data)]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  plots <- lapply(X = group.by, FUN = function(x) {
  x="ident"
    plot <- mySingleDimPlot(
      data = data[, c(
        dims, x, split.by,
        shape.by
      )], dims = dims, col.by = x, cols = cols,
      my.pt.size = pt.size, shape.by = shape.by, order = order,
      label = FALSE, alpha=alpha,cells.highlight = cells.highlight,
      cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, autosize = autosize,
      na.value = na.value, raster = raster, raster.dpi = raster.dpi
    )
    if (label) {
      plot <- dplabelClustersWithOutilnedText(plot = plot, id = x, repel = repel, size = label.size, split.by = split.by, box = label.box, color = label.color, alpha=1,max.iter=max.iter)
    }
    if (!is.null(x = split.by)) {
      plot <- plot + Seurat:::FacetTheme() + facet_wrap(
        facets = vars(!!sym(x = split.by)),
        ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
          length(x = unique(x = data[, split.by]))
        } else {
          ncol
        }
      )
    }
    plot <- if (is.null(x = orig.groups)) {
      plot + labs(title = NULL)
    } else {
      plot + CenterTitle()
    }
  })
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}










############################################################################



mySingleDimPlot <- function(data, dims, col.by = NULL, cols = NULL, pt.size = NULL, my.pt.size = my.pt.size,
                            shape.by = NULL, alpha.by = NULL, order = NULL, label = FALSE,
                            repel = FALSE, label.size = 4, cells.highlight = NULL, cols.highlight = "#DE2D26", alpha = 1,
                            sizes.highlight = 1, na.value = "grey50", raster = NULL,
                            raster.dpi = NULL, autosize=FALSE) {
     if (autosize==TRUE){
          pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
     }
     else {
          pt.size <- my.pt.size}
     if ((nrow(x = data) > 1e+05) & !isFALSE(raster)) {
          message(
               "Rasterizing points since number of points exceeds 100,000.",
               "\nTo disable this behavior set `raster=FALSE`"
          )
     }
     raster <- raster %||% (nrow(x = data) > 1e+05)
     if (!is.null(x = raster.dpi)) {
          if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) !=
              2) {
               stop("'raster.dpi' must be a two-length numeric vector")
          }
     }
     if (length(x = dims) != 2) {
          stop("'dims' must be a two-length vector")
     }
     if (!is.data.frame(x = data)) {
          data <- as.data.frame(x = data)
     }
     if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
          stop("Cannot find dimensions to plot in data")
     } else if (is.numeric(x = dims)) {
          dims <- colnames(x = data)[dims]
     }
     if (!is.null(x = cells.highlight)) {
          highlight.info <- Seurat:::SetHighlight(
               cells.highlight = cells.highlight,
               cells.all = rownames(x = data), sizes.highlight = pt.size,
               pt.size = pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||%
                    "#C3C3C3"
          )
          order <- highlight.info$plot.order
          data$highlight <- highlight.info$highlight
          col.by <- "highlight"
          pt.size <- highlight.info$size
          cols <- highlight.info$color
     }
     if (!is.null(x = order) && !is.null(x = col.by)) {
          if (typeof(x = order) == "logical") {
               if (order) {
                    data <- data[order(
                         !is.na(x = data[, col.by]),
                         data[, col.by]
                    ), ]
               }
          } else {
               order <- rev(x = c(order, setdiff(x = unique(x = data[
                    ,
                    col.by
               ]), y = order)))
               data[, col.by] <- factor(x = data[, col.by], levels = order)
               new.order <- order(x = data[, col.by])
               data <- data[new.order, ]
               if (length(x = pt.size) == length(x = new.order)) {
                    pt.size <- pt.size[new.order]
               }
          }
     }
     if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
          warning("Cannot find ", col.by, " in plotting data, not coloring plot")
          col.by <- NULL
     } else {
          col.index <- match(x = col.by, table = colnames(x = data))
          if (grepl(pattern = "^\\d", x = col.by)) {
               col.by <- paste0("x", col.by)
          } else if (grepl(pattern = "-", x = col.by)) {
               col.by <- gsub(
                    pattern = "-", replacement = ".",
                    x = col.by
               )
          }
          colnames(x = data)[col.index] <- col.by
     }
     if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
          warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
     }
     if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
          warning("Cannot find alpha variable ", alpha.by, " in data, setting to NULL",
                  call. = FALSE, immediate. = TRUE
          )
          alpha.by <- NULL
     }
     plot <- ggplot(data = data)
     plot <- if (isTRUE(x = raster)) {
          plot + geom_scattermore(mapping = aes_string(
               x = dims[1],
               y = dims[2], color = paste0("`", col.by, "`"), shape = shape.by,
               alpha = alpha.by
          ), pointsize = pt.size, pixels = raster.dpi)
     } else {
          plot + geom_point(mapping = aes_string(
               x = dims[1], y = dims[2],
               color = paste0("`", col.by, "`"), shape = shape.by
               
          ), size = pt.size, alpha = alpha)
     }
     plot <- plot + guides(color = guide_legend(override.aes = list(size = 3))) +
          labs(color = NULL, title = col.by, ) + CenterTitle()
     if (label && !is.null(x = col.by)) {
          plot <- LabelClusters(
               plot = plot, id = col.by, repel = repel,
               size = label.size
          )
     }
     if (!is.null(x = cols)) {
          if (length(x = cols) == 1 && (is.numeric(x = cols) ||
                                        cols %in% rownames(x = brewer.pal.info))) {
               scale <- scale_color_brewer(palette = cols, na.value = na.value)
          } else if (length(x = cols) == 1 && (cols %in% c(
               "alphabet",
               "alphabet2", "glasbey", "polychrome", "stepped"
          ))) {
               colors <- DiscretePalette(length(unique(data[[col.by]])),
                                         palette = cols
               )
               scale <- scale_color_manual(values = colors, na.value = na.value)
          } else {
               scale <- scale_color_manual(values = cols, na.value = na.value)
          }
          plot <- plot + scale
     }
     plot <- plot + theme_cowplot()
     return(plot)
}
