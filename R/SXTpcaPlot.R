#' @title SXTpcaPlot
#' @description Draw PCA score or loading plot for SXTpcaData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param SXTpcaData SXTpcaPlot subject for drawing plots.
#' @param plot Which class of plot you want to draw, "score" or "loading", default is "score".
#' @param loading.plot.name The name of loading plot.
#' @param score.plot.name The name of score plot.
#' @param color The colors you want to use in different classs.
#' @param shape The shapes you want to use in different classes.
#' @param text Add text or not? Default is FALSE.
#' @param ellipse Draw ellipse or not?
#' @param path Directory you want to output results.
#' @return PCA score and loading plot.
#' @export

##plot functioin for SXTpcaData
SXTpcaPlot <- function(SXTpcaData = SXTpcaData,
                       plot = "score",
                       loading.plot.name = "PCA loading",
                       score.plot.name = "PCA score",
                       color = c("palegreen","royalblue","firebrick1","tan1","deepskyblue",
                                 "cyan","gray48", "chocolate4","darkmagenta","indianred1"),
                       shape = c(17,19,15,18,2,8,11,13,12,14),
                       cexlab = 1.3,
                       cexaxis = 1.3,
                       cexa = 1.3,
                       cextext = 1,
                       width = 7,
                       height = 7,
                       text = FALSE,
                       ellipse = FALSE,
                       path = NULL) {

  # browser()
  if (is.null(path)) path <- getwd()
  ## get data
  subject <- SXTpcaData[["subject"]]
  qc <- SXTpcaData[["qc"]]
  sample.pca = SXTpcaData[["sample.pca"]]
  info = SXTpcaData[["info"]]
  QC = SXTpcaData[["QC"]]
  int.pca <- sample.pca


  loading <- summary(int.pca)$rotation
  pov <- summary(int.pca)$importance[2, ]
  sd <- summary(int.pca)$importance[1, ]
  cp <- summary(int.pca)$importance[3, ]
  pc <- int.pca$x

  pc1 <- pov[1]
  pc2 <- pov[2]
  pc3 <- pov[3]

  x <- pc[, 1]
  y <- pc[, 2]
  z <- pc[, 3]

  xmin <- 1.2 * min(x)
  xmax <- 1.2 * max(x)
  ymin <- 1.2 * min(y)
  ymax <- 1.2 * max(y)
  zmin <- 1.2 * min(z)
  zmax <- 1.2 * max(z)

  name <- colnames(subject)
  label <- list()
  for (i in 1:length(info)) {
    label[[i]] <- match(as.character(info[[i]]), name)
    label[[i]] <- label[[i]][!is.na(label[[i]])]
  }

  legend <- NULL
  for (i in 1:length(label)) {
    legend[label[[i]]] <- names(info)[i]
  }

  if (QC) {
    legend[(ncol(subject)+1):(ncol(subject)+ncol(qc))] <- "QC"
  }

  colour <- NULL
  if (length(color) < length(info))
    stop("Color list is not enough")

  colourlist <- color
  for (i in 1:length(label)) {
    colour[label[[i]]] <- colourlist[i]
  }
  if (QC) {
    colour[(ncol(subject)+1):(ncol(subject)+ncol(qc))] <- colourlist[length(info) + 1]
  }

  pcha <- NULL
  if (length(shape) < length(info))
    stop("shape list is not enough")
  pchalist <- shape
  for (i in 1:length(label)) {
    pcha[label[[i]]] <- pchalist[i]
  }
  if (QC) {
    pcha[(ncol(subject)+1):(ncol(subject)+ncol(qc))] <- pchalist[length(info) + 1]
  }

  if (plot == "loading") {
    #laoding plot
    pdf(file.path(path, paste(loading.plot.name,".pdf", sep = "")),
        width = width,
        height = height)
    plot(
      loading[, 1],
      loading[, 2],
      pch = 20,
      xlab = "Component 1",
      ylab = "Component 2",
      cex.lab = cexlab,
      cex.axis = cexaxis
    )
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    plot(
      loading[, 2],
      loading[, 3],
      pch = 20,
      xlab = "Component 2",
      ylab = "Component 3",
      cex.lab = cexlab,
      cex.axis = cexaxis
    )
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    plot(
      loading[, 1],
      loading[, 3],
      pch = 20,
      xlab = "Component 1",
      ylab = "Component 3",
      cex.lab = cexlab,
      cex.axis = cexaxis
    )
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    #loading plot 3d
    scatterplot3d(
      loading[, 1],
      loading[, 2],
      loading[, 3],
      xlab = "Component 1",
      ylab = "Component 2",
      zlab = "Component 3",
      angle = 40,
      pch = 20,
      box = FALSE,
      cex.symbol = 1,
      cex.lab = 1.3,
      cex.axis = 0.8
    )
    dev.off()
  }

  if (plot == "score") {
    #PCA 2D
    pdf(file.path(path, paste(score.plot.name,".pdf", sep = "")),
        width = width,
        height = height)
    #t1plot
    par(mar = c(5, 5, 4, 2))
    plot(
      y,
      ylim = c(ymin, ymax),
      col = colour,
      pch = pcha,
      xlab = "QC injection order",
      ylab = paste("PC1:", pc1, sep = ""),
      cex = cexa,
      cex.axis = cexaxis,
      cex.lab = cexlab
    )

    abline(h = 0, lty = 2)


    #t1 vs t2 plot
    par(mar = c(5, 5, 4, 2))
    plot(
      x,
      y,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      col = colour,
      pch = pcha,
      xlab = paste("PC1:", pc1, sep = ""),
      ylab = paste("PC2:", pc2, sep = ""),
      cex = cexa,
      cex.axis = cexaxis,
      cex.lab = cexlab
    )
    if (text)
    {
      text(x, y, name, pos = 4, cex = cextext)
    }

    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    if (ellipse)
    {
      lines(ellipse::ellipse(
        0,
        scale = c(sd(x), sd(y)),
        centre = c(mean(x), mean(y))
      ), lty = 2)
    }

    if (QC) {
      legend(
        "topleft",
        c(names(info), "QC"),
        pch = pchalist[1:(length(info) + 1)],
        col = colourlist[1:(length(info) + 1)],
        bty = "n",
        cex = 1.1
      )
    }
    else{
      legend(
        "topleft",
        names(info),
        pch = pchalist[1:length(info)],
        col = colourlist[1:length(info)],
        bty = "n",
        cex = 1.1
      )
    }

    #t2 vs t3 plot

    par(mar = c(5, 5, 4, 2))
    plot(
      y,
      z,
      xlim = c(ymin, ymax),
      ylim = c(zmin, zmax),
      col = colour,
      pch = pcha,
      xlab = paste("PC2:", pc2, sep = ""),
      ylab = paste("PC3:", pc3, sep = ""),
      cex = cexa,
      cex.axis = cexaxis,
      cex.lab = cexlab
    )
    if (text)
    {
      text(y, z, name, pos = 4, cex = cextext)
    }
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
    if (ellipse)
    {
      lines(ellipse::ellipse(
        0,
        scale = c(sd(y), sd(z)),
        centre = c(mean(y), mean(z))
      ), lty = 2)
    }


    if (QC) {
      legend(
        "topleft",
        c(names(info), "QC"),
        pch = pchalist[1:(length(info) + 1)],
        col = colourlist[1:(length(info) + 1)],
        bty = "n",
        cex = 1.1
      )
    }
    else{
      legend(
        "topleft",
        names(info),
        pch = pchalist[1:length(info)],
        col = colourlist[1:length(info)],
        bty = "n",
        cex = 1.1
      )
    }

    #t1 vs t3 plot
    plot(x, z,
         xlim = c(xmin, xmax),
         ylim = c(zmin, zmax),
         col = colour,
         pch = pcha,
         xlab = paste("PC1:", pc1, sep = ""),
         ylab = paste("PC3:", pc3, sep = ""),
         cex = 1.3,
         cex.axis = 1.3,
         cex.lab = 1.3
    )
    if (text)
    {
      text(x, z, name, pos = 4, cex = cextext)
    }
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
    if (ellipse)
    {
      lines(ellipse::ellipse(
        0,
        scale = c(sd(x), sd(z)),
        centre = c(mean(x), mean(z))
      ), lty = 2)
    }


    if (QC) {
      legend(
        "topleft",
        c(names(info), "QC"),
        pch = pchalist[1:(length(info) + 1)],
        col = colourlist[1:(length(info) + 1)],
        bty = "n",
        cex = 1.1
      )
    }
    else{
      legend(
        "topleft",
        names(info),
        pch = pchalist[1:length(info)],
        col = colourlist[1:length(info)],
        bty = "n",
        cex = 1.1
      )
    }

    #PCA 3D
    scatterplot3d(x, y, z,
                  color = colour,
                  xlab = paste("PC1:", pc1, sep = ""),
                  ylab = paste("PC2:", pc2, sep = ""),
                  zlab = paste("PC3:", pc3, sep = ""),
                  angle = 40,
                  pch = pcha,
                  box = FALSE,
                  cex.symbol = 1,
                  cex.lab = 1.3,
                  cex.axis = 1.3,
                  xlim = c(xmin, xmax),
                  ylim = c(ymin, ymax),
                  zlim = c(zmin, zmax)
    )

    if (QC) {
      legend(
        "topleft",
        c(names(info), "QC"),
        pch = pchalist[1:(length(info) + 1)],
        col = colourlist[1:(length(info) + 1)],
        bty = "n",
        cex = 1.5
      )
    }
    else{
      legend(
        "topleft",
        names(info),
        pch = pchalist[1:length(info)],
        col = colourlist[1:length(info)],
        bty = "n",
        cex = 1.5
      )
    }
    dev.off()
  }
}