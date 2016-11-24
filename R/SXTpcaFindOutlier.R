#' @title SXTpcaFindOutlier
#' @description Find which samples are in outside of 95\% confidence interval in PCA score plot.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param SXTpcaData SXTpcaData subject for outliers finding.
#' @param CI Confidence interval.
#' @param output.plot Output PCA score plot or not.Default is TRUE.
#' @param plot.name The name of score plot.
#' @param path Directory you want to output results.
#' @return It return a outlier index.
#' @export

##find outlier functions using SXTpcaData
SXTpcaFindOutlier <- function(SXTpcaData = SXTpcaData,
                              CI = 0.95,
                              output.plot = TRUE,
                              plot.name = "PCA score plot for outliers",
                              path = NULL) {
  # browser()
  if (is.null(path))
    path <- getwd()
  sample.pca <- SXTpcaData[["sample.pca"]]
  data <- SXTpcaData[["subject"]]

  loading <- summary(sample.pca)$rotation
  pov <- summary(sample.pca)$importance[2,]
  sd <- summary(sample.pca)$importance[1,]
  cp <- summary(sample.pca)$importance[3,]
  pc <- sample.pca$x

  pc1 <- round(pov[1], 2)
  pc2 <- round(pov[2], 2)
  pc3 <- round(pov[3], 2)

  x <- pc[, 1]
  y <- pc[, 2]
  z <- pc[, 3]

  xmin <- 1.2 * min(x)
  xmax <- 1.2 * max(x)
  ymin <- 1.2 * min(y)
  ymax <- 1.2 * max(y)
  zmin <- 1.2 * min(z)
  zmax <- 1.2 * max(z)

  ellipse.data <-
    SXTellipse(cor(x, y),
               scale = c(sd(x), sd(y)),
               centre = c(mean(x), mean(y)))
  data.for.plot <- ellipse.data[["ellipse.data"]]
  t = ellipse.data[["t"]]
  a = ellipse.data[["a"]]
  d = ellipse.data[["d"]]
  scale = ellipse.data[["scale"]]
  centre = ellipse.data[["centre"]]

  ## right point
  point1 <- c(t * scale[1] * cos(0) + centre[1], t * scale[2] *
                cos(-d) + centre[2])
  ## top point
  point2 <- c(t * scale[1] * cos(pi / 2) + centre[1], t * scale[2] *
                cos(pi / 2 - d) + centre[2])
  ## left point
  point3 <- c(t * scale[1] * cos(pi) + centre[1], t * scale[2] *
                cos(pi - d) + centre[2])

  ellipse.a <- sqrt(sum((point1 - centre) ^ 2))
  ellipse.b <- sqrt(sum((point2 - centre) ^ 2))
  ellipse.c <- sqrt(ellipse.a ^ 2 - ellipse.b ^ 2)

  ## get the focus points

  lm.reg <- lm(c(point1[2], centre[2]) ~ c(point1[1], centre[1]))
  a <- lm.reg$coefficients[[2]]
  b <- lm.reg$coefficients[[1]]

  foo.f.a <- a ^ 2 + 1
  foo.f.b <-  2 * a * (b - centre[2]) - 2 * centre[1]
  foo.f.c <- centre[1] ^ 2 + (b - centre[2]) ^ 2 - ellipse.c ^ 2

  foo.f <- function(x,
                    a = foo.f.a,
                    b = foo.f.b,
                    c = foo.f.c)
  {
    a * x ^ 2 + b * x + c
  }
  result1 <-
    uniroot(
      foo.f,
      c(0, 10000),
      a = foo.f.a,
      b = foo.f.b,
      c = foo.f.c,
      tol = 0.0001
    )
  result2 <-
    uniroot(
      foo.f,
      c(-10000, 0),
      a = foo.f.a,
      b = foo.f.b,
      c = foo.f.c,
      tol = 0.0001
    )

  p1 <- c(result1$root, foo.f(x = result1$root))
  p2 <- c(result2$root, foo.f(x = result2$root))

  x1 <- data.for.plot[, 1]
  y1 <- data.for.plot[, 2]
  ellipse.standard <-
    mean(sqrt((x1 - p1[1]) ^ 2 + (y1 - p1[2]) ^ 2) + sqrt((x1 - p2[1]) ^ 2 + (y1 -
                                                                                p2[2]) ^ 2))

  distance <-
    sqrt((x - p1[1]) ^ 2 + (y - p1[2]) ^ 2) + sqrt((x - p2[1]) ^ 2 + (y - p2[2]) ^
                                                     2)

  outlier.index <- which(distance > ellipse.standard)
  if (length(outlier.index) > 0) {
    cat(names(x)[outlier.index], " are outliers!!!\n")
  }

  if (output.plot) {
    colour <- rep(NA, length(x))
    colour[outlier.index] <- "firebrick1"
    colour[is.na(colour)] <- "black"
    pdf(file.path(path, paste(plot.name, ".pdf", sep = "")),
        width = 7,
        height = 7)
    par(mar = c(5,5,4,2))
    plot(
      x,
      y,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      xlab = paste("PC1:", pc1),
      ylab = paste("PC2:", pc2),
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      col = colour,
      main = "PCA score plot for outliers"
    )
    if (length(outlier.index) > 0) {
      text(x = x[outlier.index],
           y = y[outlier.index],
           names(x)[outlier.index],
           pos = 4)
    } else
      (text(x = x,
            y = y,
            names(x),
            pos = 4))
    abline(h = 0, lty = 2)
    abline(v = 0, lty = 2)
    lines(data.for.plot, lty = 2)
    dev.off()
  }

  return(outlier.index)

}