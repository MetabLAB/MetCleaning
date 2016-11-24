#' @title PLSanalysis
#' @description PLS analysis for MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param color Color for different class.
#' @param shape Shape for different class.
#' @param scale.method Whihch scale methd you want to use? "auto" or "pareto",
#' defaulit is "auto".
#' @param path Work directory.
#' @param text Add text in PCA score plot or not? Deefault is FALSE.
#' @param ellipse Add ellipse in PCA score plot or not? Deefault is TRUE.
#' @param xlim1_ylim1 The x and y axis limitation. Default is NULL.
#' @return PLS score plot: PLS score plot.
#' @return permutation test plot: Permutation test plot.
#' @export

PLSanalysis <- function(MetFlowData = MetFlowData,
                        #used data
                        log.scale = 10,
                        scalemethod = "auto",
                        plsmethod = "plsreg",
                        path = NULL,
                        width = 7,
                        height = 7,
                        text = FALSE,
                        ellipse = TRUE,
                        color = c("palegreen",
                                  "firebrick1",
                                  "royalblue",
                                  "yellow",
                                  "black",
                                  "cyan",
                                  "gray48"),
                        shape = c(17, 19, 15, 18, 2, 8, 11),
                        cexa = 1,
                        xlim1 = NULL,
                        ylim1 = NULL)
#parameter setting
{
  # browser()
  # requireNamespace("pls")
  if (is.null(path)) {
    path <- getwd()
  } else{
    dir.create(path)
  }
  options(warn = -1)

  subject <- MetFlowData[["subject"]]
  qc <- MetFlowData[["qc"]]
  subject.info <- MetFlowData[["subject.info"]]
  group <- subject.info[, "group"]
  group.unique <- sort(unique(group))
  subject.name <- subject.info[, 1]

  if (is.null(qc)) {
    QC <- FALSE
  }

  info <- list()
  for (i in 1:length(group.unique)) {
    info[[i]] <- subject.name[which(group == group.unique[i])]
  }

  names(info) <- group.unique

  int <- t(subject)
  index <- NULL
  for (i in 1:length(info)) {
    index1 <- as.character(info[[i]])
    index <- c(index, index1)
  }
  if (length(which(index == "")) != 0)  {
    index <- index[-which(index == "")]
  }


  index <- index[!is.na(index)]
  index <- match(index, rownames(int))
  index <- index[!is.na(index)]
  int <- int[index, ]

  #######
  name <- rownames(int)
  # browser()
  Y <- NULL
  label <- list()
  for (i in 1:length(info)) {
    label[[i]] <- match(as.character(info[[i]]), name)
    label[[i]] <- label[[i]][!is.na(label[[i]])]
    Y[label[[i]]] <- i - 1
  }

  ##log transformation
  if (log.scale == FALSE) {
    int <- int
  }

  if (log.scale == "e") {
    int <- log(int + 1)
  }

  if (log.scale != FALSE & log.scale != "e") {
    int <- log(int + 1, as.numeric(log.scale))
  }

  int.scale <- SXTscale(int, method = scalemethod)
  # int.Y<-SXTscale(Y,method=scalemethod)
  int.Y <- Y
  # browser()
  ncompa <- min(nrow(int), ncol(int)) - 1

  if (plsmethod == "plsr") {
    pls1 <-
      pls::plsr(
        int.Y ~ int.scale,
        scale = FALSE,
        validation = "CV",
        ncomp = ncompa,
        method = "oscorespls"
      )
    save(pls1, file = "pls1")

    #########select the number of compents#################
    msep <- pls::MSEP(pls1)
    save(msep, file = file.path(path, "msep"))
    msep <- msep$val[, ,]

    yesnot <- "y"
    while (yesnot == "y" | yesnot == "") {
      comps.number <- readline("How many comps do you want to see?")
      while (!exists("comps.number") | comps.number == "")
      {
        cat("You must give a comps number to continute!!!\n")
        comps.number <-
          readline("How many comps do you want to see?")
      }
      comps.number <- as.numeric(comps.number)
      par(mar = c(5,5,4,2))
      plot(
        x = c(1:comps.number),
        y = msep[1, 2:(comps.number + 1)],
        type = "b",
        col = "tomato",
        pch = 19,
        xlab = "ncomp",
        ylab = "MSEP",
        cex.lab = 1.5,
        cex.axis = 1.3,
        ylim= c(0.8 * min(msep[c(1, 2), 2:(comps.number + 1)]),
                1.2 * max(msep[c(1, 2), 2:(comps.number + 1)]))
      )
      points(
        x = c(1:comps.number),
        y = msep[2, 2:(comps.number + 1)],
        type = "b",
        pch = 17,
        col = "grey"
      )
      legend(
        "top",
        legend = c("CV", "adjCV"),
        col = c("tomato", "grey"),
        pch = c(19, 19),
        bty = "n",
        cex = 1.3,
        pt.cex = 1.3
      )
      yesnot <- readline("Do you want to see the next plot? (y/n)")
    }

    pdf(file.path(path, "MSEP plot.pdf"))
    par(mar = c(5,5,4,2))
    plot(
      x = c(1:comps.number),
      y = msep[1, 2:(comps.number + 1)],
      type = "b",
      col = "tomato",
      pch = 19,
      xlab = "ncomp",
      ylab = "MSEP",
      cex.lab = 1.5,
      cex.axis = 1.3,
      ylim= c(0.99 * min(msep[c(1, 2), 2:(comps.number + 1)]),
              1.1 * max(msep[c(1, 2), 2:(comps.number + 1)]))
    )
    points(
      x = c(1:comps.number),
      y = msep[2, 2:(comps.number + 1)],
      type = "b",
      pch = 17,
      col = "grey"
    )
    legend(
      "top",
      legend = c("CV", "adjCV"),
      col = c("tomato", "grey"),
      pch = c(19, 19),
      bty = "n",
      cex = 1.3,
      pt.cex = 1.3
    )
    dev.off()

    number <-
      readline("Please type number and press Enter to continute:")
    while (!exists("number") | number == "") {
      cat("You must give a number to continute!!!\n")
      number <-
        readline("Please type comps number value and press Enter  to continute:")
    }
    number <- as.numeric(number)

    ##################construct final pls model###################
    pls2 <-
      pls::plsr(
        int.Y ~ int.scale,
        scale = FALSE,
        validation = "CV",
        ncomp = number,
        method = "oscorespls"
      )
    save(pls2, file = file.path(path, "pls2"))
    vip <- SXTvip(pls2)
    save(vip, file = file.path(path, "vip"))

    scores <- scores(pls2)
    x <- scores[, 1]
    y <- scores[, 2]
    if (number > 2) {
      z <- scores[, 3]
      zmin <- 1.2 * min(z)
      zmax <- 1.2 * max(z)
    }

    xmin <- 1.2 * min(x)
    xmax <- 1.2 * max(x)
    ymin <- 1.2 * min(y)
    ymax <- 1.2 * max(y)
  }

  else {
    # browser()
    dummy <- SXTdummy(Y)
    # int.dummy<-SXTscale(dummy,method=scalemethod)
    int.dummy <- dummy
    # ncompa = nrow(int.scale) - 1
    ncompa <- min(nrow(int), ncol(int)) - 1
    pls1 <- plsdepot::plsreg1(int.scale, Y, comps = ncompa)
    save(pls1, file = file.path(path, "pls1"))
    #########select the number of compents#################
    Q2cum <- pls1$Q2[, 5]
    Q2cum[is.nan(Q2cum)] <- 1
    yesnot <- "y"
    while (yesnot == "y" | yesnot == "") {
      comps.number <- readline("How many comps do you want to see?")
      while (!exists("comps.number") |
             comps.number == "") {
        cat("You must give a comps number to continute!!!\n")
        comps.number <-
          readline("How many comps do you want to see?")
      }
      comps.number <- as.numeric(comps.number)
      par(mar = c(5,5,4,2))
      a <-
        barplot(
          Q2cum[1:comps.number],
          xlab = "ncomp",
          ylab = "Q2cum",
          cex.lab = 1.5,
          cex.axis = 1.3
        )
      abline(h = 0)
      points(
        a,
        Q2cum[1:comps.number],
        type = "b",
        col = "tomato",
        pch = 20,
        cex = 2
      )
      yesnot <- readline("Do you want to see the next plot? (y/n)")
    }
    pdf(file.path(path, "Q2cum plot.pdf"),
        width = 7,
        height = 7)
    par(mar = c(5,5,4,2))
    a <- barplot(
      Q2cum[1:comps.number],
      xlab = "ncomp",
      ylab = "Q2cum",
      cex.lab = 1.5,
      cex.axis = 1.3
    )
    abline(h = 0)
    points(
      a,
      Q2cum[1:comps.number],
      type = "b",
      col = "tomato",
      pch = 20,
      cex = 2
    )
    dev.off()

    number <-
      readline("Please type number and press Enter to continute:")
    while (!exists("number") |
           number == "") {
      cat("You must give a number to continute!!!\n")
      number <-
        readline("Please type comps number value and press Enter to continute:")
    }
    number <- as.numeric(number)

    ##################construct final pls model###################
    cat(paste(
      "Construct PLS model with all peaks using",
      number,
      "comps ...",
      "\n"
    ))
    pls2 <- plsdepot::plsreg1(int.scale, Y, comps = number)
    save(pls2, file = file.path(path, "pls2"))
    pls.temp <- plsdepot::plsreg2(int.scale, int.dummy, comps = number)
    vip <- pls.temp$VIP
    Q2cum <- pls2$Q2[, 5]
    R2cum <- cumsum(pls2$R2)

    write.csv(cbind(R2cum, Q2cum),
              file.path(path, "R2Q2.csv"),
              row.names = F)

    ##draw barplot of Q2cum, R2Xcum and R2Ycum
    Q2R2 <- cbind(R2cum, Q2cum)
    colnames(Q2R2) <- c("R2cum", "Q2cum")
    pdf(file.path(path, "Q2R2cum.pdf"),
        width = 8,
        height = 6)
    par(mar = c(5,5,4,2))
    barplot(
      t(Q2R2),
      beside = T,
      col = c("royalblue", "tomato"),
      cex.lab = 1.5,
      cex.axis = 1.3,
      cex.names = 1.5
    )
    abline(h = 0)
    legend(
      "topleft",
      legend = c("R2Ycum", "Q2cum"),
      pch = 15,
      col = c("royalblue", "tomato"),
      cex = 1.3,
      pt.cex = 1.3,
      bty = "n"
    )
    dev.off()

    save(vip, file = file.path(path, "vip"))
    save(Q2cum, file = file.path(path, "Q2cum"))
    save(R2cum, file = file.path(path, "R2cum"))

    x <- pls2$x.scores[, 1]
    y <- pls2$x.scores[, 2]
    if (number > 2) {
      z <- pls2$x.scores[, 3]
      zmin <- 1.2 * min(z)
      zmax <- 1.2 * max(z)
    }

    xmin <- 1.2 * min(x)
    xmax <- 1.2 * max(x)
    ymin <- 1.2 * min(y)
    ymax <- 1.2 * max(y)
  }

  legend <- NULL
  for (i in 1:length(label)) {
    legend[label[[i]]] <- names(info)[i]
  }

  colour <- NULL
  colourlist <- color
  for (i in 1:length(label)) {
    colour[label[[i]]] <- colourlist[i]
  }


  pcha <- NULL
  pchalist <- shape
  for (i in 1:length(label)) {
    pcha[label[[i]]] <- pchalist[i]
  }

  if (is.null(xlim1)) {
    xlim = c(xmin, xmax)
  } else {
    xlim = xlim1
  }
  if (is.null(ylim1)) {
    ylim = c(ymin, ymax)
  } else {
    ylim = ylim1
  }

  #PLS 2D
  pdf(file.path(path, "pls score plot.pdf"),
      width = width,
      height = height)
  par(mar = c(5,5,4,2))
  plot(
    x,
    y,
    xlim = xlim,
    ylim = ylim,
    col = colour,
    pch = pcha,
    xlab = "t[1]",
    ylab = "t[2]",
    cex = cexa,
    cex.axis = 1.3,
    cex.lab = 1.35
  )
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  if (text) {
    text(x, y, rownames(int), pos = 4)
  }
  if (ellipse) {
    lines(ellipse::ellipse(
      0,
      scale = c(sd(x), sd(y)),
      centre = c(mean(x), mean(y))
    ), lty = 2)
  }

  legend(
    "topleft",
    names(info),
    pch = pchalist[1:length(info)],
    col = colourlist[1:length(info)],
    bty = "n",
    cex = 1.5
  )
  dev.off()

  #PLS 3D
  if (number > 2) {
    pdf(file.path(path, "pls score plot 3d.pdf"),
        width = width,
        height = height)
    scatterplot3d::scatterplot3d(
      x,
      y,
      z,
      color = colour,
      xlab = "t[1]",
      ylab = "t[2]",
      zlab = "t[3]",
      angle = 50,
      pch = pcha,
      box = FALSE,
      cex.symbol = cexa,
      cex.lab = 1.5,
      cex.axis = 1.3,
      xlim = xlim,
      ylim = ylim,
      zlim = c(zmin, zmax)
    )
    legend(
      "topleft",
      names(info),
      pch = pchalist[1:length(info)],
      col = colourlist[1:length(info)],
      bty = "n",
      cex = 1.5
    )
    dev.off()
  }
  # browser()
  PLSpermutation(
    data = t(subject),
    log.scale = log.scale,
    path = path,
    info = info,
    repeats = 200,
    ncomp = number,
    scalemethod = scalemethod
  )
  cat("\n")
  cat("PLS analysis is done\n")
}
