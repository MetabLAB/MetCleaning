#' @title RSDoverview
#' @description Evaluate the RSD information before and after processing.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData.before MetFlowData.before.
#' @param MetFlowData.after MetFlowData.after.
#' @param path Work directory.
#' @return RSD comparation.

RSDoverview <- function(MetFlowData.before = MetFlowData1,
                        MetFlowData.after = MetFlowData2,
                        path = NULL) {
  if (is.null(path)) {
    path <- getwd()
  } else{
    dir.create(path)
  }

  hasQC <- MetFlowData.before[["hasQC"]]
  if (hasQC == 'no') {
    stop("Data has no QC!!!")
  }
  ##before
  subject.bef <- MetFlowData.before[["subject"]]
  qc.bef <- MetFlowData.before[["qc"]]
  subject.info.bef <- MetFlowData.before[["subject.info"]]
  qc.info.bef <- MetFlowData.before[["qc.info"]]

  ##split data
  data.bef <- SplitBatch(MetFlowData = MetFlowData.before)
  subject.bef1 <- data.bef[[1]]
  qc.bef1 <- data.bef[[2]]

  ##after
  subject.aft <- MetFlowData.after[["subject"]]
  qc.aft <- MetFlowData.after[["qc"]]
  subject.info.aft <- MetFlowData.after[["subject.info"]]
  qc.info.aft <- MetFlowData.after[["qc.info"]]

  ##split data
  data.aft <- SplitBatch(MetFlowData = MetFlowData.after)
  subject.aft1 <- data.aft[[1]]
  qc.aft1 <- data.aft[[2]]


  pdf(
    file.path(path, "RSD distribution in different batch.pdf"),
    width = 14,
    height = 7
  )
  layout(matrix(c(1:2), ncol = 2))
  par(mar = c(5, 5, 4, 2))
  for (i in 1:length(qc.bef1)) {
    #before
    QC <- qc.bef1[[i]]
    QC.rsd <- apply(QC, 1, function(x) {
      sd(x) * 100 / mean(x)
    })
    par(mar = c(5, 5, 4, 2))
    colours1 <- rep(NA, length(QC.rsd))
    colours1[QC.rsd > 30] <- "tomato"
    colours1[is.na(colours1)] <- "grey"

    plot(
      QC.rsd,
      xlab = "Feature index",
      ylab = "Relative Standard Deviation (RSD, %)",
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      col = colours1,
      main = paste("Batch", i, "QC (Before)"),
      cex.main = 1.3
    )
    abline(h = 30, lty = 2)
    legend("topleft",
           paste("RSD<30%: ", round(sum(QC.rsd < 30) / length(QC.rsd), 4) * 100, "%"),
           bty = "n",
           cex = 1.3)

    #after
    QC <- qc.aft1[[i]]
    QC.rsd <- apply(QC, 1, function(x) {
      sd(x) * 100 / mean(x)
    })
    par(mar = c(5, 5, 4, 2))
    colours1 <- rep(NA, length(QC.rsd))
    colours1[QC.rsd > 30] <- "tomato"
    colours1[is.na(colours1)] <- "grey"

    plot(
      QC.rsd,
      xlab = "Feature index",
      ylab = "Relative Standard Deviation (RSD, %)",
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      col = colours1,
      main = paste("Batch", i, "QC (After)"),
      cex.main = 1.3
    )
    abline(h = 30, lty = 2)
    legend("topleft",
           paste("RSD<30%: ", round(sum(QC.rsd < 30) / length(QC.rsd), 4) * 100, "%"),
           bty = "n",
           cex = 1.3)

  }
  dev.off()

  pdf(
    file.path(path, "RSD distribution in all batches.pdf"),
    width = 14,
    height = 7
  )
  layout(matrix(c(1:2), ncol = 2))
  par(mar = c(5, 5, 4, 2))
  # before
  qc.rsd <- apply(qc.bef, 1, function(x) {
    sd(x) * 100 / mean(x)
  })
  colours1 <- rep(NA, length(qc.rsd))
  colours1[qc.rsd > 30] <- "tomato"
  colours1[is.na(colours1)] <- "grey"
  plot(
    qc.rsd,
    xlab = "Feature index",
    ylab = "Relative Standard Deviation (RSD, %)",
    cex.lab = 1.3,
    cex.axis = 1.3,
    pch = 19,
    col = colours1,
    cex.main = 1.3,
    main = "Before"
  )
  abline(h = 30, lty = 2)
  legend("topleft",
         paste("RSD<30%: ", round(sum(qc.rsd < 30) / length(qc.rsd), 4) * 100, "%"),
         bty = "n",
         cex = 1.3)

  # after
  qc.rsd <- apply(qc.aft, 1, function(x) {
    sd(x) * 100 / mean(x)
  })
  colours1 <- rep(NA, length(qc.rsd))
  colours1[qc.rsd > 30] <- "tomato"
  colours1[is.na(colours1)] <- "grey"
  plot(
    qc.rsd,
    xlab = "Feature index",
    ylab = "Relative Standard Deviation (RSD, %)",
    cex.lab = 1.3,
    cex.axis = 1.3,
    pch = 19,
    col = colours1,
    cex.main = 1.3,
    main = "After"
  )
  abline(h = 30, lty = 2)
  legend("topleft",
         paste("RSD<30%: ", round(sum(qc.rsd < 30) / length(qc.rsd), 4) * 100, "%"),
         bty = "n",
         cex = 1.3)

  dev.off()
  layout(1)
}