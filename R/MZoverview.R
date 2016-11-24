#' @title MZOverview
#' @description Evaluate the missing values or zero value information
#' in MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData: MetFlowData.
#' @param path Work Directory.
#' @param what Missing value ("mv") or zero values ("zero").
#' @return Batch i feature MV/zero ratio.csv: Batch i feature MV/zero ratio
#' distribution.
#' @return Batch i QC feature MV/zero ratio.pdf: Batch i QC feature MV.zero
#' ratio distribution.
#' @return Batch i Subject feature MV/zero ratio.pdf: Batch i Subject feature
#'  MV/zero ratio distribution.
#' @return Sample MV/zero distribution.csv: Sample MV/zero ratio distribution.
#' @export

MZoverview <- function(MetFlowData = MetFlowData,
                       path = NULL,
                       what = "mv") {
  # browser()

  if (is.null(path)) {
    path <- getwd()
  }
  else {
    dir.create(path)
  }
  # browser()
  subject <- MetFlowData[["subject"]]
  qc <- MetFlowData[["qc"]]
  tags <- MetFlowData[["tags"]]
  subject.order <- as.numeric(MetFlowData[["subject.order"]])
  qc.order <- as.numeric(MetFlowData[["qc.order"]])

  if (what == "zero" & (sum(is.na(qc)) + sum(is.na(subject))) != 0){
    stop("Please impute MV first.")
  }

  ##variable for MV record
  data <- SplitBatch(MetFlowData = MetFlowData)
  subject1 <- data[[1]]
  qc1 <- data[[2]]
  subject.info1 <- data[[3]]
  qc.info1 <- data[[4]]

  s.number.all <- NULL
  s.per.all <- NULL
  q.number.all <- NULL
  q.per.all <- NULL

  b <- list()
  for (i in 1:length(subject1)) {
    temp.subject <- subject1[[i]]
    temp.qc <- qc1[[i]]

    if (what == "mv") {
      s.number <- sum(is.na(temp.subject))
      q.number <- sum(is.na(temp.qc))
    } else {
      s.number <- sum(temp.subject == 0)
      q.number <- sum(temp.qc == 0)
    }

    s.per <-
      round(s.number * 100 / (nrow(temp.subject) * ncol(temp.subject)), 2)
    q.per <-
      round(q.number * 100 / (nrow(temp.qc) * ncol(temp.qc)), 2)

    ## record
    s.number.all[i] <- s.number
    s.per.all[i] <- s.per

    q.number.all[i] <- q.number
    q.per.all[i] <- q.per

    #ratio in each feature and sample
    if (what == "mv") {
      s.feature.per <-
        apply(subject1[[i]], 1, function(x) {
          sum(is.na(x)) * 100 / ncol(subject1[[i]])
        })
      s.sample.per <-
        apply(subject1[[i]], 2, function(x) {
          sum(is.na(x)) * 100 / nrow(subject1[[i]])
        })

      q.feature.per <-
        apply(qc1[[i]], 1, function(x) {
          sum(is.na(x)) * 100 / ncol(qc1[[i]])
        })
      q.sample.per <-
        apply(qc1[[i]], 2, function(x) {
          sum(is.na(x)) * 100 / nrow(qc1[[i]])
        })
    } else {
      s.feature.per <-
        apply(subject1[[i]], 1, function(x) {
          sum(x == 0) * 100 / ncol(subject1[[i]])
        })
      s.sample.per <-
        apply(subject1[[i]], 2, function(x) {
          sum(x == 0) * 100 / nrow(subject1[[i]])
        })

      q.feature.per <-
        apply(qc1[[i]], 1, function(x) {
          sum(x == 0) * 100 / ncol(qc1[[i]])
        })
      q.sample.per <-
        apply(qc1[[i]], 2, function(x) {
          sum(x == 0) * 100 / nrow(qc1[[i]])
        })
    }

    ## plot
    pdf(file.path(
      path,
      paste(
        "Batch",
        i,
        "feature",
        ifelse(what == "mv", "MV", "zero") ,
        "distribution.pdf"
      )
    ), width = 14)
    layout(matrix(c(1:2), ncol = 2))
    par(mar = c(5, 5, 4, 2))
    #subject
    plot(
      s.feature.per,
      xlab = "Feature index",
      ylab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      main = paste(
        ifelse(what == "mv", "MV", "zero"),
        "ratio distribution (Subject)"
      ),
      cex.main = 1.3
    )
    abline(
      h = 30,
      lty = 2,
      col = "firebrick1",
      lwd = 2
    )
    #QC
    plot(
      q.feature.per,
      xlab = "Feature index",
      ylab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      main = paste(
        ifelse(what == "mv", "MV", "zero"),
        "ratio distribution (QC)"
      ),
      cex.main = 1.3
    )
    abline(
      h = 30,
      lty = 2,
      col = "firebrick1",
      lwd = 2
    )

    #subject
    hist(
      s.feature.per,
      xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
      cex.lab = 1.3,
      cex.axis = 1.3,
      main = paste(
        "Histogram of",
        ifelse(what == "mv", "MV", "zero"),
        "ratio distribution (Subject)"
      ),
      cex.main = 1.3
    )
    #QC
    hist(
      q.feature.per,
      xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
      cex.lab = 1.3,
      cex.axis = 1.3,
      main = paste(
        "Histogram of",
        ifelse(what == "mv", "MV", "zero"),
        "ratio distribution (QC)"
      ),
      cex.main = 1.3
    )
# browser()
    #subject
    plot(
      x = sort(s.feature.per),
      y = c(1:length(s.feature.per)) * 100 / length(s.feature.per),
      type = "l",
      xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
      ylab = "Cumulative feature percentage (%)",
      cex.lab = 1.3,
      cex.axis = 1.3,
      lwd = 2,
      main = paste(
        "Cumulative",
        ifelse(what == "mv", "MV", "zero"),
        "ratio (Subject)"
      )
    )
    a <-
      round(sum(s.feature.per < 50) * 100 / length(s.feature.per), 2)
    abline(
      v = 50,
      lty = 2,
      col = "firebrick1",
      lwd = 2
    )
    legend(
      "topleft",
      legend = paste(a, "%"),
      title = paste(ifelse(what == "mv", "MV", "zero"), "ratio < 30%"),
      bty = "n"
    )
    #QC
    plot(
      x = sort(q.feature.per),
      y = c(1:length(q.feature.per)) * 100 / length(q.feature.per),
      type = "l",
      xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
      ylab = "Cumulative feature percentage (%)",
      cex.lab = 1.3,
      cex.axis = 1.3,
      lwd = 2,
      main = paste(
        "Cumulative",
        ifelse(what == "mv", "MV", "zero"),
        "ratio (QC)"
      )
    )
    a <-
      round(sum(q.feature.per < 50) * 100 / length(q.feature.per), 2)
    abline(
      v = 50,
      lty = 2,
      col = "firebrick1",
      lwd = 2
    )
    legend(
      "topleft",
      legend = paste(a, "%"),
      title = paste(ifelse(what == "mv", "MV", "zero"), "ratio < 30%"),
      bty = "n"
    )
    dev.off()

    b[[i]] <- data.frame(s.feature.per, q.feature.per)
    colnames(b[[i]]) <-
      c(
        paste(
          "Batch",
          i,
          "Subject",
          ifelse(what == "mv", "MV", "zero"),
          "MV ratio"
        ),
        paste(
          "Batch",
          i,
          "QC",
          ifelse(what == "mv", "MV", "zero"),
          "MV ratio"
        )
      )
  }

  if (length(b) > 1) {
    b1 <- b[[1]]
    for (i in 2:length(b)) {
      b1 <- cbind(b1, b[[i]])
    }
  } else{
    b1 <- b[[1]]
  }
  b1 <- cbind(tags[, "name"], b1)
  colnames(b1)[1] <- "Feature name"
  write.csv(b1,
            file.path(path,
                      paste(
                        "feature",
                        ifelse(what == "mv", "MV", "zero") , "ratio.csv"
                      )))
# browser()
  ## whole data
  if (what == "mv") {
    q.sample.ratio <-
      apply(qc, 2, function(x) {
        sum(is.na(x)) * 100 / nrow(qc)
      })
    s.sample.ratio <-
      apply(subject, 2, function(x) {
        sum(is.na(x)) * 100 / nrow(subject)
      })
  } else {
    q.sample.ratio <-
      apply(qc, 2, function(x) {
        sum(x == 0) * 100 / nrow(qc)
      })
    s.sample.ratio <-
      apply(subject, 2, function(x) {
        sum(x == 0) * 100 / nrow(subject)
      })
  }
  pdf(file.path(path, "Sample MV distribution.pdf"))
  par(mar = c(5, 5, 4, 2))
  plot(
    subject.order,
    s.sample.ratio,
    xlab = "Injection order",
    ylab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
    cex.lab = 1.3,
    cex.axis = 1.3,
    pch = 19,
    main = paste(ifelse(what == "mv", "MV", "zero"), "ratio distribution"),
    cex.main = 1.3,
    col = "royalblue",
    ylim = c(0, max(c(
      s.sample.ratio, q.sample.ratio
    ))),
    xlim = c(1, max(qc.order))
  )

  abline(
    h = 50,
    lty = 2,
    col = "firebrick1",
    lwd = 2
  )

  #add text
  idx <- which(s.sample.ratio >= 50)
  if (length(idx) >= 1) {
    text(
      x = subject.order[idx],
      y = s.sample.ratio[idx],
      labels = colnames(subject)[idx],
      pos = 1
    )
  }
  points(qc.order, q.sample.ratio, pch = 19, col = "firebrick1")

  #add text
  idx <- which(q.sample.ratio >= 50)
  if (length(idx) >= 1) {
    text(
      x = qc.order[idx],
      y = q.sample.ratio[idx],
      labels = colnames(qc)[idx],
      pos = 1
    )
  }
  v.index <- qc.order[which(diff(qc.order) == 1)]
  for (i in v.index) {
    abline(v = i, lty = 2, lwd = 2)
  }
  dev.off()
  c <- data.frame(c(q.sample.ratio, s.sample.ratio))
  colnames(c) <- paste(ifelse(what == "mv", "MV", "zero"), "ratio")
  # browser()
  write.csv(c, file.path(path,
                         paste(
                           "Sample",
                           ifelse(what == "mv", "MV", "zero"),
                           "distribution.csv"
                         )))
  # browser()
  ## MV information in each batch
  batch.info <-
    data.frame(s.number.all,
               s.per.all,
               q.number.all,
               q.per.all)
  colnames(batch.info) <-
    c(
      paste("Subject", ifelse(what == "mv", "MV", "zero"), "number"),
      paste("Subject", ifelse(what == "mv", "MV", "zero"), "percentage(%)"),
      paste("QC", ifelse(what == "mv", "MV", "zero"), "number"),
      paste("QC", ifelse(what == "mv", "MV", "zero"), "percentage(%)")
    )
  write.csv(batch.info,
            file.path(path,
                      paste(
                        ifelse(what == "mv", "MV", "zero"),
                        "information in each batch.csv"
                      )))


}