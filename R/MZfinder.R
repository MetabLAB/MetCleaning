#' Find feature and sample outliers according to MV/zero ratio.
#'
#' @title MZfinder
#' @description Find feature and sample outliers according to MV/zero ratio.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param obs.per.cutoff The observation MV/zero ratio cutoff.
#' @param var.per.cutoff The variable MV/zero ratio cutoff.
#' @return Return a MZfinderData contains MetFlowData,feature.remove, qc.remove and subject.remove.
#' @export

##remove peaks whose MV ratio > threshold in QC or subject.
MZfinder <- function(MetFlowData = MetFlowData,
                     obs.per.cutoff = 0.5,
                     var.per.cutoff = 0.5,
                     what = "mv",
                     path = NULL) {
  # browser()
  options(warn = -1)
  if (is.null(path))
    path <- getwd()
  hasQC <- MetFlowData[["hasQC"]]
  if (what == "mv") {
    path1 <- file.path(path, "mv finder")
  } else {
    path1 <- file.path(path, "zero finder")
  }
  dir.create(path1)

  qc <- MetFlowData[["qc"]]
  subject <- MetFlowData[["subject"]]
  qc.info <- MetFlowData[["qc.info"]]
  subject.info <- MetFlowData[["subject.info"]]
  tags <- MetFlowData[["tags"]]
  subject.order <- as.numeric(MetFlowData[["subject.order"]])
  qc.order <- as.numeric(MetFlowData[["qc.order"]])

  subject.name <- subject.info[, 1]
  if (hasQC != "no") {
    qc.name <- qc.info[, 1]
  }
  else {
    qc.name <- NULL
  }

  data <- SplitBatch(MetFlowData = MetFlowData)
  subject1 <- data[[1]]
  qc1 <- data[[2]]
  subject.info1 <- data[[3]]
  qc.info1 <- data[[4]]

  var.index <- list()
  if (hasQC == "yes") {
    ##remove peak whose MV/zero ratio more than 50%
    for (i in 1:length(subject1)) {
      if (what == "mv") {
        SXTMinifracData <- SXTMinifrac(
          subject1[[i]],
          qc1[[i]],
          filter.item = "mv",
          filter.rule = "union",
          minifrac.variable = var.per.cutoff,
          minifrac.observation = 0
        )
      } else {
        SXTMinifracData <- SXTMinifrac(
          subject1[[i]],
          qc1[[i]],
          filter.item = "zero",
          filter.rule = "union",
          minifrac.variable = var.per.cutoff,
          minifrac.observation = 0
        )
      }

      var.index[[i]] <- SXTMinifracData[["var.index"]]
    }
  }
  else {
    for (i in 1:length(subject1)) {
      if (what == "mv") {
        SXTMinifracData <- SXTMinifrac(
          subject1[[i]],
          filter.item = "mv",
          minifrac.variable = var.per.cutoff,
          minifrac.observation = 0
        )
      } else {
        SXTMinifracData <- SXTMinifrac(
          subject1[[i]],
          filter.item = "zero",
          minifrac.variable = var.per.cutoff,
          minifrac.observation = 0
        )
      }
      var.index[[i]] <- SXTMinifracData[["var.index"]]
    }
  }

  var.index1 <- var.index[[1]]
  if (length(var.index) > 1) {
    for (i in 2:length(var.index)) {
      var.index1 <- intersect(var.index1, var.index[[i]])
    }
  }

  subject1 <- lapply(subject1, function(x) {
    x[var.index1,]
  })

  if (hasQC == "yes") {
    qc1 <- lapply(qc1, function(x) {
      x[var.index1,]
    })
  }
  else {
    qc1 <- NULL
  }

  if (hasQC == "yes") {
    subject2 <- subject1[[1]]
    qc2 <- qc1[[1]]

    if (length(var.index) > 1) {
      for (i in 2:length(subject1)) {
        subject2 <- cbind(subject2, subject1[[i]])
        qc2 <- cbind(qc2, qc1[[i]])
      }
    }
  }
  else {
    subject2 <- subject1[[1]]
    if (length(var.index) > 1) {
      for (i in 2:length(subject1)) {
        subject2 <- cbind(subject2, subject1[[i]])
      }
    }
    qc <- NULL
  }

  ##remove sample whose MV/zero more than 50%
  if (hasQC == "yes") {
    q.sample.per <- apply(qc2, 2, function(x) {
      ifelse(what == "mv", sum(!is.na(x)), sum(x != 0)) / nrow(qc2)
    })
    qc.index <- which(q.sample.per > obs.per.cutoff)
  }
  else {
    qc.index <- NULL
  }
  s.sample.per <-
    apply(subject2, 2, function(x) {
      ifelse(what == "mv", sum(!is.na(x)), sum(x != 0)) / nrow(subject2)
    })
  subject.index <- which(s.sample.per > obs.per.cutoff)

  ## output some information about feature and samples
  feature.remove <- setdiff(c(1:nrow(subject)), var.index1)
  if (hasQC == "yes") {
    qc.remove <- setdiff(c(1:ncol(qc)), qc.index)
  }
  else {
    qc.remove <- NULL
  }
  subject.remove <- setdiff(c(1:ncol(subject)), subject.index)

  if (length(feature.remove) == 0) {
    cat("No features should be removed.\n")
  }
  else {
    write.csv(tags[feature.remove,],
              file.path(path1, "Feature to remove information.csv"))
  }

  if (hasQC != "no") {
    if (length(qc.remove) == 0) {
      cat("No QC should be removed.\n")
    }
    else {
      cat(colnames(qc2)[qc.remove], "sholud be removed!!!\n")
    }
  }
  # browser()
  if (length(subject.remove) == 0) {
    cat("No subject should be removed.\n")
  }
  else {
    cat(colnames(subject2)[subject.remove], "sholud be removed!!!\n")
  }
  options(warn = 0)
  MZfinderData <- list(
    MetFlowData = MetFlowData,
    feature.remove = feature.remove,
    qc.remove = qc.remove,
    subject.remove = subject.remove
  )
  class(MZfinderData) <- "MZfinderData"
  return(MZfinderData)
}