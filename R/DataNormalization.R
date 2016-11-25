#' @title DataNormalization
#' @description Normalize dataset.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param path Workdirectory.
#' @param method: Normalization method, mean, median, total svr or loess,
#' default is svr. Please see the details.
#' @param dimension1: Remain dimension or not. Default is TRUE.
#' @param svrloess You can see MetNormalizer for help.
#' @return The normalization results can be got from help of
#' \code{\link{MetNormalizer}}.
#' @seealso \code{\link{MetNormalizer}}
#' @export
#' @details Data normalization is a useful method to reduce unwanted variances
#' in metabolomics data. The most used normalization methods is
#' \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:2303e537-2864-4df5-91f0-c3ae731711af}{sample scale method}
#' (median, mean and total intensity). In large
#' scale metabolomics, QC based normalization methods are more useful.
#' \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:abe41368-d08d-4806-871f-3aa035d21743}{LOESS}
#' and \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:2303e537-2864-4df5-91f0-c3ae731711af}{SVR}
#' normalization are two most used normalization
#' methods.


### Data normalization for MetFlowData
DataNormalization <- function(MetFlowData = MetFlowData,
                              path = NULL,
                              method = "svr",
                              dimension1 = TRUE,
                              ## parameters for loess normalization
                              optimization = TRUE,
                              begin = 0.5,
                              end = 1,
                              step = 0.2,
                              ##svr parameters
                              multiple = 5,
                              threads = 2,
                              rerun.loess = TRUE,
                              rerun.svr = TRUE,
                              peakplot = TRUE,
                              datastyle = "tof",
                              user = "other") {
  options(warn = -1)
  if (is.null(path)) {
    path <- getwd()
  }
  # browser()
  path1 <- file.path(path, "7Normalization result")
  dir.create(path1)

  qc <- MetFlowData[["qc"]]
  subject <- MetFlowData[["subject"]]
  tags <- MetFlowData[["tags"]]
  subject.info <- MetFlowData[["subject.info"]]
  qc.info <- MetFlowData[["qc.info"]]

  if (sum(is.na(MetFlowData[["qc"]]))+sum(is.na(MetFlowData[["subject"]])) != 0)
  {
    stop("Plase impute MV in sampe first!!!")
  }

  ##mean normalization
  if (method == "mean") {
    qc1 <- apply(qc, 2, function(x) {
      x / mean(x)
    })
    subject1 <- apply(subject, 2, function(x) {
      x / mean(x)
    })
    if (dimension1)
    {
      qc.mean <- apply(qc, 1, mean)
      qc2 <- qc1 * qc.mean
      subject2 <- subject1 * qc.mean
      MetFlowData[["qc"]] <- qc2
      MetFlowData[["subject"]] <- subject2
    }
    else {
      MetFlowData[["qc"]] <- qc1
      MetFlowData[["subject"]] <- subject1
    }
  }

  #median normalization
  if (method == "median") {
    qc1 <- apply(qc, 2, function(x) {
      x / median(x)
    })
    subject1 <- apply(subject, 2, function(x) {
      x / median(x)
    })
    if (dimension1)
    {
      qc.median <- apply(qc, 1, median)
      qc2 <- qc1 * qc.median
      subject2 <- subject1 * qc.median
      MetFlowData[["qc"]] <- qc2
      MetFlowData[["subject"]] <- subject2
    }
    else {
      MetFlowData[["qc"]] <- qc1
      MetFlowData[["subject"]] <- subject1
    }
  }

  ##total intensity normalization
  if (method == "total") {
    qc1 <- apply(qc, 2, function(x) {
      x / sum(x)
    })
    subject1 <- apply(subject, 2, function(x) {
      x / sum(x)
    })
    if (dimension1)
    {
      qc.mean <- apply(qc, 1, mean)
      qc2 <- qc1 * qc.mean
      subject2 <- subject1 * qc.mean
      MetFlowData[["qc"]] <- qc2
      MetFlowData[["subject"]] <- subject2
    }
    else {
      MetFlowData[["qc"]] <- qc1
      MetFlowData[["subject"]] <- subject1
    }
  }
  # browser()
  ## split batch

  data <- SplitBatch(MetFlowData = MetFlowData)

  subject1 <- data[[1]]
  qc1 <- data[[2]]
  subject.info1 <- data[[3]]
  qc.info1 <- data[[4]]

  ##SVR normalization
  if (method == "svr") {
    for (i in 1:length(qc1)) {
      cat(paste("Batch", i))
      cat("\n")
      cat("------------------\n")
      MetFlowData[["data.normalization"]] <- "yes"
      tags <- MetFlowData[["tags"]]
      data <- cbind(tags, qc1[[i]], subject1[[i]])
      sample.info <- rbind(subject.info1[[i]], qc.info1[[i]])

      path2 <- file.path(path1, paste("Batch", i, "normalization"))
      dir.create(path2)
      write.csv(data, file.path(path2, "data.csv"), row.names = F)
      write.csv(sample.info,
                file.path(path2, "worklist.csv"),
                row.names = F)
      MetNormalizer(
        normalization.method = "svr",
        peakplot = peakplot,
        multiple = multiple,
        rerun.svr = rerun.svr,
        datastyle = datastyle,
        dimension1 = dimension1,
        user = user,
        path = path2
      )
    }
    # browser()
    ## replace subject and qc
    for (i in 1:length(qc1)) {
      path.for.data <- file.path(path1, paste("Batch", i, "normalization"))
      load(file.path(path.for.data, "svr normalization result/data svr nor"))
      subject1[[i]] <- t(sample.nor)
      qc1[[i]] <- t(QC.nor)
    }

    subject2 <- subject1[[1]]
    qc2 <- qc1[[1]]

    if (length(qc1) > 1) {
      for (i in 2:length(subject1)) {
        subject2 <- cbind(subject2, subject1[[i]])
        qc2 <- cbind(qc2, qc1[[i]])
      }
    }
    MetFlowData[["subject"]] <- subject2
    MetFlowData[["qc"]] <- qc2
  }

  ## LOESS normalization
  if (method == "loess") {
    for (i in 1:length(qc1)) {
      cat(paste("Batch", i))
      cat("\n")
      cat("------------------\n")
      tags <- MetFlowData[["tags"]]
      data <- cbind(tags, subject1[[i]], qc1[[i]])
      sample.info <- rbind(subject.info1[[i]], qc.info1[[i]])
      path2 <- file.path(path1, paste("Batch", i, "normalization"))
      dir.create(path2)
      write.csv(data, file.path(path2, "data.csv"), row.names = F)
      write.csv(sample.info,
                file.path(path2, "worklist.csv"),
                row.names = F)
      MetNormalizer(
        normalization.method = "loess",
        peakplot = peakplot,
        begin = begin,
        end = end,
        step = step,
        rerun.loess = return.loess,
        dimension1 = dimension1,
        datastyle = datastyle,
        user = user,
        path = path2
      )
      load(file.path(path2, "loess normalization result/data loess nor"))
      qc1[[i]] <- t(QC.nor)
      subject1[[i]] <- t(sample.nor)
    }
    subject2 <- subject1[[1]]
    qc2 <- qc1[[1]]

    if (length(qc1) > 1) {
      for (i in 2:length(subject1)) {
        subject2 <- cbind(subject2, subject1[[i]])
        qc2 <- cbind(qc2, qc1[[i]])
      }
    }
    MetFlowData[["subject"]] <- subject2
    MetFlowData[["qc"]] <- qc2
  }
  MetFlowData[["normalization"]] <- "yes"
  MetFlowData[["normalization.method"]] <- method
  options(warn = 0)
  return(MetFlowData)

}
