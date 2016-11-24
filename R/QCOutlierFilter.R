#' @title QCOutlierFilter
#' @description Using PCA to find QC outliers.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param CI confidence interval.
#' @param path Work directory.
#' @return MetFlowData whose qc outliers have been removed.
#' @export
#' @seealso \code{\link{SubjectOutlierFilter}}

QCOutlierFilter <- function(MetFlowData = MetFlowData,
                            CI = 0.95,
                            path = NULL) {
  # browser()
  QCOutlierFiderData <- QCOutlierFinder(MetFlowData = MetFlowData,
                                        CI = CI,
                                        path = path)

  metData <- QCOutlierFiderData[[1]]
  obs.remove <- QCOutlierFiderData[[2]]

  data <- SplitBatch(MetFlowData = metData)
  qc1 <- data[[2]]

  for (i in 1:length(qc1)) {
    cat(paste("Batch", i))
    cat("\n")
    cat("-----------------------\n")
    temp.qc <- qc1[[i]]
    temp.idx <- obs.remove[[i]]
    if (length(temp.idx) != 0) {
      cat("QC shoulde be removed are:\n")
      cat(temp.idx)
      cat("\n")
      temp.idx <-
        readline(
          "Which QC you want to remove(please type the index of QC sample,
          and separate them using comma,
          if you don't want to remove any QC, please type n):"
        )
      # browser()
      if (temp.idx == "n") {
        temp.qc <- temp.qc
      } else {
        temp.idx <- strsplit(temp.idx, split = ",")
        temp.idx <- as.numeric(temp.idx)
        temp.idx <- as.numeric(temp.idx)
        temp.qc <- temp.qc[, -temp.idx]
      }
    } else {
      temp.qc <- temp.qc
    }
    qc1[[i]] <- temp.qc
  }


  qc2 <- qc1[[1]]
  if (length(qc1) > 1) {
    for (i in 2:length(qc1)) {
      qc2 <- cbind(qc2, qc1[[i]])
    }
  }

  ##remove QC information who have been removed from data
  qc.name <- colnames(qc2)
  qc.info <- metData[["qc.info"]]
  qc.index <- which(is.na(match(qc.info[, 1], qc.name)))

  if (length(qc.index) != 0) {
    qc.info <- qc.info[-qc.index, ]
  }

  metData[["qc.info"]] <- qc.info
  metData[["qc"]] <- qc2
  metData[["qc.outlier.filter"]] <- "yes"
  return(metData)

  }
