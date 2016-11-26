#' @title SubjectOutlierFilter
#' @description Using PCA to filter subject outliers.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @return MetFlowData whose subject outliers have been removed.
#' @export
#' @seealso \code{\link{QCOutlierFilter}}

SubjectOutlierFilter <- function(MetFlowData = MetFlowData,
                                 CI = 0.95,
                                 path = NULL) {
  # browser()
  SubjectOutlierFiderData <- SubjectOutlierFinder(MetFlowData = MetFlowData,
                                        CI = CI,
                                        path = path)

  metData <- SubjectOutlierFiderData[[1]]
  obs.remove <- SubjectOutlierFiderData[[2]]

  data <- SplitBatch(MetFlowData = metData)
  subject1 <- data[[1]]

  for (i in 1:length(subject1)) {
    cat(paste("Batch",i))
    cat("\n")
    cat("-------------------------------------------\n")
    temp.subject <- subject1[[i]]
    temp.idx <- obs.remove[[i]]
    if (length(temp.idx) != 0) {
      # browser()
       cat("Subject shoulde be removed are:")
      cat(temp.idx)
      cat("\n")
      temp.idx <-
        readline(
          "Which subject you want to remove(please type the index of subject sample,
          and separate them using comma,
          if you don't want to remove any subject, please type n):"
        )

      if (temp.idx == "n") {
        temp.subject <- temp.subject
      } else {
        temp.idx <- strsplit(temp.idx, split = ",")
        temp.idx <- as.numeric(temp.idx[[1]])
        temp.idx <- as.numeric(temp.idx)
        temp.subject <- temp.subject[, -temp.idx]
      }
    } else {
      temp.subject <- temp.subject
    }
    subject1[[i]] <- temp.subject
  }


  subject2 <- subject1[[1]]
  if (length(subject1) > 1) {
    for (i in 2:length(subject1)) {
      subject2 <- cbind(subject2, subject1[[i]])
    }
  }

  ##remove subject information who have been removed from data
  subject.name <- colnames(subject2)
  subject.info <- metData[["subject.info"]]
  subject.order <- metData[["subject.order"]]
  subject.index <- which(is.na(match(subject.info[, 1], subject.name)))

  if (length(subject.index) != 0) {
    subject.info <- subject.info[-subject.index, ]
    subject.order <- subject.order[-subject.index]
  }

  metData[["subject.info"]] <- subject.info
  metData[["subject"]] <- subject2
  metData[["subject.order"]] <- subject.order
  metData[["subject.outlier.filter"]] <- "yes"
  return(metData)

  }
