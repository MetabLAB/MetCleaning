#' @title SplitBatch
#' @description Split MetFlowData accoding to different batch.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @return Return a data (list), subject, qc, subject.info and qc.info.
#' @export

SplitBatch <- function(MetFlowData = MetFlowData) {
  ## split batch
  # browser()
  hasQC <- MetFlowData[["hasQC"]]
  qc <- MetFlowData[["qc"]]
  subject <- MetFlowData[["subject"]]
  qc.info <- MetFlowData[["qc.info"]]
  subject.info <- MetFlowData[["subject.info"]]
  tags <- MetFlowData[["tags"]]

  subject.name <- as.character(subject.info[, 1])
  if (hasQC != "no") {
    qc.name <- as.character(qc.info[, 1])
  }
  else {
    qc.name <- NULL
  }

  subject.batch <- as.numeric(subject.info[, 4])
  if (hasQC != "no") {
    qc.batch <- as.numeric(qc.info[, 4])
  }
  else {
    qc.batch <- NULL
  }

  subject1 <- list()
  qc1 <- list()
  subject.info1 <- list()
  qc.info1 <- list()
  # browser()
  for (i in 1:length(unique(subject.batch))) {
    subject.name.for.this.batch <- subject.name[subject.batch == i]
    if (hasQC == "no") {
      qc.name.for.this.batch <- NULL
    }
    else {
      qc.name.for.this.batch <- qc.name[qc.batch == i]
    }
    subject.index.for.this.batch <-
      match(subject.name.for.this.batch, colnames(subject))
    subject.index.for.this.batch <-
      subject.index.for.this.batch[!is.na(subject.index.for.this.batch)]
    if (hasQC == "no") {
      qc.index.for.this.batch <- NULL
    }
    else {
      qc.index.for.this.batch <-
        match(qc.name.for.this.batch, colnames(qc))
      qc.index.for.this.batch <-
        qc.index.for.this.batch[!is.na(qc.index.for.this.batch)]
    }

    subject1[[i]] <- subject[, subject.index.for.this.batch]
    subject.info1[[i]] <- subject.info[subject.batch == i, ]
    if (hasQC == "no") {
      qc1[[i]] <- NULL
      qc.info1[[i]] <- NULL
    }
    else {
      qc1[[i]] <- qc[, qc.index.for.this.batch]
      qc.info1[[i]] <- qc.info[qc.batch == i, ]
    }
  }

  data <- list(subject1, qc1, subject.info1, qc.info1)
  return(data)
}