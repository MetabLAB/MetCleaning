#' @title PrintMethod
#' @description Print method for class in MetProcesser.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param myclass The name of class object in MetProcesser.
#' @return Print some information of class object in your screen.
#' @export
#' @seealso \code{\link{print}}

## MetFlowData
print.MetFlowData <- function(MetFlowData) {
  hasQC <- MetFlowData[["hasQC"]]
  subject <- MetFlowData[["subject"]]
  qc <- MetFlowData[["qc"]]
  tags <- MetFlowData[["tags"]]
  subject.info <- MetFlowData[["subject.info"]]
  qc.info <- MetFlowData[["qc.info"]]
  subject.order <- MetFlowData[["subject.order"]]
  qc.order <- MetFlowData[["qc.order"]]

  subject.batch <- as.numeric(subject.info[,4])
  subject.name <- subject.info[,1]

  if(hasQC == "no") {qc.batch <- NULL; qc.name <- NULL}
  else {
  qc.batch <- as.numeric(qc.info[,4])
  qc.name <- qc.info[,1]
  }

  unique.batch <- unique(subject.batch)
  summary.subject.batch <- table(subject.batch)
  if(hasQC == "no") {summary.qc.batch <- NULL}
  else {
  summary.qc.batch <- table(qc.batch)
  }
  summary.batch <- cbind(summary.subject.batch, summary.qc.batch)

  ## batch information
  cat("There are", length(unique.batch), ifelse(length(unique.batch)==1, "batch", "batches"))
  cat("\n")

  for (i in 1:length(unique.batch)){
    cat(paste("Batch",i,":\n"))
    cat(paste("Subject sample number:",summary.batch[i,1]))
    cat("\n")
    if (hasQC == "no") {cat(paste("QC sample number:",0))}
    else {
    cat(paste("QC sample number:",summary.batch[i,2]))
    }
    cat("\n")
  }

  ## peak information
  cat("---------------\n")
  cat(paste("Peak number:",nrow(subject)))
  cat("\n")
  cat("The tags information contains:\n")
  cat(colnames(tags))
  cat('\n')

  ## subject sample information
  cat("---------------\n")
  cat("Subject sample info:\n")
  cat(paste("Subject sample number:", ncol(subject)))
  cat("\n")
  cat("The Subject information contains:\n")
  cat(colnames(subject.info))
  cat('\n')

  ## QC sample information
  if (hasQC != "no") {
  cat("---------------\n")
  cat("QC sample info:\n")
  cat(paste("QC sample number:", ncol(qc)))
  cat("\n")
  cat("---------------\n")
  }

  ## other processing information
  cat("MV imputation:",MetFlowData[["mv.imputation"]])
  cat("\n")
  cat("Imputation method:",MetFlowData[["imputation.method"]])
  cat("\n")
  cat("Zero filter:",MetFlowData[["zero.filter"]])
  cat("\n")
  cat("Zero filter criteria:",MetFlowData[["zero.filter.criteria"]])
  cat("\n")
  cat("QC outlier filter:",MetFlowData[["qc.outlier.filter"]])
  cat("\n")
  cat("Normalization:",MetFlowData[["normalization"]])
  cat("\n")
  cat("Normalization method:",MetFlowData[["normalization.method"]])
  cat("\n")
  cat("Data integration:",MetFlowData[["data.integration"]])
  cat("\n")
  cat("Data integration method:",MetFlowData[["data.integration.method"]])
  cat("\n")
  cat("Has IS:",MetFlowData[["hasIS"]])
  cat("\n")
  cat("Has QC:",MetFlowData[["hasQC"]])
  cat("\n")
  cat("Peak identification:",MetFlowData[["peak.identification"]])
  cat("\n")
}



###print method for SXTMinfracData
print.SXTMinifracData <- function(SXTMinifracData) {
  data <- SXTMinifracData[["data"]]
  var.index <- SXTMinifracData[["var.index"]]
  obs.index <- SXTMinifracData[["obs.index"]]
  filter.item <- SXTMinifracData[["filter.item"]]
  filter.rule <- SXTMinifracData[["filter.rule"]]
  minifrac.variable <- SXTMinifracData[["minifrac.variable"]]
  minifrac.observation <- SXTMinifracData[["minifrac.observation"]]

  ##filter information
  cat("filter.item:",filter.item)
  cat("\n")
  cat("filter.rule:",filter.rule)
  cat("\n")
  cat("minifrac.variable:",minifrac.variable)
  cat("\n")
  cat("minifrac.observation:",minifrac.observation)
  cat("\n")
  cat("\n")

  cat("There are",length(data), "data.\n")
  cat("-------------------------------\n")

  for (i in 1:length(data)) {
    temp.data <- data[[i]]
    cat("Data",i)
    cat("\n")
    cat(nrow(temp.data),"rows\n")
    cat(ncol(temp.data),"columns\n")
    cat("\n")
    }

}

##print method for SXTpcaData
print.SXTpcaData <- function(SXTpcaData) {
info <- SXTpcaData[["info"]]
QC <- SXTpcaData[["QC"]]
scale.method <- SXTpcaData[["scale.method"]]
  cat("There are",length(info), "class in this PCA analysis:\n")
  cat(names(info))
  cat("\n")
  cat("---------------------------------------\n")
  cat("QC are contained in this PCA analysis?",QC)
}

