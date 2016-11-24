#' @title SubjectOutlierFinder
#' @description Using PCA to find subject outliers.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param CI Confidence interva.
#' @param path Work directory.
#' @return SubjectOutlierFinderData contains MetFlowData and subject outlier index (obs.remove).
#' @export
#' @seealso \code{\link{SubjectOutlierFilter}}

#subject outlier filtering according to zero ratio and PCA
SubjectOutlierFinder <- function(MetFlowData = MetFlowData,
                                 CI = 0.95,
                                 path = NULL){
  options(warn = -1)
  # browser()
  if (is.null(path)) {
    path <- getwd()
  }else{
    dir.create(path)
  }
  subject <- MetFlowData[["subject"]]
  subject.info <- MetFlowData[["subject.info"]]
  tags <- MetFlowData[["tags"]]

data <- SplitBatch(MetFlowData = MetFlowData)
subject1 <- data[[1]]

obs.remove <- list()
  ## PCA analysis
  for (i in 1:length(subject1)) {
    info <- list("Subject" = colnames(subject1[[i]]))
    SXTpcaData <- SXTpca(subject = subject1[[i]], info = info, QC = F, scale.method = "auto")
    index2 <- SXTpcaFindOutlier(SXTpcaData = SXTpcaData,
                                CI = CI,
                                plot.name = paste("Batch",i,"outliers"),
                                output.plot = TRUE,
                                path = path)
    obs.remove[[i]] <- index2
  }
SubjectOutlierFinderData <- list(MetFlowData = MetFlowData,
                            obs.remove = obs.remove)
class(SubjectOutlierFinderData) <- "SubjectOutlierFinderData"
options(warn = 0)
return(SubjectOutlierFinderData)
}
