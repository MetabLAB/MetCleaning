#' @title ExportData
#' @description Export MetFlowData as csv.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param data.name The name of the data you want to output. Default is
#' "data_new".
#' @param subject.info.name The name for subject information you want to
#' output. Default is "subject.info".
#' @param qc.info.name The name for QC information you want to output.
#' Default is "qc.info".
#' @param path Work directory.
#' @return Write csv data.
#' @seealso \code{\link{ImportData}}
#' @export

ExportData <- function(MetFlowData = MetFlowData,
                       data.name = "data_new",
                       subject.info.name = "subject.info",
                       qc.info.name = "qc.info",
                       path = NULL){
  if (is.null(path)) {path <- getwd()}
  else {dir.create(path)}

  subject <- MetFlowData[["subject"]]
  qc <- MetFlowData[["qc"]]
  tags <- MetFlowData[["tags"]]
  subject.info <- MetFlowData[["subject.info"]]
  qc.info <- MetFlowData[["qc.info"]]

  write.csv(cbind(tags, subject, qc), file.path(path,paste(data.name,".csv", sep = "")), row.names = F)
  write.csv(subject.info, file.path(path,paste(subject.info.name,".csv", sep = "")), row.names = F)
  write.csv(qc.info, file.path(path,paste(qc.info.name,".csv", sep = "")), row.names = F)
}