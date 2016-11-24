#' Use internal standards to correct data.
#'
#' @title IScorrection
#' @description Use internal standards to correct data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param mzerror mz tolerance.
#' @param rterror rt tolerance.
#' @param IS IS name.
#' @param plot.output Output plot or not?
#' @param path work directory.
#' @param dimensioon1 Keep dimension or not?
#' @return Return a new MetFlowData after IS correction.
#' @export
### IS signal correction for MetFlowData
IScorrction <- function(MetFlowData = MetFlowData,
                        mzerror = 15,
                        rterror = 30,
                        rt.unit.is.second = T,
                        IS = "IS.csv",
                        plot.output = TRUE,
                        path = NULL,
                        dimension1 = TRUE) {
  if (is.null(path)) path <- getwd()

  subject <- MetFlowData[["subject"]]
  qc <- MetFlowData[["qc"]]


  GetIS(MetFlowData = MetFlowData,
        mzerror = mzerror,
        rterror = rterror,
        rt.unit.is.second = rt.unit.is.second,
        IS = IS,
        plot.output = plot.output,
        path = path)
  load(file.path(path,"newIS"))
  newIS <- newIS

  if (all(newIS[,"find.or.not"] == "No")) stop("Don't find IS in your data!!!")
  IS.area <- newIS[,-c(1:8),drop = F]
  IS.area.name <- colnames(IS.area)
  if (nrow(IS.area) > 1) IS.area <- as.numeric(apply(IS.area, 2, function(x) {mean(x)}))
  IS.area <- as.numeric(IS.area)

  qc.index <- match(colnames(qc), IS.area.name)
  IS.area.for.qc <- IS.area[qc.index]
  subject.index <- match(colnames(subject), IS.area.name)
  IS.area.for.subject <- IS.area[subject.index]

  ##IS signal correction
  qc <- t(t(qc)/IS.area.for.qc)
  subject <- t(t(subject)/IS.area.for.subject)

  if (dimension1) {
    IS.mean <- mean(IS.area)
    qc <- qc * IS.mean
    subject <- subject * IS.mean
  }

  MetFlowData[["qc"]] <- qc
  MetFlowData[["subject"]] <- subject
  MetFlowData[["IS.correction"]] <- "yes"

  return(MetFlowData)
}
