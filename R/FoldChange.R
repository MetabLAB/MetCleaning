#' @title FoldChange
#' @description Calculate fold change for MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param to Which the ratio you want? default is case/control.
#' @param ratio Which ratio you want to use to calculate fold change. median ot mean.
#' @return MetFlowData which has been added fold change information in tags.
#' @export

FoldChange <- function(MetFlowData = MetFlowData,
                       to = c("case", "control"),
                       ratio = "median") {
  # browser()
  subject <- MetFlowData[["subject"]]
  subject.info <- MetFlowData[["subject.info"]]
  group <- subject.info[, "group"]
  tags <- MetFlowData[["tags"]]

  if (length(unique(group)) != 2) {
    stop("No two class data!!!")
  }
  group2.index <- which(group == to[1])
  group1.index <- which(group == to[2])

  if (ratio == "median") {
    foldchange <-
      apply(subject, 1, function(x) {
        median(x[group2.index] + 1) / median(x[group1.index] + 1)
      })
  }
  if (ratio == "mean") {
    foldchange <-
      apply(subject, 1, function(x) {
        mean(x[group2.index] + 1) / mean(x[group1.index] + 1)
      })
  }

  if (any(colnames(tags) == "foldchange")) {
    tags[, "foldchange"] <- foldchange
  }
  else {
    tags <- cbind(tags, foldchange)
  }

  MetFlowData[["tags"]] <- tags
  MetFlowData[["foldchange"]] <-
    paste(to[1], " vs ", to[2], "(", ratio, ")", sep = "")
  return(MetFlowData)
}
