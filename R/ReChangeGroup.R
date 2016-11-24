#' @title ReChangeGroup
#' @description Change the group information in MetFlowData for statistical
#' analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param new.groupNew Group information name. Default is new.group.csv. It can
#'  be use the sample.information which is only changed the group column.
#' @return Return a standard MetProcesser data which is changed the
#' group informatio.
#' @export

ReChangeGroup <- function(MetFlowData = MetFlowData,
                          new.group = "new.group.csv") {
  # browser()
  new.group <- read.csv(new.group, stringsAsFactors = F)

  subject.info <- MetFlowData[['subject.info']]
  subject <- MetFlowData[['subject']]
  subject.order <- MetFlowData[['subject.order']]

  ##which sample you want to remove from the dataset
  remove.name <- new.group[,1][which(is.na(new.group[,"group"]))]
  if(length(remove.name) != 0) {
  cat("The samples you want to remove from dataset are:\n")
  cat(remove.name)
  right <- readline("Right(y) or wrong(n)?")
  if (right == "n") {
    cat("Please change your new group information again!")
    return(MetFlowData)}
  }

# browser()

  ##remove the NA from new.group inforamtion
  new.group <- new.group[which(!is.na(new.group[,"group"])), ]
  new.subject.info <- new.group[new.group[, "class"] == "Subject", ]
  new.subject.name <- new.subject.info[, 1]

  ##remove samples from MetFlowData
  if (length(remove.name) !=0) {
    remove.idx <- match(remove.name, subject.info[,1])
    subject <- subject[ ,-remove.idx]
    subject.info <- subject.info[-remove.idx,]
    subject.order <- subject.order[-remove.idx]
  }

  subject.name <- subject.info[, 1]

  ##change group information
  index <- match(subject.name, new.subject.name)
  new.subject.info <- new.subject.info[index, ]

  MetFlowData[["subject.info"]] <- new.subject.info
  MetFlowData[["subject"]] <- subject
  MetFlowData[["subject,order"]] <- subject.order
  return(MetFlowData)
}