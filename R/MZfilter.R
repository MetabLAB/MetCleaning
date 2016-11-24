#' @title MZfilter
#' @description Filter feature and samples according to MV/zero ratio.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param var.per.cutoff The MV/zero ratio cutoff value of features,
#' default is 0.5. It means that for a feature, non-MV/zero ratio must be
#' larger than obs.per.cutoff.
#' @param obs.per.cutoff The MV/zero ratio cutoff value of samples,
#' default is 0.5. It means that for a sample, non-MV ratio/zero must be larger
#' than obs.per.cutoff.
#' @param what Filter missing values ("mv") or zero values ("zero")?
#' @param path Work directory.
#' @return Return a MetFlowData which has been filtered according MV/zero ratio.
#' @export

MZfilter <- function(MetFlowData = MetFlowData,
                     obs.per.cutoff = 0.5,
                     var.per.cutoff = 0.5,
                     what = "mv",
                     path = NULL) {
  if (is.null(path)) {
    path <- getwd()
  }
  else {
    dir.create(path)
  }
  options(warn = -1)
  # browser()
  MZfinderData <- MZfinder(
    MetFlowData = MetFlowData,
    obs.per.cutoff = obs.per.cutoff,
    var.per.cutoff = var.per.cutoff,
    what = what,
    path = path
  )

  MetFlowData <- MZfinderData[[1]]
  feature.remove <- MZfinderData[[2]]
  qc.remove <- MZfinderData[[3]]
  subject.remove <- MZfinderData[[4]]

  qc <- MetFlowData[["qc"]]
  subject <- MetFlowData[["subject"]]
  qc.info <- MetFlowData[["qc.info"]]
  subject.info <- MetFlowData[["subject.info"]]
  tags <- MetFlowData[["tags"]]
  hasQC <- MetFlowData[["hasQC"]]

  if (length(feature.remove) != 0) {
    if (hasQC == "yes")  {
      qc <- qc[-feature.remove, ]
    }
    subject <- subject[-feature.remove, ]
    tags <- tags[-feature.remove, ]
  }
# browser()
  if (!is.null(qc.remove) & length(qc.remove) != 0) {
    cat("QC shoulde be removed are: \n")
    cat(qc.remove)
    cat("\n")
    cat("\n")
    qc.remove <-
      readline(
        "Which QC you want to remove(please type the index of QC sample,
         and separate them using comma,
         if you don't want to remove any QC, please type n):"
      )
    # browser()
    if (qc.remove == "n") {
      qc <- qc
    } else {
      qc.remove <- strsplit(qc.remove, split = ",")
      qc.remove <- as.numeric(qc.remove[[1]])
      qc.remove <- as.numeric(qc.remove)
      qc.remove.name <- colnames(qc)[qc.remove]
      qc <- qc[, -qc.remove]
    }
  }

  if (!is.null(subject.remove) & length(subject.remove) != 0) {
    cat("Subject shoulde be removed are: \n")
    cat(subject.remove)
    cat("\n")
    subject.remove <-
      readline(
        "Which subject you want to remove(please type the index of subject sample,
        and separate them using comma,
        if you don't want to remove any subject, please type n):"
      )
    if (subject.remove == "n") {
      subject <- subject
    }
    else {
      subject.remove <- strsplit(subject.remove, split = ",")
      subject.remove <- as.numeric(subject.remove[[1]])
      subject.remove <-
        as.numeric(subject.remove)
      subject <- subject[, -subject.remove]
    }
  }


  subject.name <- colnames(subject)
  if (hasQC == 'yes') {
    qc.name <- colnames(qc)
  }
  else {
    qc.name <- NULL
  }
# browser()
  subject.index <-
    which(is.na(match(subject.info[, 1], subject.name)))
  if (length(subject.index) != 0) {
    subject.info <- subject.info[-subject.index,]
    MetFlowData[["subject.order"]] <-
      MetFlowData[["subject.order"]][-subject.index]
  }

  if (hasQC == 'yes') {
    qc.index <- which(is.na(match(qc.info[, 1], qc.name)))
    if (length(qc.index) != 0) {
      qc.info <- qc.info[-qc.index,]
      MetFlowData[["qc.order"]] <-
        MetFlowData[["qc.order"]][-qc.index]
    }
  }

  MetFlowData[["subject.info"]] <- subject.info
  MetFlowData[["qc.info"]] <- qc.info
  MetFlowData[["subject"]] <- subject
  MetFlowData[["qc"]] <- qc
  MetFlowData[["tags"]] <- tags
  MetFlowData[["mv.filter"]] <- "yes"
  MetFlowData[["mv.filter.criteria"]] <-
    paste("varibale:",
          var.per.cutoff,
          "observation:",
          obs.per.cutoff)
  options(warn = 0)
  return(MetFlowData)
  }
