#' @title DataOverview
#' @description Give a overview of MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData
#' @param feature.distribution Draw a rt vs mz vs intensity plot or not. Default is TRUE.
#' @return Data overview_RT vs mz vs intensity.pdf:  A RT vs mz vs intensity plot.
#' @return Data overview.txt: A overview information for MetFlowData.
#' @export

### DataOverview for MeeFlowData
DataOverview <- function(MetFlowData = MetFlowData,
                         feature.distribution = TRUE,
                         path = NULL) {
  if (is.null(path)) {
    path <- getwd()
  }
  else{
    dir.create(path)
  }
  # browser()
  hasQC <- MetFlowData[["hasQC"]]
  subject <- MetFlowData[["subject"]]
  qc <- MetFlowData[["qc"]]
  tags <- MetFlowData[["tags"]]
  subject.info <- MetFlowData[["subject.info"]]
  qc.info <- MetFlowData[["qc.info"]]
  subject.order <- as.numeric(MetFlowData[["subject.order"]])
  qc.order <- as.numeric(MetFlowData[["qc.order"]])

  #### a 3D figure: mz vs RT vs intensity
  if (feature.distribution) {
    mz <- as.numeric(tags[, "mz"])
    rt <- as.numeric(tags[, "rt"])
    if (hasQC != "no") {
      sample <- data.frame(qc, subject)
    }
    else {
      sample <- subject
    }

    int.log <-
      log(apply(sample, 1, function(x) {
        mean(x, na.rm = T)
      }) + 10, 10)


    rt.mz.int <- data.frame(rt, mz, int.log)

    library(ggplot2)
    par(mar = c(5, 5, 4, 2))
    rt.mz.int <-
      ggplot(data = rt.mz.int, aes(x = rt, y = mz, colour = int.log)) + geom_point(alpha = 0.3) +
      scale_color_gradient(low = "green", high = "red") +
      labs(x = "Retention time (RT)", y = "Mass to charge ratio (m/z)", colour = "log10(intensity)") +
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14)) +
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10)) +
      theme(legend.title = element_text(size = 14)) +
      theme(legend.text = element_text(size = 10)) +
      ggtitle("RT vs mz vs Intensity") +
      theme_bw()

    ggsave(
      filename = file.path(path, "Data overview_RT vs mz vs intensity.pdf"),
      plot = rt.mz.int,
      width = 8,
      height = 6
    )

  }


  data <- SplitBatch(MetFlowData = MetFlowData)
  subject1 <- data[[1]]
  qc1 <- data[[2]]
  subject.info1 <- data[[3]]
  qc.info1 <- data[[4]]

  file.name.for.txt <- file.path(path, "Data overview.txt")

  ## begin output data overview
  cat("The overview of data\n", file = file.name.for.txt, append = F)
  cat("---------------\n", file = file.name.for.txt, append = T)
  cat('\n', file = file.name.for.txt, append = T)
  cat('\n', file = file.name.for.txt, append = T)

  ## batch information
  cat("There are",
      length(subject1),
      "batches",
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  ##
  subject.number <- unlist(lapply(subject1, length))
  names(subject.number) <- paste("Batch", 1:length(subject1))
  if (hasQC == "no") {
    qc.number <- NULL
  }
  else {
    qc.number <- unlist(lapply(qc1, length))
    names(qc.number) <- paste("Batch", 1:length(subject1))
  }

  cat("Subject number in each Batch:\n",
      file = file.name.for.txt,
      append = T)
  options(warn = -1)
  cat("\n", file = file.name.for.txt, append = T)
  cat(subject.number, file = file.name.for.txt, append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat(qc.number, file = file.name.for.txt, append = T)
  cat("\n", file = file.name.for.txt, append = T)
  options(warn = 0)

  ## peak information
  cat("---------------\n", file = file.name.for.txt, append = T)
  cat(paste("Peak number:", nrow(subject)),
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)

  cat("The tags information contains:\n",
      file = file.name.for.txt,
      append = T)
  cat(colnames(tags), file = file.name.for.txt, append = T)
  cat('\n', file = file.name.for.txt, append = T)
  cat('\n', file = file.name.for.txt, append = T)

  ## subject sample information
  cat("---------------\n", file = file.name.for.txt, append = T)
  cat("Subject sample info:\n", file = file.name.for.txt, append = T)
  cat(paste("Subject sample number:", ncol(subject)),
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)

  cat("Subject sample name:\n", file = file.name.for.txt, append = T)
  cat(colnames(subject), file = file.name.for.txt, append = T)
  cat("\n", file = file.name.for.txt, append = T)

  if (hasQC != "no") {
    ## QC sample information
    cat("\n", file = file.name.for.txt, append = T)
    cat("---------------\n", file = file.name.for.txt, append = T)
    cat("QC sample info:\n", file = file.name.for.txt, append = T)
    cat(paste("QC sample number:", ncol(qc)),
        file = file.name.for.txt,
        append = T)
    cat("\n", file = file.name.for.txt, append = T)

    cat("QC sample name:\n", file = file.name.for.txt, append = T)
    cat(colnames(qc), file = file.name.for.txt, append = T)
    cat("\n", file = file.name.for.txt, append = T)
  }
  else {
    cat("QC sample info:\n", file = file.name.for.txt, append = T)
    cat("No QC in data", file = file.name.for.txt, append = T)
    cat("\n", file = file.name.for.txt, append = T)
  }

  ## other processing information
  cat("\n", file = file.name.for.txt, append = T)
  cat("---------------\n", file = file.name.for.txt, append = T)
  cat("Some processing information\n",
      file = file.name.for.txt,
      append = T)
  cat("MV imputation:",
      MetFlowData[["mv.imputation"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat("Imputation method:",
      MetFlowData[["imputation.method"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat("Zero filter:",
      MetFlowData[["zero.filter"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat("Zero filter criteria:",
      MetFlowData[["zero.filter.criteria"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat("QC outlier filter:",
      MetFlowData[["qc.outlier.filter"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat("Normalization:",
      MetFlowData[["normalization"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat("Normalization method:",
      MetFlowData[["normalization.method"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat("Data integration:",
      MetFlowData[["data.integration"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
  cat("Data integration method:",
      MetFlowData[["data.integration.method"]],
      file = file.name.for.txt,
      append = T)
  cat("\n", file = file.name.for.txt, append = T)
}