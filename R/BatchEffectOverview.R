#' Using PCA score plot to view the batch effect.
#'
#' @title BatchEffectOverview
#' @description Using PCA score plot to view the batch effect.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData.before MetFlowData before normalization or integration.
#' @param MetFlowData.before MetFlowData after normalization or integration.
#' @param path work directory
#' @return Give PCA score plot for QC and subject.
#' @export

### Batch effect for multiple batch datasets
BatchEffectOverview <- function(MetFlowData.before = MetFlowData1,
                                MetFlowData.after = MetFlowData2,
                                path = NULL) {
  options(warn = -1)
# browser()
  if (is.null(path)) {
    path <- getwd()
  } else {
    dir.create(path)
  }

  hasQC <- MetFlowData.before[["hasQC"]]
  subject.info <- MetFlowData.before[["subject.info"]]
  qc.info <- MetFlowData.before[["qc.info"]]

  subject.bef <- MetFlowData.before[["subject"]]
  qc.bef <- MetFlowData.before[["qc"]]

  subject.aft <- MetFlowData.after[["subject"]]
  qc.aft <- MetFlowData.after[["qc"]]

  if ((sum(is.na(subject.bef)) + sum(is.na(qc.bef))) != 0) {
    stop("The data has MV, please do MV imputation first!!!")
  }

  data.bef <- SplitBatch(MetFlowData = MetFlowData.before)
  subject.bef1 <- data.bef[[1]]
  qc.bef1 <- data.bef[[2]]
  subject.info.bef1 <- data.bef[[3]]
  qc.info.bef1 <- data.bef[[4]]

  data.aft <- SplitBatch(MetFlowData = MetFlowData.after)
  subject.aft1 <- data.aft[[1]]
  qc.aft1 <- data.aft[[2]]
  subject.info.aft1 <- data.aft[[3]]
  qc.info.aft1 <- data.aft[[4]]

  if(hasQC != "no"){
  ##PCA analysis  for QC samples
  qc.name <- qc.info[,1]
  qc.batch <- as.numeric(qc.info[,4])

  qc.info <- list()
  for (i in 1:length(qc.bef1)) {
    qc.info[[i]] <- colnames(qc.bef1[[i]])
  }


  names(qc.info) <- paste("Batch", c(1:length(qc.info)),sep = "")
  ##before
  qc.SXTpcaData <- SXTpca(subject = qc.bef,
                          info = qc.info,
                          QC = F,
                          scale.method = "auto")

  SXTpcaPlot(SXTpcaData = qc.SXTpcaData,
             score.plot.name = "Before Batch effect in QC PCA",
             ellipse = T,
             path = path)
  ##after
  qc.SXTpcaData <- SXTpca(subject = qc.aft,
                          info = qc.info,
                          QC = F,
                          scale.method = "auto")

  SXTpcaPlot(SXTpcaData = qc.SXTpcaData,
             score.plot.name = "After Batch effect in QC PCA",
             ellipse = T,
             path = path)
}
  ##PCA analysis for subject samples
  subject.info <- list()
  for (i in 1:length(subject.bef1)) {
    subject.info[[i]] <- colnames(subject.bef1[[i]])
  }
  names(subject.info) <- paste("Batch", c(1:length(subject.info)),sep = "")

  ##before
  subject.SXTpcaData <- SXTpca(subject = subject.bef,
                               info = subject.info,
                               QC = F,
                               scale.method = "auto")

  SXTpcaPlot(SXTpcaData = subject.SXTpcaData,
             score.plot.name = "Before Batch effect in Subject PCA",
             ellipse = T,
             path = path)

  ##after
  subject.SXTpcaData <- SXTpca(subject = subject.aft,
                               info = subject.info,
                               QC = F,
                               scale.method = "auto")

  SXTpcaPlot(SXTpcaData = subject.SXTpcaData,
             score.plot.name = "After Batch effect in Subject PCA",
             ellipse = T,
             path = path)

  if (hasQC != "no"){
  ## QC sum intensity distribution
  colours.list <- c("palegreen","royalblue","firebrick1","tan1","deepskyblue",
                            "cyan","gray48", "chocolate4","darkmagenta","indianred1")
  colours <- NULL
  for (i in 1:length(qc.bef1)) {
  colours[match(colnames(qc.bef1[[i]]), colnames(qc.bef))] <- colours.list[i]
  }

  pdf(file.path(path, "QC total intensity distribution.pdf"), width = 14, height = 7)
  layout(matrix(c(1:2), ncol = 2))
  par(mar = c(5,5,4,2))
  ##before
  plot(colSums(qc.bef), col = colours, pch = 19, xlab = "QC injection order", ylab = "Total intensity",
       cex.lab = 1.3, cex.axis = 1.3, cex = 1.3,
       ylim = c(0.5*min(colSums(qc.bef)), 1.5*max(colSums(qc.bef))),
       main = "Before")

  qc.scale <- t(apply(qc.bef, 1, function(x) {(x-mean(x))/sd(x)}))
  boxplot(qc.scale, col = colours, xlab = "QC index", ylab = "Intensity (auto scaled)",
          cex.lab = 1.3, cex.axis = 1.3, notch = F, outline = F, main = "Before")

  ##after
  plot(colSums(qc.aft), col = colours, pch = 19, xlab = "QC injection order", ylab = "Total intensity",
       cex.lab = 1.3, cex.axis = 1.3, cex = 1.3,
       ylim = c(0.5*min(colSums(qc.bef)), 1.5*max(colSums(qc.bef))),
       main = "After")

  qc.scale <- t(apply(qc.aft, 1, function(x) {(x-mean(x))/sd(x)}))
  boxplot(qc.scale, col = colours, xlab = "QC index", ylab = "Intensity (auto scaled)",
          cex.lab = 1.3, cex.axis = 1.3, notch = F, outline = F, main = "After")

  dev.off()
  }
  options(warn = -1)
}