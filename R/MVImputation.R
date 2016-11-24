#' @title MVimputation
#' @description Impute MV in data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData: MetFlowData.
#' @param imputation.method: Which imputation method you want to use? It
#' contains "knn", "rf" (missForest), "mean", "median", "zero", "minium",
#' "bpca" (BPCA), "svd" (SVD) and "ppca" (PPCA). Default is "knn". The detial of
#' this method can be find in detail and reference paperes.
#' @param knn You can see them from help of impute package: impute.knn.
#' @param missForest You can see them from help of missForest package: missForest.
#' @param BPCA_PPCA_SVD You can see them from help of pcaMethods: bpca, ppca, svd.
#' @return Return a MetFlowData whose MVs have been imputated.
#' @export
#' @seealso The MV imputation methods can see in
#' \code{\link[impute]{impute.knn}}, \code{\link[missForest]{missForest}},
#' \code{\link[pcaMethods]{bpca}}, \code{\link[pcaMethods]{ppca}} and
#' \code{\link[pcaMethods]{svf}}.
#' @references The MV imputation in metabolomics data can see in
#' \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:c9d05d0f-e945-43d0-bb4a-50ea0f90338e}{Guida's} paper.

MVimputation <- function(MetFlowData = MetFlowData,
                         ##MV imputation method
                         imputation.method = "knn",
                         # knn parameters
                         k = 10,
                         rowmax = 0.5,
                         colmax = 0.8,
                         maxp = 1500,
                         rng.seed = 362436069,
                         # missForest parameters
                         maxiter = 10,
                         ntree = 100,
                         decreasing = FALSE,
                         replace = TRUE,
                         classwt = NULL,
                         cutoff = NULL,
                         strata = NULL,
                         sampsize = NULL,
                         nodesize = NULL,
                         maxnodes = NULL,
                         xtrue = NA,
                         parallelize = 'no',
                         #BPCA PPCA, and SVD parameters
                         nPcs = 2,
                         maxSteps = 100,
                         threshold = 1e-04) {
  options(warn = -1)
  # browser()
  #### MV imputation
  if ((sum(is.na(MetFlowData[["subject"]])) + sum(is.na(MetFlowData[["qc"]]))) == 0)
  {
    warning("MVs have been imputed!!!")
    return(MetFlowData)
  }
  qc <- MetFlowData[["qc"]]
  subject <- MetFlowData[["subject"]]
  qc.info <- MetFlowData[["qc.info"]]
  subject.info <- MetFlowData[["subject.info"]]
  tags <- MetFlowData[["tags"]]

  subject.name <- subject.info[, 1]
  qc.name <- qc.info[, 1]

  subject.batch <- subject.info[, 4]
  qc.batch <- qc.info[, 4]

  data <- SplitBatch(MetFlowData =  MetFlowData)
  subject1 <- data[[1]]
  qc1 <- data[[2]]

  var.index <- list()
  for (i in 1:length(subject1)) {
    temp <- cbind(qc1[[i]], subject1[[i]])
    temp <- SXTMVimputation(
      data = temp,
      method = imputation.method,
      # knn parameters
      k = k,
      rowmax = rowmax,
      colmax = colmax,
      maxp = maxp,
      rng.seed = rng.seed,
      # missForest parameters
      maxiter = maxiter,
      ntree = ntree,
      decreasing = decreasing,
      replace = replace,
      classwt = classwt,
      cutoff = cutoff,
      strata = strata,
      sampsize = sampsize,
      nodesize = nodesize,
      maxnodes = maxnodes,
      xtrue = xtrue,
      parallelize = parallelize,
      #BPCA PPCA, and SVD parameters
      nPcs = nPcs,
      maxSteps = maxSteps,
      threshold = threshold
    )
    qc1[[i]] <- temp[,(1:ncol(qc1[[i]]))]
    subject1[[i]] <- temp[,-c(1:ncol(qc1[[i]]))]
  }

  subject2 <- subject1[[1]]
  qc2 <- qc1[[1]]

  if (length(subject1) > 1) {
  for (i in 2:length(subject1)) {
    subject2 <- cbind(subject2, subject1[[i]])
    qc2 <- cbind(qc2, qc1[[i]])
  }
  }

  MetFlowData[["subject"]] <- subject2
  MetFlowData[["qc"]] <- qc2
  MetFlowData[["mv.imputation"]] <- "yes"
  MetFlowData[["imputation.method"]] <- imputation.method
  options(warn = 0)
  return(MetFlowData)
}
