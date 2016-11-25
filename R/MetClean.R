#' @title MetClean
#' @description A whole work flow for high throughput MS based metabolomics data
#'cleaning.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param parameters All the parameters can be found in others functions.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @details The manual of MetCleaning can be found in \href{https://github.com/jaspershen/MetCleaning/blob/master/vignettes/MetCleaning.pdf}{github}
#' or in \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:d8b5fcff-c725-4689-a97d-7ff106322fb6}{my library}.
#' @examples
#' \dontrun{
#' #load the demo data
#' data(data, package = "MetCleaning")
#' data(sample.information, package = "MetCleaning")
#' ##create a folder for MetClean demo
#' dir.create("Demo for MetClean")
#' setwd("Demo for MetClean")
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#' #run MetClean
#' MetClean(#ImportData para
#' data = "data.csv",
#' sample.information = "sample.information.csv",
#' polarity = "positive",
#' #DataNormalization
#' method = "svr",
#' threads = 2)
#' }

MetClean <- function(#ImportData para
                   data = "data.csv",
                   sample.information = "sample.information.csv",
                   polarity = "positive",
                   hasIS = "no",
                   hasQC = "yes",
                   #MVFilter para
                   obs.mv.cutoff = 0.5,
                   var.mv.cutoff = 0.5,
                   #MVimputation
                   imputation.method = "knn",
                   k = 10,
                   rowmax = 0.5,
                   colmax = 0.5,
                   maxp = 1500,
                   #ZeroFilter para
                   obs.zero.cutoff = 0.5,
                   var.zero.cutoff = 0.5,
                   #DataNormalization
                   method = "svr",
                   threads = 2,
                   #PeakIdentification
                   hmdb.matching = TRUE,
                   show = 5,
                   mass.tolerance = 30,
                   mz.tolerance = 30,
                   rt.tolerance = 180,
                   #DataOverview para
                   met.plot = TRUE,
                   path = NULL,
                   worklist.from = "manual"){

  # browser()
  if (is.null(path)) {
    path <- getwd()
  } else{
    dir.create(path)
  }

  path.inter <- file.path(path, "intermediate")
  dir.create(path.inter)

  options(warn = -1)

  #----------------------------------------------------------------------------
  #read data
  cat("Importing data...\n")
  if (all(dir("intermediate") != "met.data.raw")) {
    met.data <- ImportData(data = data,
                           sample.information = sample.information,
                           polarity = polarity,
                           hasIS = hasIS,
                           hasQC = hasQC,
                           worklist.from = worklist.from)
    #save data
    met.data.raw <- met.data
    save(met.data.raw, file = file.path(path.inter, "met.data.raw"))
  } else {
    load(file.path(path.inter, "met.data.raw"))
    met.data <- met.data.raw
    }

  # browser()
  batch <- unique(met.data[["subject.info"]][,4])
  subject <- met.data[["subject"]]

  #----------------------------------------------------------------------------
  if (sum(is.na(subject)) != 0){
  #观察缺失值分布
  MZoverview(MetFlowData = met.data,
             path = file.path(path,"1MV overview"),
             what = "mv")

  cat("---------------------------------------------------------------------\n")
  cat("Missing values filter...\n")
  #对缺失值进行筛选
  met.data <- MZfilter(MetFlowData = met.data,
                       obs.per.cutoff = obs.mv.cutoff,
                       var.per.cutoff = var.mv.cutoff,
                       what = "mv",
                       path = file.path(path,"2MV filter"))
  #save data
  met.data.mv.filter <- met.data
  save(met.data.mv.filter, file = file.path(path.inter, "met.data.mv.filter"))

  cat("---------------------------------------------------------------------\n")
  cat("Missing values imputation...\n")
  #对缺失值进行补齐
  met.data <- MVimputation(MetFlowData = met.data,
                           ##MV imputation method
                           imputation.method = imputation.method,
                           # knn parameters
                           k = k,
                           rowmax = rowmax,
                           colmax = colmax,
                           maxp = maxp)
  #save data
  met.data.mv.imputation <- met.data
  save(met.data.mv.imputation, file = file.path(path.inter, "met.data.mv.imputation"))
  }

  subject <- met.data[["subject"]]

  #观察零值分布
  MZoverview(MetFlowData = met.data,
             what = "zero",
             path = file.path(path,"3Zero overview"))

  cat("---------------------------------------------------------------------\n")
  cat("Zero filter...\n")
  #对零进行筛选
  met.data <- MZfilter(MetFlowData = met.data,
                       obs.per.cutoff = obs.zero.cutoff,
                       var.per.cutoff = var.zero.cutoff,
                       what = "zero",
                       path = file.path(path,"4Zero filter"))
  #save data
  met.data.zero.filter <- met.data
  save(met.data.zero.filter, file = file.path(path.inter, "met.data.zero.filter"))

  cat("---------------------------------------------------------------------\n")
  cat("Peak identification...\n")
  #物质鉴定 peak identification
  if (any(dir() == "peak identification")) {
    if (all(dir("intermediate") != "met.data.peak.iden")) {
      met.data <- PeakIdentification(MetFlowData = met.data,
                                     ##parameters for matching
                                     mz.tolerance = mz.tolerance,
                                     rt.tolerance = rt.tolerance)
      #save data
      met.data.peak.iden <- met.data
      save(met.data.peak.iden, file = file.path(path.inter, "met.data.peak.iden"))
    }else{
      load(file.path(path.inter, "met.data.peak.iden"))
      met.data <- met.data.peak.iden
    }
  }

  if (hmdb.matching) {
    #物质鉴定 mass identification

    if (all(dir("intermediate") != "met.data.mass.iden")) {
      met.data <- MassIdentification(MetFlowData = met.data,
                                     mass.tolerance = mass.tolerance,
                                     polarity = "positive",
                                     show = 5)
      #save data
      met.data.mass.iden <- met.data
      save(met.data.mass.iden, file = file.path(path.inter, "met.data.mass.iden"))
    }else{
      load(file.path(path.inter, "met.data.mass.iden"))
      met.data <- met.data.mass.iden
    }
  }

  cat("---------------------------------------------------------------------\n")
  cat("QC outlier filtering...\n")

   #标出QC的outlier
  if (all(dir("intermediate") != "met.data.qc.outlier.filter")) {
  met.data <- QCOutlierFilter(MetFlowData = met.data,
                              CI = 0.95,
                              path = "5QC outlier filter")
  met.data.qc.outlier.filter <- met.data
  save(met.data.qc.outlier.filter,
       file = file.path(path.inter,"met.data.qc.outlier.filter"))
  } else {
    load(file.path(path.inter, "met.data.qc.outlier.filter"))
    met.data <- met.data.qc.outlier.filter
  }

  cat("---------------------------------------------------------------------\n")
  cat("Data normalization...\n")
  ##Data Normalization
  if (length(batch) > 1) {
    peak.plot = FALSE
  }else{
      peak.plot = TRUE
    }
  if (all(dir("intermediate") != "met.data.nor")) {
  met.data <- DataNormalization(MetFlowData = met.data,
                                method = method,
                                threads = threads,
                                path = path,
                                peakplot = peak.plot)
  #save data
  met.data.nor <- met.data
  save(met.data.nor, file = file.path(path.inter, "met.data.nor"))
  }else{
    load(file.path(path.inter, "met.data.nor"))
    met.data <- met.data.nor
  }

  cat("---------------------------------------------------------------------\n")
  cat("Subject outlier filtering...\n")
  #标出subject的outlier
  if (all(dir("intermediate") != "met.data.subject.outlier.filter")) {
  met.data <- SubjectOutlierFilter(MetFlowData = met.data,
                                   CI = 0.95,
                                   path = file.path("6Subject outlier finder"))
  met.data.subject.outlier.filter <- met.data
  save(met.data.subject.outlier.filter,
       file = file.path(path.inter, "met.data.subject.outlier.filter"))
  } else {
    load(file.path(path.inter, "met.data.subject.outlier.filter"))
    met.data <- met.data.subject.outlier.filter
  }


  cat("---------------------------------------------------------------------\n")
  cat("Data integration...\n")

  if(length(batch) > 1) {
  ## data integration
  met.data <- DataIntegration(MetFlowData = met.data)
  #save data
  met.data.integration <- met.data
  save(met.data.integration, file = file.path(path.inter, "met.data.integration"))
  }

  if (length(batch) > 1) {
    #batch effect
    BatchEffectOverview(MetFlowData.before = met.data.zero.filter,
                        MetFlowData.after = met.data,
                        path = file.path(path, "8Batch effect"))

    if(met.plot) {
      cat("---------------------------------------------------------------------\n")
      cat("Metabolite plot...\n")
    #metabolite plot
    MetabolitePlot(MetFlowData.before = met.data.zero.filter,
                   MetFlowData.after = met.data,
                   path = file.path(path, "9metabolite plot"))
    }
  }

  cat("---------------------------------------------------------------------\n")
  cat("RSD overview...\n")
  #RSD分布
  RSDoverview(MetFlowData.before = met.data.zero.filter,
              MetFlowData.after = met.data,
              path = file.path(path,"10RSD overview"))

  #data overview
  DataOverview(MetFlowData = met.data,
               feature.distribution = TRUE,
               path = file.path("11Data overview"))

  #ouput data
  ExportData(MetFlowData = met.data,
             data.name = "data_after_pre",
             subject.info.name = "subject.info",
             qc.info.name = "qc.info")

  #save data
  met.data.after.pre <- met.data
  save(met.data.after.pre, file = file.path(path, "met.data.after.pre"))
  cat("---------------------------------------------------------------------\n")
  cat("MetClean is done!!!\n")
  options(warn = 0)
}



##RSDoverview function
RSDoverview <- function(MetFlowData.before = MetFlowData1,
                        MetFlowData.after = MetFlowData2,
                        path = NULL) {
  if (is.null(path)) {
    path <- getwd()
  }else{
    dir.create(path)
  }

  hasQC <- MetFlowData.before[["hasQC"]]
  if (hasQC == 'no') {stop("Data has no QC!!!")}
  ##before
  subject.bef <- MetFlowData.before[["subject"]]
  qc.bef <- MetFlowData.before[["qc"]]
  subject.info.bef <- MetFlowData.before[["subject.info"]]
  qc.info.bef <- MetFlowData.before[["qc.info"]]

  ##split data
  data.bef <- SplitBatch(MetFlowData = MetFlowData.before)
  subject.bef1 <- data.bef[[1]]
  qc.bef1 <- data.bef[[2]]

  ##after
  subject.aft <- MetFlowData.after[["subject"]]
  qc.aft <- MetFlowData.after[["qc"]]
  subject.info.aft <- MetFlowData.after[["subject.info"]]
  qc.info.aft <- MetFlowData.after[["qc.info"]]

  ##split data
  data.aft <- SplitBatch(MetFlowData = MetFlowData.after)
  subject.aft1 <- data.aft[[1]]
  qc.aft1 <- data.aft[[2]]


  pdf(file.path(path, "RSD distribution in different batch.pdf"),width = 14, height = 7)
  layout(matrix(c(1:2), ncol = 2))
  par(mar = c(5,5,4,2))
  for (i in 1:length(qc.bef1)) {
    #before
    QC <- qc.bef1[[i]]
    QC.rsd <- apply(QC, 1, function(x) {sd(x)*100/mean(x)})
    par(mar = c(5,5,4,2))
    colours1 <- rep(NA, length(QC.rsd))
    colours1[QC.rsd>30] <- "firebrick1"
    colours1[is.na(colours1)] <- "black"

    plot(QC.rsd, xlab = "Feature index", ylab = "Relative Standard Deviation (RSD, %)",
         cex.lab = 1.3, cex.axis = 1.3, pch = 19, col = colours1, main = paste("Batch",i,"QC (Before)"),
         cex.main = 1.3)
    abline(h = 30, lty = 2)
    legend("topleft",paste("RSD<30%: ",round(sum(QC.rsd < 30)/length(QC.rsd),4)*100,"%"), bty = "n",
           cex = 1.3)

    #after
    QC <- qc.aft1[[i]]
    QC.rsd <- apply(QC, 1, function(x) {sd(x)*100/mean(x)})
    par(mar = c(5,5,4,2))
    colours1 <- rep(NA, length(QC.rsd))
    colours1[QC.rsd>30] <- "firebrick1"
    colours1[is.na(colours1)] <- "black"

    plot(QC.rsd, xlab = "Feature index", ylab = "Relative Standard Deviation (RSD, %)",
         cex.lab = 1.3, cex.axis = 1.3, pch = 19, col = colours1, main = paste("Batch",i,"QC (After)"),
         cex.main = 1.3)
    abline(h = 30, lty = 2)
    legend("topleft",paste("RSD<30%: ",round(sum(QC.rsd < 30)/length(QC.rsd),4)*100,"%"), bty = "n",
           cex = 1.3)

  }
  dev.off()

  pdf(file.path(path, "RSD distribution in all batches.pdf"), width = 14, height = 7)
  layout(matrix(c(1:2),ncol = 2))
  par(mar=c(5,5,4,2))
  # before
  qc.rsd <- apply(qc.bef, 1, function(x) {sd(x)*100/mean(x)})
  colours1 <- rep(NA, length(qc.rsd))
  colours1[qc.rsd>30] <- "firebrick1"
  colours1[is.na(colours1)] <- "black"
  plot(qc.rsd, xlab = "Feature index", ylab = "Relative Standard Deviation (RSD, %)",
       cex.lab = 1.3, cex.axis = 1.3, pch = 19, col = colours1,cex.main = 1.3,
       main = "Before")
  abline(h = 30, lty = 2)
  legend("topleft",paste("RSD<30%: ",round(sum(qc.rsd < 30)/length(qc.rsd),4)*100,"%"), bty = "n",
         cex = 1.3)

  # after
  qc.rsd <- apply(qc.aft, 1, function(x) {sd(x)*100/mean(x)})
  colours1 <- rep(NA, length(qc.rsd))
  colours1[qc.rsd>30] <- "firebrick1"
  colours1[is.na(colours1)] <- "black"
  plot(qc.rsd, xlab = "Feature index", ylab = "Relative Standard Deviation (RSD, %)",
       cex.lab = 1.3, cex.axis = 1.3, pch = 19, col = colours1,cex.main = 1.3,
       main = "After")
  abline(h = 30, lty = 2)
  legend("topleft",paste("RSD<30%: ",round(sum(qc.rsd < 30)/length(qc.rsd),4)*100,"%"), bty = "n",
         cex = 1.3)

  dev.off()
  layout(1)
}