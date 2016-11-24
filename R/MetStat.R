#' @title MetStat
#' @description A whole work flow for high throughput MS based metabolomics
#' data statistical analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param parameters All the parameters can be found in others functions.
#' @return All the results can be got form other functions and instruction.
#' @export
#' @examples
#' \dontrun{
#' ## load the demo data
#'data(MetFlowData, package = "MetProcesser")
#'data(new.group, package = "MetProcesser")
#'##create a folder for MetStat demo
#'dir.create("Demo for MetStat")
#'setwd("Demo for MetStat")
#'## export the demo data as csv
#'write.csv(new.group, "new.group.csv", row.names = FALSE)
#'## run MetStat
#'MetStat(MetFlowData = MetFlowData, new.group = TRUE)
#' }

MetStat <- function(MetFlowData = MetFlowData,
                    new.group = TRUE,
                    rsd.cutoff = 30,
                    #transformation para
                    log.scale = FALSE,
                    #PCA analysis para
                    QC = TRUE,
                    scale.method = "auto",
                    #PLS analysis para
                    plsmethod = "plsr",
                    #FoldChange para
                    fc = TRUE,
                    to = c("case", "control"),
                    ratio = "median",
                    #UnivariateTest para
                    test.method = "t",
                    adjust.method = "fdr",
                    #MarkerSelection para
                    foldchange = "foldchange",
                    p = "p",
                    vip = "vip",
                    foldchange.cutoff = c(4 / 3, 3 / 4),
                    p.cutoff = 0.05,
                    vip.cutoff = 0,
                    #VolcanoPlot para
                    x = "foldchange",
                    y = "p",
                    z = "vip",
                    col = c("black", "firebrick1"),
                    #heatmap para
                    variable = "all",
                    Group = c("control", "case"),
                    #MarkerShow para
                    beeswarm = TRUE,
                    pls.analysis = TRUE,
                    pca.analysis = TRUE,
                    uni.test = TRUE,
                    path = NULL) {
  options(warn = -1)
  # browser()
  if (is.null(path)) {
    path = getwd()
  } else{
    path <- dir.create(path)
  }

  path.inter <- file.path(path, "intermediate")
  dir.create(path.inter)

  ##RSD filtering
  MetFlowData <-
    RSDfilter(MetFlowData = MetFlowData, rsd.cutoff = rsd.cutoff)
  #save data
  met.data.rsd.filter <- MetFlowData
  save(met.data.rsd.filter,
       file = file.path(path.inter, "met.data.rsd.filter"))

  ##ReChangeGroup
  if (new.group) {
    cat("---------------------------------------------------------------------\n")
    cat("Change group information\n")
    if (all(dir() != "new.group.csv"))
      stop("No new.group information!!!")
    met.data <- ReChangeGroup(MetFlowData = MetFlowData)
    #save data
    met.data.new.group <- met.data
    save(met.data.new.group,
         file = file.path(path.inter, "met.data.new.group"))
  }

  if (pca.analysis) {
    cat("---------------------------------------------------------------------\n")
    cat("PCA analysis\n")
    ##PCA analysis
    PCAanalysis(
      MetFlowData = met.data,
      QC = QC,
      log.scale = log.scale,
      scale.method = scale.method,
      path = file.path(path, "PCA analysis")
    )
  }

  if (pls.analysis) {
    cat("---------------------------------------------------------------------\n")
    cat("PLS analysis\n")
    ##PLS analysis
    PLSanalysis(
      MetFlowData = met.data,
      log.scale = log.scale,
      scalemethod = scale.method,
      path = file.path(path, "PLS analysis"),
      plsmethod = plsmethod
    )

    load(file.path(path, "PLS analysis", "vip"))
    vip <- apply(vip, 2, mean)
  }

  ##VIP
  tags <- met.data[["tags"]]
  tags <- data.frame(tags, vip)
  met.data[["tags"]] <- tags

  #save data
  met.data.vip <- met.data
  save(met.data.vip, file = file.path(path.inter, "met.data.vip"))

  cat("---------------------------------------------------------------------\n")
  cat("Heat map...\n")
  HeatMap(
    MetFlowData = met.data,
    log.scale = log.scale,
    variable = variable,
    Group = Group,
    scale.method = scale.method,
    path = file.path(path, "heat map")
  )

  if (fc) {
    ##FoldChange
    cat("---------------------------------------------------------------------\n")
    met.data <- FoldChange(MetFlowData = met.data,
                           to = to,
                           ratio = ratio)
    #save data
    met.data.fc <- met.data
    save(met.data.fc, file = file.path(path.inter, "met.data.fc"))
  }

  if (uni.test) {
    ##UnivariateTest
    cat("---------------------------------------------------------------------\n")
    cat("Univariate test\n")
    met.data <- UnivariateTest(
      MetFlowData = met.data,
      test.method = test.method,
      adjust.method = adjust.method,
      log.scale = log.scale
    )
    #save data
    met.data.uni.test <- met.data
    save(met.data.uni.test,
         file = file.path(path.inter, "met.data.uni.test"))
  }

  ##MarkerSelection
  cat("---------------------------------------------------------------------\n")
  cat("Marker selection\n")
  met.data <- MarkerSelection(
    MetFlowData = met.data,
    foldchange = foldchange,
    p = p,
    vip = vip,
    foldchange.cutoff = foldchange.cutoff,
    p.cutoff = p.cutoff,
    vip.cutoff = vip.cutoff,
    path = file.path(path, "marker selection")
  )
  #save data
  met.data.after.stat <- met.data
  save(met.data.after.stat,
       file = file.path(path.inter, "met.data.after.stat"))

  ##Export data
  ExportData(
    MetFlowData = met.data,
    data.name = "data.after.stat",
    subject.info.name = "subject.info",
    qc.info.name = "qc.info"
  )

  #save data
  met.data.after.stat <- met.data
  save(met.data.after.stat,
       file = file.path(path.inter, "met.data.after.stat"))

  ##VolcanoPlot
  cat("---------------------------------------------------------------------\n")
  cat("Draw volcano plot\n")
  VolcanoPlot(
    MetFlowData = met.data,
    x = x,
    y = y,
    z = z,
    col = c("black", "firebrick1"),
    foldchange.cutoff = foldchange.cutoff,
    p.cutoff = p.cutoff,
    vip.cutoff = vip.cutoff,
    path = file.path(path, "marker selection")
  )

  ##MarkerShow
  cat("---------------------------------------------------------------------\n")
  cat("Marker show\n")
  MarkerShow(
    MetFlowData = met.data,
    beeswarm = beeswarm,
    path = file.path(path, "marker selection")
  )
  options(warn = 0)

  cat("---------------------------------------------------------------------\n")
  cat("MetStat is done!!!\n")
}