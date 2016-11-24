#' @title MarkerShow
#' @description Draw boxplot for markers.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param path The directory you want to write results.
#' @param beeswarm Do you want draw beeswarm on the boxplot? Deafult is TRUE.
#' @return Box plot of markers.
#' @export

MarkerShow <- function(MetFlowData = MetFlowData,
                       beeswarm = TRUE,
                       path = NULL) {
  # browser()
  options(warn = -1)
  if (is.null(path)) {
    path <- getwd()
  } else{
    dir.create(path)
  }

  library(beeswarm)
  tags <- MetFlowData[["tags"]]
  subject <- MetFlowData[["subject"]]
  subject.info <- MetFlowData[["subject.info"]]
  if (all(colnames(tags) != "is.marker")) {
    stop("Please select marker first(use MarkerSelection function)!!!")
  }

  is.marker <- tags[, "is.marker"]
  marker.index <- which(is.marker == "yes")
  feature.name <- tags[, "name"]
  group <- subject.info[, "group"]
  subject.name <- subject.info[, 1]
  group.unique <- sort(unique(group))

  info <- list()
  for (i in 1:length(group.unique)) {
    info[[i]] <- subject.name[which(group == group.unique[i])]
  }
  names(info) <- group.unique

  for (i in marker.index) {
    temp.data <- list()
    for (j in 1:length(info)) {
      temp.data[[j]] <-
        as.numeric(subject[i,])[match(info[[j]], subject.name)]
    }
    names(temp.data) <- names(info)
    pdf(file.path(path, paste(feature.name[i], "boxplot.pdf")))
    par(mar = c(5, 5, 4, 2))
    boxplot(
      temp.data,
      xlab = "Group",
      ylab = "Intensity",
      cex.lab = 1.5,
      cex.axis = 1.3
    )
    if (beeswarm) {
      beeswarm(temp.data,
               pch = 19,
               add = T,
               col = "grey")
    }
    dev.off()
  }
  options(warn = 0)
}
