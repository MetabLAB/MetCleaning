#' Select markers using fold change and p values.
#'
#' @title MarkerSelection
#' @description Select markers using fold change and p values.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param path The directory you want to write results.
#' @param foldchange which foldchange information your want to use?
#' @param p which p value information your want to use?
#' @param foldchange.cutoff The cutoff value of fold change.
#' @param p.cutoff The cutoff of p value.
#' @return MetFlowData: MetFlowData which has been added is.marker information into tags.
#' @return maker.information.csv: The marker information.
#' @export

MarkerSelection <- function(MetFlowData = MetFlowData,
                            foldchange = "foldchange",
                            p = "p.correct",
                            vip = "vip",
                            foldchange.cutoff = c(4/3, 3/4),
                            p.cutoff = 0.05,
                            vip.cutoff = 1,
                            path = NULL) {
   # browser()
  if (is.null(path)) {path <- getwd()}
  else {dir.create(path)}

  tags <- MetFlowData[["tags"]]
  if (any(colnames(tags) == "is.marker")) {
    warning("Markers have been selected!!!")
  }

  f.cutoff1 <- as.numeric(foldchange.cutoff[1])
  f.cutoff2 <- as.numeric(foldchange.cutoff[2])
  p.cutoff <- as.numeric(p.cutoff)
  vip.cutoff <- as.numeric(vip.cutoff)

  ##foldchange selection
  if (foldchange != FALSE) {
  if (all(colnames(tags) != "foldchange")) {
    stop("Please calculate fold change first!!!")
  }
    foldchange <- as.numeric(tags[,foldchange])
    marker.index1 <- which(foldchange > f.cutoff1 | foldchange < f.cutoff2)
  }
  else {
    marker.index1 <- c(1:nrow(tags))
  }

  ##vip selection
  if (vip != FALSE) {
    if (all(colnames(tags) != "vip")) {
      stop("Please calculate VIP first!!!")
    }
    vip <- as.numeric(tags[,"vip"])
    marker.index2 <- which(vip > vip.cutoff)
  } else {
    marker.index2 <- c(1:nrow(tags))
  }

  ##p selection
  if (p != FALSE) {
  if (all(colnames(tags) != "p")) {
    stop("Please calculate p value first!!!")
  }
    p <- as.numeric(tags[,p])
    marker.index3 <- which(p < p.cutoff)
  } else {
    marker.index3 <- c(1:nrow(tags))
  }

  marker.index <- intersect(intersect(marker.index1, marker.index2), marker.index3)
  if (length(marker.index) == 0) {stop("No marker are selected, please change canditios and try again!!!")}
  cat(paste("There are",length(marker.index), "variables are selected as marker.\n"))
  marker.info <- tags[marker.index,]
  write.csv(marker.info, file.path(path,"marker.info.csv"), row.names = F)
  is.marker <- rep(NA, nrow(tags))
  is.marker[marker.index] <- "yes"
  is.marker[is.na(is.marker)] <- "no"

  if (any(colnames(tags) == "is.marker")) {
    tags[,"is.marker"] <- is.marker
  }else{
  tags <- cbind(tags, is.marker)
  }

  MetFlowData[["tags"]] <- tags
  MetFlowData[["marker.selection.condition"]] <-
  paste("foldchange:",foldchange,"fc:",
        paste(foldchange.cutoff,collapse = ","),"p:",p,"pc:",
        p.cutoff,"vip:",vip,"vc:",vip.cutoff)
  return(MetFlowData)
}