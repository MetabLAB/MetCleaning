#' @title VolcanoPlot
#' @description Draw volcano plot.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param path The directory you want to write results.
#' @param foldchange_p_foldchange.cutoff_p.cutoff See MarkerSeletion.
#' @param col Colour for markers and non-markers.
#' @return volcano plot.
#' @export

VolcanoPlot <- function(MetFlowData = MetFlowData,
                        x = "foldchange",
                        y = "p",
                        z = "vip",
                        col = c("black", "firebrick1"),
                        foldchange.cutoff = c(4/3, 3/4),
                        p.cutoff = 0.05,
                        vip.cutoff = 1,
                        path = NULL){
  # browser()
  if(is.null(path)) {path <- gewtd()}
  else{dir.create(path)}

  tags <- MetFlowData[["tags"]]

  # if (foldchange %in% colnames(tags)) {
  #   cat(paste("Use", foldchange, "as fold change.\n"))
  #   cat("\n")
  # }else {
  #   stop(paste("No", foldchange,"!!!"))
  # }
  #
  # if (p %in% colnames(tags)) {
  #   cat(paste("Use", p, "as p value.\n"))
  #   cat("\n")
  # }else {
  #   stop(paste("No", p,"!!!"))
  # }
  #
  # if (vip %in% colnames(tags)) {
  #   cat(paste("Use", vip, "as VIP value.\n"))
  #   cat("\n")
  # }else {
  #   stop(paste("No", vip,"!!!"))
  # }

  x1 <- as.numeric(tags[,x])
  y1 <- as.numeric(tags[,y])
  if (z != FALSE) {
    z1 <- as.numeric(tags[,z])
  }

  f.cutoff1 <- as.numeric(foldchange.cutoff[1])
  f.cutoff2 <- as.numeric(foldchange.cutoff[2])
  p.cutoff <- as.numeric(p.cutoff)
  vip.cutoff <- as.numeric(vip.cutoff)

  ##x transformation
  if (x == "foldchange"){
    x1 <- log(x1, 2)
    x.lab1 <- "log2(Fold change)"
  }
  if (x == "p") {
    x1 <- -log(x1, 10)
    x.lab1 <- "-log10(p value)"
  }
  if (x == "vip") {
    x1 <- x1
    x.lab1 <- "VIP"
  }

  ##y transformation
  if (y == "foldchange"){
    y1 <- log(y1, 2)
    y.lab1 <- "log2(Fold change)"
  }
  if (y == "p") {
    y1 <- -log(y1, 10)
    y.lab1 <- "-log10(p value)"
  }
  if (y == "vip") {
    y1 <- y1
    y.lab1 <- "VIP"
  }

  ##z transformation
  if (z != FALSE) {
    if (z == "foldchange"){
      z1 <- log(z1, 2)
    }
    if (z == "p") {
      z1 <- -log(z1, 10)
    }
    if (z == "vip") {
      z1 <- z1
    }
  }

  if ("is.marker" %in% colnames(tags)) {
    marker.index <- which(tags[,"is.marker"] == "yes")
    if (length(marker.index) == 0) {stop("No marker are selected, please change canditios and try again!!!")}
  } else {
    stop("Please select marker first (Using MarkerSelection function)!!!")
  }

  colour <- rep(NA, length(p))
  colour[marker.index] <- col[2]
  colour[is.na(colour)] <- col[1]

  if (z != FALSE) {
    temp1 <- c(min(z1), max(z1))
    temp2 <- c(0.5, 1.3)
    lm.reg <- lm(temp2 ~ temp1)
    cexa <- lm.reg[["coefficients"]][2] * z1 +  lm.reg[["coefficients"]][1]
  } else {
    cexa <- 1
  }
# browser()
  pdf(file.path(path,"Volcano plot.pdf"))
  par(mar = c(5,5,4,2))
  plot(x1, y1, xlab = x.lab1, ylab = y.lab1,
       cex.lab = 1.3, cex.axis = 1.3, col = colour, pch  = 19,
       cex = cexa)
  # abline(h = -log(p.cutoff,10), lty = 2)
  # abline(v = log(f.cutoff1,2), lty = 2)
  # abline(v = log(f.cutoff2,2), lty = 2)
  vip.mean <- (min(z1)+max(z1))/2
  cexa.mean <- round(lm.reg[["coefficients"]][2]*vip.mean+lm.reg[["coefficients"]][1],1)
  if (z != FALSE) {
    legend("topleft", legend = c(round(min(z1),1), round(vip.mean,1), round(max(z1))),
           title = "VIP", bty = "n", cex = 1.3,
           pt.cex = c(0.5,cexa.mean,1.3), pch = 19)
  }

  dev.off()
}