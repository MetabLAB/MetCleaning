#' Get internal standatds from data.
#'
#' @title GetIS
#' @description Get internal standatds from data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData for internal standatd getting.
#' @param path Directory for outputing results.
#' @param mzerror mz tolerance for IS.
#' @param rterror rt tolerance for IS.
#' @param rt.unit.is.second RT's unit is second or not? Default is T.
#' @param IS IS data name for reading.
#' @param plot.ouput Output plot or not? Default is TRUE.
#' @return newIS: A new IS data for IS.
#' @export

GetIS <- function(MetFlowData = MetFlowData,
                  mzerror = 15,
                  rterror = 30,
                  rt.unit.is.second = T,
                  IS = "IS.csv",
                  plot.output = TRUE,
                  path = NULL
                  ) {
   browser()
  if (is.null(path)) {path <- getwd()}
  ## get data information
  tags <- MetFlowData[["tags"]]
  subject <- MetFlowData[["subject"]]
  qc <- MetFlowData[["qc"]]
  subject.order <- as.numeric(MetFlowData[["subject.order"]])
  qc.order <- as.numeric(MetFlowData[["qc.order"]])

  mz <- as.numeric(tags[, "mzmed"])
  rt <- as.numeric(tags[, "rtmed"])

  ## get IS information
  is <- read.csv(IS, stringsAsFactors = F)
  is.name <- is[, 1]
  is.mz <- as.numeric(is[, 2])
  is.rt <- as.numeric(is[, 3])

  ## change RT unit
  if (!rt.unit.is.second) {
    rt <- rt * 60
    is.rt <- is.rt * 60
  }

  ## begin find IS in data
  # New IS information file
  newIS <- matrix(NA, ncol=ncol(is)+ncol(qc)+ncol(subject))
  colnames(newIS) <- c(colnames(is), colnames(qc), colnames(subject))

  for (i in 1:length(is.name)) {
    mz.different <- (is.mz[i] - mz)*10^6/is.mz[i]
    rt.different <- is.rt[i] - rt
    index.for.is <- which(abs(mz.different) <= mzerror & abs(rt.different) <= rterror)
    ## not find
    if (length(index.for.is) == 0)
    {cat(paste(is.name[i], "is not found in data!!!"));cat("\n")
      newIS <- rbind(newIS, c(is[i,], rep(NA, ncol(qc) + ncol(subject))))
      }

    ## find only one
    if (length(index.for.is) == 1) {
      intensity.in.subject <- as.numeric(subject[index.for.is,])
      intensity.in.qc <- as.numeric(qc[index.for.is,])
      intensity <- c(intensity.in.qc, intensity.in.subject)

      newIS <- rbind(newIS, c(is[i,], intensity))
    }
    ## find more than one
    if (length(index.for.is) > 1) {
      intensity.in.subject <- subject[index.for.is,]
      intensity.in.qc <- qc[index.for.is,]
      intensity <- cbind(intensity.in.qc, intensity.in.subject)

      temp <- is[rep(i,length(index.for.is)),]
      temp[,1] <- paste(temp[,1],c(1:nrow(temp)), sep="_")
      newIS <- rbind(cbind(temp, intensity))

    }
  }

  newIS <- newIS[-1,]
  newIS <- apply(newIS, 2, as.character)

  if (all(is.na(newIS[,-c(1:3)]))) {
    cat("\n")
    cat("Don't find any IS in you data!!!\n")
    }

  ## add some ststistical results
  newIS.tags <- newIS[,c(1:3)]
  newIS.qc <- newIS[,-c(1:3)][,c(1:ncol(qc))]
  newIS.subject <- newIS[,-c(1:(3 + ncol(qc)))]

  qc.rsd1 <- round(apply(newIS.qc, 1, function(x) {if (all(is.na(as.numeric(x))))
    {return(NA)} else {sd(as.numeric(x), na.rm = T)*100/mean(as.numeric(x), na.rm = T)}}),2)

  qc.rsd2 <- round(apply(newIS.qc, 1, function(x) {if (all(is.na(as.numeric(x))))
  {return(NA)} else {x <- as.numeric(x); x[is.na(x)] <- 0; return((sd(x))*100/mean(x))}}),2)

  subject.rsd1 <- round(apply(newIS.subject, 1, function(x) {if (all(is.na(as.numeric(x))))
  {return(NA)} else {sd(as.numeric(x), na.rm = T)*100/mean(as.numeric(x), na.rm = T)}}),2)

  subject.rsd2 <- round(apply(newIS.subject, 1, function(x) {if (all(is.na(as.numeric(x))))
  {return(NA)} else {x <- as.numeric(x); x[is.na(x)] <- 0; return((sd(x))*100/mean(x))}}),2)

  rsd.info <- cbind(qc.rsd1, qc.rsd2, subject.rsd1, subject.rsd2)
# browser()
  find.or.not <- apply(newIS.qc, 1, function(x) {if(all(is.na(as.numeric(x)))) {return("No")} else {return("YES")}})

  add.info <- cbind(rsd.info, find.or.not)
  newIS <- cbind(newIS.tags, add.info, newIS.qc, newIS.subject)
  write.csv(newIS, file.path(path,"newIS.csv"), row.names = F)
  save(newIS, file = file.path(path,"newIS"))
  ##output some result
  if (plot.output) {
  for (i in 1:nrow(newIS)) {
    if (find.or.not[i] == "No") {next}

    intensity.in.qc <- as.numeric(newIS.qc[i,])
    intensity.in.subject <- as.numeric(newIS.subject[i,])
    ## MV or not?
    mv.in.subject.index <- which(is.na(intensity.in.subject))
    mv.in.qc.index <- which(is.na(intensity.in.qc))
    ## change MV to zero
    intensity.in.qc[mv.in.qc.index] <- 0
    intensity.in.subject[mv.in.subject.index] <- 0
browser()
    # IS in QC samples
    pdf(file.path(path, paste("IS ",i," in QC sample.pdf",sep = "")), height = 6,width = 8)
    colour.qc <- rep(NA, length(qc.order))
    colour.qc[mv.in.qc.index] <- "firebrick1"
    colour.qc[is.na(colour.qc)] <- "black"

    plot(qc.order, intensity.in.qc, xlab="QC order", ylab = "Intensity",
         main = sprintf("%s\nRSD no MV: %s%s\n RSD MV to zero: %s%s",
                        as.character(newIS.tags[i,1]), qc.rsd1[i], "%",qc.rsd2[i], "%"),
         cex.lab=1.3,cex.axis=1.3, cex.main = 1, pch = 19, col = colour.qc)
    abline(h = mean(intensity.in.qc), lty = 2, col = "firebrick1")

    legend("bottomright", legend = "MV in QC",
           col = "firebrick1", pch = 19, pt.cex = 1.3,cex = 1)
    dev.off()

    # IS in subject samples
    pdf(file.path(path, paste("IS ",i," in subject sample.pdf",sep="")), width = 8,height = 6)

    colour.subject <- rep(NA, length(subject.order))
    colour.subject[mv.in.subject.index] <- "firebrick1"
    colour.subject[is.na(colour.subject)] <- "black"

    plot(subject.order, intensity.in.subject, xlab="Subject order", ylab = "Intensity",
         main = sprintf("%s\nRSD no MV: %s%s\n RSD MV to zero: %s%s",
                        as.character(newIS.tags[i,1]), subject.rsd1[i], "%",subject.rsd2[i], "%"),
         cex.lab=1.3,cex.axis=1.3, cex.main = 1, pch = 19, col = colour.subject)
    abline(h = mean(intensity.in.subject), lty = 2, col = "firebrick1")

    legend("bottomright", legend = "MV in Subject",
           col = "firebrick1", pch = 19, pt.cex = 1.3,cex = 1)
    dev.off()
  }
  }
  }
