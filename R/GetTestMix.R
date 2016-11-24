#' Get internal standatds from data.
#'
#' @title GetTestMix
#' @description Test mix information.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param data data name for text mixture.
#' @param test.mix.is Internal standards information in test information.
#' @param test.mix.data Test mixture data from XCMS.
#' @param rterror rt tolerance for IS.
#' @param mzerror mz tolerance for IS.
#' @export

GetTestMix <- function(data = "Test mixture.csv",
                       test.mix.is = "test.mix.is.csv",
                       test.mix.info = "test.mix.info.csv",
                       mzerror = 15,
                       rterror = 30,
                       rt.unit.is.second = T,
                       plot.output = TRUE,
                       path = NULL) {

   # browser()
  if (is.null(path)) {path <- getwd()}

  ## get data information
  data <- read.csv(data, stringsAsFactors = F)
  is <- read.csv(test.mix.is, stringsAsFactors = F)
  test.mix.info <- read.csv(test.mix.info, stringsAsFactors = F)

  tags <- data[,-grep("Test_mix", colnames(data))]
  sample <- data[,grep("Test_mix", colnames(data))]

  ## order sample and test.mex.info
  test.mix.info <- test.mix.info[order(as.numeric(test.mix.info[,2])),]

  ## sort sample in data according to sample order
  sample.index <- grep("Test_mix", colnames(data))
  sample <- data[,sample.index]
  tags <- data[,-sample.index]
  sample.name <- colnames(sample)
  sample.order <-
    as.numeric(substr(x = sample.name, start = 9, stop = nchar(sample.name)))
  sample <- sample[,order(sample.order)]
  data <- cbind(tags, sample)



  sample.index <- grep("Test_mix", colnames(data))
  sample <- data[,sample.index]
  tags <- data[,-sample.index]

  mz <- as.numeric(tags[, "mzmed"])
  rt <- as.numeric(tags[, "rtmed"])

  sample.name <- test.mix.info[,1]
  sample.order <- as.numeric(test.mix.info[,2])
  sample.batch <- as.numeric(test.mix.info[,3])


  ### name from text.mix.info vs name from data
  if(!identical(sort(colnames(sample)), sort(sample.name))) stop("Name form data and form Test mixture are nor equal!!!")

  ## get IS information
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
  newIS <- matrix(NA, ncol = ncol(is) + ncol(sample))
  colnames(newIS) <- c(colnames(is), colnames(sample))

  for (i in 1:length(is.name)) {
    mz.different <- (is.mz[i] - mz)*10^6/is.mz[i]
    rt.different <- is.rt[i] - rt
    index.for.is <- which(abs(mz.different) <= mzerror & abs(rt.different) <= rterror)
    ## not find
    if (length(index.for.is) == 0)
    {cat(paste(is.name[i], "is not found in data!!!"));cat("\n")
      newIS <- rbind(newIS, c(is[i,], rep(NA, ncol(sample))))
    }

    ## find only one
    if (length(index.for.is) == 1) {
      intensity <- as.numeric(sample[index.for.is,])
      newIS <- rbind(newIS, c(is[i,], intensity))
    }
    ## find more than one
    if (length(index.for.is) > 1) {
      intensity <- sample[index.for.is,]

      temp <- is[rep(i,length(index.for.is)),]
      temp[,1] <- paste(temp[,1],c(1:nrow(temp)), sep="_")
      newIS <- rbind(cbind(temp, intensity))

    }
  }

  newIS <- newIS[-1,]
  newIS <- apply(newIS, 2, as.character)

  if (all(is.na(newIS[,-c(1:3)]))) {
    cat("\n")
    return("Don't find any IS in you data!!!")
  }

  ## add some ststistical results
  newIS.tags <- newIS[,c(1:3)]
  newIS.sample <- newIS[,-c(1:3)]

  sample.rsd <- round(apply(newIS.sample, 1, function(x) {if (all(is.na(as.numeric(x))))
  {return(NA)} else {sd(as.numeric(x), na.rm = T)*100/mean(as.numeric(x), na.rm = T)}}),2)


  rsd.info <- cbind(sample.rsd)
  # browser()
  find.or.not <- apply(newIS.sample, 1, function(x) {if(all(is.na(as.numeric(x)))) {return("No")} else {return("YES")}})

  add.info <- cbind(rsd.info, find.or.not)
  newtest.mix <- cbind(newIS.tags, add.info, newIS.sample)
  write.csv(newtest.mix, file.path(path,"newtest.mix.csv"), row.names = F)
  save(newtest.mix, file = file.path(path,"newtest.mix"))

   ##output some result
  if (plot.output) {
    for (i in 1:nrow(newIS)) {
      if (find.or.not[i] == "No") {next}

      intensity.in.sample <- as.numeric(newIS.sample[i,])

      # IS in samples
      pdf(file.path(path, paste("IS ",i," in Test mixture.pdf",sep = "")), height = 6,width = 8)

      plot(sample.order, intensity.in.sample, xlab="Injection order", ylab = "Intensity",
           main = sprintf("%s\nRSD %s%s", as.character(newIS.tags[i,1]), sample.rsd[i], "%"),
           cex.lab=1.3,cex.axis=1.3, cex.main = 1, pch = 19)
      abline(h = mean(intensity.in.sample), lty = 2, col = "firebrick1")

      dev.off()


    }
  }


}