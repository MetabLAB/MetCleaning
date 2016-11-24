#' @title SXTMTmatch
#' @description Match two data according to mz and RT.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param data1 First data for matching, first column must be mz and seconod column must be rt.
#' @param datad2 Second data for matching, first column must be mz and seconod column must be rt.
#' @param mz.tolerance mz tolerance for ms1 and ms2 data matching.
#' @param rt.tolerance RT tolerance for ms1 and ms2 data matching.
#' @return Return a result which give the matching result of data1 and database.
#' @export

SXTMTmatch <- function(data1,
                       data2,
                       mz.tolerance = 25,
                       rt.tolerance = 180) {
    if (nrow(data1) == 0 | nrow(data2) == 0) {
      result <- NULL
      return(result)
    }
    mz1 <- as.numeric(data1[, 1])
    rt1 <- as.numeric(data1[, 2])

    mz2 <- as.numeric(data2[, 1])
    rt2 <- as.numeric(data2[, 2])

    result <- NULL
    cat("finished: %")
    cat("\n")
    for (i in 1:length(mz1)) {
      mz.error <- abs(mz1[i] - mz2) * 10 ^ 6 / mz1[i]
      rt.error <- abs(rt1[i] - rt2)
      j <- which(mz.error <= mz.tolerance & rt.error <= rt.tolerance)
      if (length(j) != 0) {
        result1 <-
          cbind(i, j, mz1[i], mz2[j], mz.error[j], rt1[i], rt2[j], rt.error[j])
        result <- rbind(result, result1)
      }

      count <- floor((length(mz1)) * c(seq(0, 1, 0.01)))
      if (any(i == count)) {
        cat(ceiling (i * 100 / length(mz1)))
        cat(" ")
      }

    }
    cat("\n")
    if (is.null(result)) {
      cat("There are not any peak be matched\n,please change the mz or rt tolerance and try again")
      cat("\n")
    }
    else {
      number1 <- length(unique(result[, 1]))
      number2 <- length(unique(result[, 2]))
      cat(
        paste(
          "There are",
          number1,
          "peaks in data1, and",
          number2,
          "peaks in data2 are matched"
        )
      )
      cat("\n")
      colnames(result) <-
        c("Index1",
          "Index2",
          "mz1",
          "mz2",
          "mz error",
          "rt1",
          "rt2",
          "rt error")
      return(result)
    }
  }
