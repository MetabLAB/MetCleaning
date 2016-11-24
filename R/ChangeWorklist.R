#' @title ChangeWorklist
#' @description Change date in worklist from GetWorklist.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param date Date you want to change in you worklist. Default is today.
#' @return Return a new worklist.
#' @export


ChangeWorklist <-
  function(date = gsub("-", "", as.character(Sys.Date()))) {
    #used to change the date in worklist
    options(warn = -1)
    #change the date
    file <- dir()
    file <- file[!file.info(file)$isdir]

    packages <- library()[[2]][, 1]
    filestyle <- substr(file, regexpr("\\.", file)[[1]] + 1, nchar(file))

    if (filestyle == "xlsx" & all(packages != "xlsx"))
    {
      stop("You R has no xlsx packages, you must install xlsx or translate your file to csv")
    }

    if (filestyle == "xlsx" & any(packages == "xlsx"))
    {
      require(xlsx)
      worklist <-
        read.xlsx(file, 1)
    }#batch design column 1 is Sample.Name

    if (filestyle == "csv")
    {
      worklist <- read.csv(file)
    }#batch design column 1 is Sample.Name

    if (filestyle != "xlsx" & filestyle != "csv")
    {
      stop("The format of file is wrong, it must be xlsx or csv")
    }


    Data.File <- worklist[, grep("Data.File", colnames(worklist))]
    date.old <- substr(file, 1, 8)
    Data.File <- gsub(date.old, date, Data.File)
    worklist[, grep("Data.File", colnames(worklist))] <- Data.File
    worklistname <- paste(date, substr(file, 9, (nchar(file) - 5)), sep =
                            "")

    write.csv(worklist, sprintf('%s.csv', worklistname), row.names = FALSE)
    cat("The date has been changed.\n")
  }