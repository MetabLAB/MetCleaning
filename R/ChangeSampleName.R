#' Change sample name to standard name from GetWorklist.
#'
#' @title ChangeSampleName
#' @description Change sample name to standard name from GetWorklist.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param data data for analysis. Default is "data.csv"
#' @param sample.information sample information name for analysis. First column
#' is sample.name, second column is injection.order, third column is
#' class (Subject or QC), fouth column is batch information,
#' fifth column is group (control or case), latter columns are other
#' information for sample. Please see demo data in example.
#' @param polarity The polarity of data. positive, negative or none, default is
#' positive. If data is positive mode, the name of sample should be added as
#'  "POS" or "NEG" as posfix. Please see the demo data in example.
#' @param posfix The posfix of the data. For example, "POS" or "NEG". Please
#' see the demo data in examples.
#' @return Return a data whose name is standard name for MetFlowData.


ChangeSampleName <- function(data = "data.csv",
           sample.information = "sample.information.csv",
           polarity = "positive",
           posfix = NULL,
           qc.has.order = FALSE,
           output = TRUE,
           path = NULL) {
  # browser()
if (is.null(path)) path <- getwd()
  data <- read.csv(file.path(path,data), stringsAsFactors = F)
  sample.information <- read.csv(file.path(path,sample.information), stringsAsFactors = F)
  sample.information <- sample.information[!is.na(sample.information[,1]),]

  ## sort sample information according to sample order
  sample.information <- sample.information[order(as.numeric(sample.information[,2])),]
  write.csv(sample.information, file.path(path,"sample.information.csv"), row.names = F)

  sample.name <- as.character(sample.information[,1])
  injection.order <- as.numeric(sample.information[,2])
  class <- sample.information[,3]
  batch <- sample.information[,4]

  ## any sample in sample information are not in data?
  index <- match(sample.name, colnames(data))
  sample.not.in.data.index <- which(is.na(index))
  if (length(sample.not.in.data.index) != 0)
  {
    stop(paste(paste(sample.name[sample.not.in.data.index], collapse = " "),
               "in sample information are not found in data!!!"))
  }

  ## sort sample in data according to sample order
  sample <- data[,index]
  tags <- data[,-index]
  data <- cbind(tags, sample)
  index <- match(sample.name, colnames(data))

  if (!is.null(posfix)) {
    if (polarity == "positive") {sample.name <- substr(sample.name,start = 1,stop = unlist(gregexpr(".POS", sample.name))-1)}
    if (polarity == "negative") {sample.name <- substr(sample.name,start = 1,stop = unlist(gregexpr(".NEG", sample.name))-1)}
     }

  subject.index <- grep("Subject",class)
  qc.index <- grep("QC",class)

  subject.name <- sample.name[subject.index]
  qc.name <- sample.name[qc.index]

  subject.order <- injection.order[subject.index]
  qc.order <- injection.order[qc.index]

  subject.name <- paste(paste("Sample", subject.order, sep=""),subject.name,sep="_")
  if (qc.has.order) {qc.name <- paste(paste("Sample", qc.order, sep=""),qc.name,sep="_")}
   else {
     qc.name <- paste("QC", rank(qc.order), sep="")
     qc.name <- paste(paste("Sample", qc.order, sep=""),qc.name,sep="_")
     }

  sample.name[subject.index] <- subject.name
  sample.name[qc.index] <- qc.name
  if (polarity == "positive") {sample.name <- paste(sample.name, "POS", sep="_")}
  if (polarity == "negative") {sample.name <- paste(sample.name, "NEG", sep="_")}
  if (polarity == "none") {sample.name <- paste(sample.name, "NONE", sep= "_")}

  colnames(data)[index] <- sample.name
  sample.information[,1] <- sample.name
  write.csv(sample.information, file.path(path, "sample.information1.csv"), row.names = F)
  if (output) {write.csv(data, file.path(path,"data1.csv"),row.names = F)}
  return(data)
  }