#' @title ImportData
#' @description Import data and retrun a standard MetProcessor data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param data Data name for analysis. Default is "data.csv". Data are csv
#' format from XCMS, MZmine or other software.
#'  Please see the demo data in example.
#' @param sample.information Sample information name for analysis. Default is
#' "sample.information.csv". Column 1 is sample.name, column 2 is
#' injection.order, column 3 is class ("Subject" or "QC"), column 4 is
#' batch information, column 5 is group ("control" or "case"), other columns
#' are information for sample. Please see demo data in example.
#' @param polarity The polarity of data, "positive", "negative"" or "none"",
#' default is positive.
#' @param hasQC The data has QC samples or not? Default is "yes".
#' @param peak.identification The data has identification result or not?
#' Default is "no".
#' @return  Return a standard MetProcesser dataset.
#' @seealso \code{\link{ExportData}}
#' @export

ImportData <- function(data = "data.csv",
                       sample.information = "sample.information.csv",
                       polarity = "positive",
                       posfix = NULL,
                       qc.has.order = FALSE,
                       worklist.from = "manual",
                       hasIS = "no",
                       hasQC = "yes",
                       peak.identification = "no",
                       path = NULL) {
  # browser()
  if (is.null(path)) path <- getwd()
  if (worklist.from != "GetWorklist")
  {data <- ChangeSampleName(data = data,
                            sample.information = sample.information,
                            polarity = polarity,
                            posfix = posfix,
                            qc.has.order = qc.has.order,
                            output = FALSE,
                            path = path)
  sample.information <- read.csv("sample.information1.csv", stringsAsFactors = F)}

  else {data <- read.csv(file.path(path,data), stringsAsFactors = F)
  sample.information <- read.csv(file.path(path,sample.information), stringsAsFactors = F)}

  ## read sample information
  sample.information <- sample.information[!is.na(sample.information[,1]),]

  ## sort sample information according to sample order
  sample.information <- sample.information[order(as.numeric(sample.information[,2])),]
  write.csv(sample.information,file.path(path, "sample.information1.csv"), row.names = F)

  ## sort sample in data according to sample order
  sample.index <- grep("Sample", colnames(data))
  sample <- data[,sample.index]
  tags <- data[,-sample.index]
  sample.name <- colnames(sample)
  sample.order <-
    as.numeric(substr(x = sample.name, start = 7, stop = unlist(lapply(gregexpr("_", sample.name),function(x) {x[1][1]}))-1))
  sample <- sample[,order(sample.order)]
  data <- cbind(tags, sample)

  ## get subject, qc and tags
  data.name <- colnames(data)
  sample.index <- grep("Sample", data.name)
  sample <- data[,sample.index]
  tags <- data[,-sample.index]
  sample.name <- colnames(sample)

  qc.index <- grep("QC", sample.name)
  ## has QC or not?
  if(length(qc.index) == 0) {hasQC = "no"}

  if (hasQC == "no") {qc <- NULL; subject <- sample; qc.name <- NULL}
  else {qc <- sample[,qc.index]; subject <- sample[,-qc.index]; qc.name <- colnames(qc)}

  subject.name <- colnames(subject)

  ## the sample name from sample information vs the sample name from data
  a <- setdiff(sample.information[,1], sample.name)
  b <- setdiff(sample.name, sample.information[,1])
  if (length(a) != 0)
  {stop(paste(paste(a,collapse = " "),"are in sample information but not in data!!!"))}

  if (length(b) != 0)
  {stop(paste(paste(b,collapse = " "),"are in data but not in sample information!!!"))}

  ## the sample order form sample information vs the sample order form sample name
  sample.name.from.sample.information <- sample.information[,1]
  sample.order.from.sample.information <-
    substr(x = sample.name.from.sample.information,start = 7,
           stop = unlist(lapply(gregexpr("_", sample.name.from.sample.information),function(x) {x[1][1]}))-1)

  sample.different.index <- which(sample.order.from.sample.information != sample.information[,2])
  if (length(sample.different.index) != 0)
  { print(sample.information[sample.different.index,c(1,2)])
  stop("The sample order from data and form sample information are different!!!")}

  if (hasQC == "no") {qc.order <- NULL}
  else {
  qc.order <- substr(x = qc.name,start = 7,stop = unlist(lapply(gregexpr("_", qc.name),function(x) {x[1][1]}))-1)
  }

  subject.order <- substr(x = subject.name,start = 7,stop = unlist(lapply(gregexpr("_", subject.name),function(x) {x[1][1]}))-1)

  subject.name <- substr(subject.name, start = unlist(lapply(gregexpr("_", subject.name),function(x) {x[1][1]}))+1,
                          stop = unlist(lapply(gregexpr("_", subject.name),function(x) {x[2][1]}))-1)

  if (hasQC == "no") {qc.name <- NULL}
  else {
  qc.name <- substr(qc.name, start = unlist(lapply(gregexpr("_", qc.name),function(x) {x[1][1]}))+1,
                         stop = unlist(lapply(gregexpr("_", qc.name),function(x) {x[2][1]}))-1)
  }
  colnames(subject) <- subject.name
  names(subject.order) <- subject.name
  if (hasQC != "no") {
  colnames(qc) <- qc.name
  names(qc.order) <- qc.name
  }
  ## remove sample order form sample.name.from.data and sample.name.from.sample.information
  sample.name.from.data <- colnames(sample)
  sample.name.from.sample.information <- sample.information[,1]
    sample.name.from.data <-
    substr(sample.name.from.data, start = unlist(lapply(gregexpr("_", sample.name.from.data),function(x) {x[1][1]}))+1,
    stop = unlist(lapply(gregexpr("_", sample.name.from.data),function(x) {x[2][1]}))-1)

    sample.name.from.sample.information <-
      substr(sample.name.from.sample.information, start = unlist(lapply(gregexpr("_", sample.name.from.sample.information),function(x) {x[1][1]}))+1,
             stop = unlist(lapply(gregexpr("_", sample.name.from.sample.information),function(x) {x[2][1]}))-1)

    colnames(sample) <- sample.name.from.data
    sample.information[,1] <- sample.name.from.sample.information

  subject.info <- sample.information[sample.information[,3]=="Subject",]
  if (hasQC == "no") {qc.info <- NULL}
  else {
  qc.info <- sample.information[sample.information[,3]=="QC",]
  }

  if ((sum(is.na(subject)) + sum(is.na(qc))) == 0) {
    mv.imputation <- "yes"
    imputation.method <- "default"
  }
  else {
    mv.imputation <- "no"
    imputation.method <- "no"
  }

  name <- tags[,"name"]
  if(polarity == "positive")
  {tags[,"name"] <- paste(name,"POS", sep = "_")
  pol <- rep("POS", nrow(tags))
  tags <- data.frame(pol,tags)
  colnames(tags)[1] <- "polarity"}
  if(polarity == "negative")
  {tags[,"name"] <- paste(name,"NEG", sep = "_")
  pol <- rep("NEG", nrow(tags))
  tags <- data.frame(pol,tags)
  colnames(tags)[1] <- "polarity"}

  MetFlowData <- list(subject = subject,
                      qc = qc,
                      tags = tags,
                      subject.info = subject.info,
                      qc.info = qc.info,
                      subject.order = subject.order,
                      qc.order = qc.order,
                      ## some preprocessing information
                      mv.imputation = mv.imputation,
                      imputation.method = imputation.method,
                      zero.filter = "no",
                      zero.filter.criteria = "no",
                      qc.outlier.filter = "no",
                      normalization = "no",
                      normalization.method = "no",
                      data.integration = "no",
                      data.integration.method = "no",
                      hasIS = hasIS,
                      hasQC = hasQC,
                      peak.identification = peak.identification)
  class(MetFlowData) <- "MetFlowData"
  return(MetFlowData)
}

