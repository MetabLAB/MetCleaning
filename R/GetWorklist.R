#' Generate Worklist for data acquisition.
#'
#' @title GetWorklist
#' @description Generate Worklist for data acquisition.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param instrument Which instrument you use? "Agilent" or "AB", default is "Agilent".
#' @param name The name of worklist.
#' @param randommethod Which random method you want to use? "no", "position" or "injection". Default is "no".
#' @param samplenumber Sample number.
#' @param replication Replication times.
#' @param QCstep QC step.
#' @param conditionQCnumber Condition QC number.
#' @param testmixstep Test mixture step.
#' @param injectionfrom Injection order from which? Default is 1.
#' @param dir Directory.
#' @return New worklist.
#' @export

GetWorklist <- function(instrument = "Agilent",
                        name = "worklist",
                        randommethod = "no",
                        samplenumber = NULL,
                        replication = 1,
                        QCstep = 8,
                        conditionQCnumber = 10,
                        testmixstep = 32,
                        injectionfrom = 1,
                        dir = "D:\\MassHunter\\Data\\SXT\\"){
  #names is the name of the folder,plates is the used plates,
  #if AB,dir is ""
# browser()
  options(warn = -1)
  file <- dir()
  packages <- library()[[2]][, 1]
  if (instrument == "AB") {
    dir = ""
  }

  filestyle <- substr(file, regexpr("\\.", file)[[1]] + 1, nchar(file))

  if (filestyle=="xlsx"&
      all(packages != "xlsx"))
  {
    stop("Please install R package: xlsx or translate your file to csv!!!\n")
  }

  if (filestyle == "xlsx" & any(packages == "xlsx"))
  {
    library(xlsx)
    x <- read.xlsx(file[file == "batch design.xlsx"], 1)
  }#batch design column 1 is Sample.Name

  if (filestyle == "csv")
  {
    x <- read.csv(file[file == "batch design.csv"])
  }#batch design column 1 is Sample.Name

  if (filestyle != "xlsx" & filestyle != "csv")
  {
    stop("The format of file is wrong, it must be xlsx or csv!!!")
  }

  options(warn = 0)
  x <- as.character(x[, 1])
  na.number <- sum(is.na(x))
  x <- x[!is.na(x)]
  space.number <- sum(x == "")
  x <- x[x != ""]

  cat(paste("\nThere are",na.number,"NAs and",space.number,
            "spaces in your batch design, please confirm there are no error!\n")
      )

  if (is.null(samplenumber)) {
    samplenumber <- length(x)
  }
  else {
    if (samplenumber > length(x)) {
      samplenumber <- samplenumber
      warning("The sample number you set is larger than the sample in your batch design!!!\n")
    }
    else {
      samplenumber <- samplenumber
    }

  }
  x <- x[1:samplenumber]

  options(warn = -1)
  plate1 <- rep(1, 51)
  plate2 <- rep(2, 51)
  plate3 <- rep(3, 51)
  plate4 <- rep(4, 51)
  plate5 <- rep(5, 51)
  plate6 <- rep(6, 51)

  plate <- rep(c(plate1, plate2), 6)
  plate <- plate[1:(samplenumber * replication)]
  real.plate <- rep(c(plate1, plate2, plate3, plate4, plate5, plate6), 6)
  real.plate <- real.plate[1:(samplenumber * replication)]
  vial.position <- rep(c(1:51), 6)
  vial.position <- vial.position[1:(samplenumber * replication)]

sub.position <- list()
  for (i in 1:6) {
    sub.position[[i]] <- paste(LETTERS[i],c(1:9), sep = "")
  }
sub.position <- unlist(sub.position)

real.position <- list()
for (i in 1:12) {
  real.position[[i]] <- paste(paste("P",i,sep=""),sub.position, sep="-")
}

real.position <- unlist(real.position)

position <- rep(real.position[1:108], 6)

  position <- position[1:(samplenumber * replication)]
  real.position <- real.position[1:(samplenumber * replication)]


  sub.position96 <- c(paste("A",c(1:12),sep=""),paste("B",c(1:12),sep=""),paste("C",c(1:12),sep=""),paste("D",c(1:12),sep=""),
                   paste("E",c(1:12),sep=""),paste("F",c(1:12),sep=""),
                   paste("G",c(1:12),sep=""),paste("H",c(1:12),sep=""))
  real.position96 <-
    rep(c(paste("P1",sub.position96,sep="-"), paste("P2",sub.position96, sep="-"),
        paste("P3",sub.position96,sep="-"), paste("P4",sub.position96, sep="-"),
        paste("P5",sub.position96,sep="-"), paste("P6",sub.position96, sep="-")),4)
  real.position96 <- real.position96[1:(samplenumber*replication)]

  position96 <- rep(c(paste("P1",sub.position96,sep="-"), paste("P2",sub.position96, sep="-")),8)
  position96 <- position96[1:(samplenumber*replication)]
  ###repeat samples?
  if (replication==1) { x<-x} else{
    x2 <- NULL
    for ( i in 1 : replication) {
      x1 <- paste( x, i,sep="_")
      x2 <- cbind(x2,x1)
    }
    x <- x2
  }

  x <- as.character(x)

  #random position or random injection order or no
  if (instrument=="Agilent") {
    if (randommethod == "position") {
      random.order <- sample(1:(samplenumber*replication))
      x <- data.frame(random.order, x)
      x <- x[order(as.numeric(x[,1])),]
      x <- as.character(x[,-1])
      x <- cbind(x,position,real.position,position96, real.position96)
      colnames(x) <- c("Sample Name","Position in 54 plate","Real position in 54 plate",
                       "Position in 96 plate", "Real position in 96 plate")
      write.csv(x,sprintf("%s sample info.csv",name))
      x <- x[,-c(3,4,5)]
    }

    if (randommethod == "injection") {
      if (length(x) > 108)
        {warning("The sample number is larger than 108, injection order random is not commended!!!")}

      x <- cbind(x,position,real.position,position96, real.position96)
      colnames(x) <- c("Sample Name","Position in 54 plate","Real position in 54 plate",
                       "Position in 96 plate", "Real position in 96 plate")
      write.csv(x,sprintf("%s sample info.csv",name))
      random.order <- sample(1:(samplenumber*replication))
      x <- data.frame(random.order, x)
      x <- x[order(x[,1]),]
      x <- x[, -c(1, 4, 5, 6)]
    }

    if (randommethod == "no") {
      x <- cbind(x,position,real.position,position96, real.position96)
      colnames(x) <- c("Sample Name","Position in 54 plate","Real position in 54 plate",
                       "Position in 96 plate", "Real position in 96 plate")
      write.csv(x, sprintf("%s sample info.csv",name))
      x <- x[,-c(3,4,5)]
    }

  }

##AB instrument
  if (instrument == "AB") {
    if (randommethod == "position") {
      random.order <- sample(1:(samplenumber*replication))
      x <- data.frame(random.order,x)
      x <- x[order(as.numeric(x[,1])),]
      x <- x[,-1]
      x <- cbind(x, plate, vialposition, real.plate, position96, real.position96)
      colnames(x) <- c("Sample Name","Plate","Position in 54 plate","Real plate of 54 plate",
                       "Position in 96 plate","Real position 96 plate")
      write.csv(x,"sample info.csv")
      x <- x[,-c(4,5,6)]
    }

    if (randommethod=="injection"){
      x <- cbind(x,plate, vialposition, real.plate, position96, real.position96)
      colnames(x) <- c("Sample Name","Plate","Position in 54 plate","Real Plate of 54 plate",
                       "Position in 96 plate", "Real position in 96 plate")
      write.csv(x,"sample info.csv")
      random.order <- sample(1:(samplenumber*replication))
      x <- cbind(random.order, x)
      x <- x[order(as.numeric(x[, 1])), ]
      x <- x[, -c(1, 5, 6, 7)]
    }

    if (randommethod=="no") {
      x <- data.frame(x,plate, vialposition, real.plate, position96, real.position96)
      colnames(x) <- c("Sample Name","Plate","Position in 54 plate","Real Plate of 54 plate",
                       "Position in 96 plate", "Real position in 96 plate")
      write.csv(x,"sample info.csv")
      x <- x[,-c(4,5,6)]
    }

  }

  #now x column 1 is Sample Name, column 2 is Sample Position
  if ( instrument == "Agilent") {
    Blank <- c("Blank","Vial1")
    Test.mix <- c("Test_mix","Vial2")
    QC <- c("QC","Vial3")
    Blank.QC <- rbind(Blank,QC)
  }
  else {
    Blank <- c("Blank", "1", "52")
    Test.mix <- c("Test_mix", "1", "53")
    QC <- c("QC", "1", "54")
    Blank.QC <- rbind(Blank, QC)

  }

  #insert Blank and QC
  x <- lapply(seq(1,nrow(x),by = QCstep), function(y) if(y + QCstep - 1 <= (samplenumber*replication)) {x[y:(y + QCstep-1),]}
  else {x[y:nrow(x),]})

  x <- lapply(x, function(y) rbind(Blank.QC,y))

  x2 <- NULL
  for (i in 1:length(x)) {
    x1 <- x[[i]]
    x2 <- rbind(x2,x1)
  }
  x <- x2
  x <- rbind(x, Blank.QC)
  #insert Test.mix
  if (testmixstep == 0) {x = x}
  else {
  x <- lapply(seq(1,nrow(x), by = testmixstep), function(y) if ( y+testmixstep-1<=nrow(x)) {x[y:(y+testmixstep-1),]} else {x[y:nrow(x),]} )
  x <- lapply(x, function(y) rbind(Test.mix,y))

  x3 <- NULL
  for (i in 1:length(x)) {
    x1 <- x[[i]]
    x3 <- rbind(x3,x1)
  }
  x <- x3
  x <- rbind(x, Test.mix)
  }

  if ( instrument == "Agilent") {
    colnames(x) <- c('Sample.Name', "Sample.Position")
  }
  if (instrument == "AB") {
    colnames(x) <- c('Sample.Name',"Plate.Position","Vial.Postion")
  }

  if (instrument == "Agilent")
  {
    x <- rbind(matrix(rep(Blank,3),ncol=2,byrow = TRUE),matrix(rep(QC,conditionQCnumber),ncol=2,byrow=TRUE),x,matrix(rep(Blank,3),ncol=2,byrow=TRUE))
  }
  if (instrument == "AB") {
    x <- rbind(matrix(rep(Blank,3),ncol = 3,byrow=TRUE),matrix(rep(QC,conditionQCnumber),ncol=3,byrow=TRUE),x,matrix(rep(Blank,3),ncol=3,byrow=TRUE))
  }

  Blank.number <- length(grep("Blank",x[,1]))
  x[,1][grep("Blank",x[,1])] <- paste("Blank",c(1:Blank.number),sep="")

  Test.mix.number <- length(grep("Test_mix",x[,1]))
  x[,1][grep("Test_mix",x[,1])] <- paste("Test_mix",c(1:Test.mix.number),sep="")

  QC.number <- length(grep("QC",x[,1]))
  x[,1][grep("QC",x[,1])][1:conditionQCnumber] <- paste("Condition_QC",c(1:conditionQCnumber),sep="")
  x[,1][grep("QC",x[,1])][conditionQCnumber+1:QC.number] <- paste("QC",c(1:(QC.number-conditionQCnumber)),sep="")
  first <- which(x[,1]=="QC1")
  last <- which(x[,1]==sprintf("QC%s",QC.number-conditionQCnumber))


  before.info <- x[1:(first-1),]
  Data.File1 <- before.info[,1]

  after.info <- x[(last+1):nrow(x),]
  Data.File5 <- after.info[,1]

  middle.info <- x[first:last, ]
  middle.info <- cbind(middle.info, c(1:nrow(middle.info)))
  Sample.QC <-
    middle.info[setdiff(1:nrow(middle.info), grep("Blank", middle.info[, 1])), ]#remove Blank in middle.info
  Sample.QC <-
    Sample.QC[setdiff(1:nrow(Sample.QC), grep("Test_mix", Sample.QC[, 1])), ]#remove Test.mix in middle.info

  middle.blank <-
    middle.info[grep("Blank", middle.info[, 1]), ,drop = FALSE]#column 3 is number
  middle.testmix <- middle.info[grep("Test_mix", middle.info[, 1]), ,drop = FALSE]

  Data.File2 <- paste("Sample",c(injectionfrom:(nrow(Sample.QC)+injectionfrom-1)),sep="")
  Data.File2 <- paste(Data.File2, Sample.QC[,1], sep="_")

  if (instrument == "Agilent") {
    Data.File2 <- cbind(Data.File2, Sample.QC[, 3])
    Data.File3 <- middle.blank[, c(1, 3)]
    Data.File4 <- middle.testmix[, c(1, 3)]
  }
  if (instrument == "AB") {
    Data.File2 <- cbind(Data.File2, Sample.QC[, 4])
    Data.File3 <- middle.blank[, c(1, 4)]
    Data.File4 <- middle.testmix[, c(1, 4)]
  }

  Data.File2 <- rbind(Data.File2,Data.File3,Data.File4)
  Data.File2 <- Data.File2[order(as.numeric(Data.File2[,2])),]
  Data.File2 <- Data.File2[,-2]

  Data.File <- c(Data.File1, Data.File2, Data.File5)
  name.POS <- paste(name,"POS")
  name.NEG <- paste(name,"NEG")
  Data.File.POS <- paste(dir,name.POS,"\\",Data.File,sep="")
  Data.File.NEG <- paste(dir,name.NEG,"\\",Data.File,sep="")
  Data.File.POS <- paste(Data.File.POS,"POS",sep="_")
  Data.File.NEG <- paste(Data.File.NEG,"NEG",sep="_")

  if (instrument == "Agilent")
  {Data.File.POS <- paste(Data.File.POS,"d",sep=".")
    Data.File.NEG<-paste(Data.File.NEG,"d",sep=".")}

  x.POS <- cbind(x,Data.File.POS)
  x.NEG <- cbind(x,Data.File.NEG)
  write.csv(x.POS, sprintf("%s POS.csv",name),row.names=FALSE)
  write.csv(x.NEG,sprintf("%s NEG.csv",name),row.names=FALSE)

  cat("Worklist is generated.\n")
}