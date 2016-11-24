VIP <- function(MetFlowData = MetFlowData,
                   #used data
                   log.scale = 10,
                   scalemethod="auto",
                   plsmethod = "plsreg")
  #parameter setting
{
  # browser()
  options(warn = -1)

  subject <- MetFlowData[["subject"]]
  qc <- MetFlowData[["qc"]]
  subject.info <- MetFlowData[["subject.info"]]
  group <- subject.info[,"group"]
  group.unique <- sort(unique(group))
  subject.name <- subject.info[,1]

  if (is.null(qc)) {QC <- FALSE}

  info <- list()
  for (i in 1:length(group.unique)) {
    info[[i]] <- subject.name[which(group == group.unique[i])]
  }

  names(info) <- group.unique

  #load needed packages
  need.packages1 <- c("pls","plsdepot")

  packages <- library()[[2]][,1]
  for (i in need.packages1) {
    if (!any(packages == i)) {install.packages(i)}
  }

  require(plsdepot)
  require(pls)

  int <- t(subject)
  index <- NULL
  for (i in 1:length(info)) {
    index1 <- as.character(info[[i]])
    index <- c(index,index1)
  }
  if (length(which(index==""))!=0)  {index<-index[-which(index=="")]}

  index <- index[!is.na(index)]
  index <- match(index, rownames(int))
  index <- index[!is.na(index)]
  int <- int[index, ]

  #######
  name <- rownames(int)
  # browser()
  Y <- NULL
  label <- list()
  for (i in 1:length( info )) {
    label[[i]] <- match(as.character(info[[i]]),name)
    label[[i]] <- label[[i]][!is.na(label[[i]])]
    Y[label[[i]]] <- i-1
  }
  # browser()
  int <- log(int+1, log.scale)
  int.scale <- SXTscale(int,method = scalemethod)
  # int.Y<-SXTscale(Y,method=scalemethod)
  int.Y <- Y
  ncompa <- nrow(int) - 1

  if (plsmethod == "plsr") {
    pls1 <- plsr(int.Y~int.scale,scale = FALSE,validation = "CV",ncomp = ncompa,method = "oscorespls")
    save(pls1,file = "pls1")

    #########select the number of compents#################
    msep <- MSEP(pls1)
    save(msep,file = file.path(path, "msep"))
    msep <- msep$val[,,]

    yesnot <- "y"
    while (yesnot == "y"|yesnot == "") {
      comps.number <-readline("How many comps do you want to see? ")
      while (!exists("comps.number")|comps.number=="")
      {cat("You must give a comps number to continute!!!\n")
        comps.number <- readline("How many comps do you want to see? ")}
      comps.number <- as.numeric(comps.number)
      plot(x = c(1:comps.number),y=msep[1,2:(comps.number+1)],type="b",col="firebrick1",pch=20,
           xlab = "ncomp",ylab = "MSEP",cex.lab=1.3,cex.axis=1.3)
      points(x=c(1:comps.number),y=msep[2,2:(comps.number+1)],type="b",pch=2)
      legend("top",legend = c("CV","adjCV"),col = c("firebrick1","black"),pch=c(20,2),
             bty = "n",cex = 1.3,pt.cex = 1.3)
      yesnot <- readline("Do you want to see the next plot? (y/n)")
    }

    pdf(file.path(path, "MSEP plot.pdf"))
    plot(x=c(1:comps.number),y=msep[1,2:(comps.number+1)],type="b",col="firebrick1",pch=20,
         xlab="ncomp",ylab="MSEP",cex.lab=1.3,cex.axis=1.3)
    points(x=c(1:comps.number),y=msep[2,2:(comps.number+1)],type="b",pch=2)
    legend("top",legend = c("CV","adjCV"),col = c("firebrick1","black"),pch=c(20,2),
           bty = "n",cex = 1.3,pt.cex = 1.3)
    dev.off()

    number<-readline("Please type number and press Enter  to continute:  ")
    while (!exists("number")|number=="") {cat("You must give a number to continute!!!\n")
      number<-readline("Please type comps number value and press Enter  to continute: ")}
    number<-as.numeric(number)

    ##################construct final pls model###################
    pls2 <- plsr(int.Y~int.scale, scale = FALSE,validation = "CV",ncomp = number,method = "oscorespls")
    save(pls2, file = file.path(path, "pls2"))
    vip <- SXTvip(pls2)
    vip <- apply(vip, 1, mean)
  }

  else {
    # browser()
    require(SXTdummy)
    dummy <- SXTdummy(Y)
    # int.dummy<-SXTscale(dummy,method=scalemethod)
    int.dummy <- dummy
    # ncompa = nrow(int.scale) - 1
    ncompa <- min(nrow(int), ncol(int))
    pls1 <- plsreg1(int.scale,Y, comps = ncompa)
    save(pls1,file = file.path(path, "pls1"))
    #########select the number of compents#################
    Q2cum <- pls1$Q2[,5]
    Q2cum[is.nan(Q2cum)] <- 1
    yesnot <- "y"
    while (yesnot=="y"|yesnot=="") {
      comps.number<-readline("How many comps do you want to see? ")
      while (!exists("comps.number")|comps.number=="") {cat("You must give a comps number to continute!!!\n")
        comps.number<-readline("How many comps do you want to see? ")}
      comps.number<-as.numeric(comps.number)
      barplot(Q2cum[1:comps.number],xlab="ncomp",ylab="Q2cum",cex.lab=1.3,cex.axis=1.3)
      a <- barplot(Q2cum[1:comps.number],xlab="ncomp",ylab="Q2cum",cex.lab=1.3,cex.axis=1.3)
      abline(h=0)
      points(a,Q2cum[1:comps.number],type="b",col="red",pch=20,cex=2)
      yesnot <- readline("Do you want to see the next plot? (y/n)")
    }

    number <- readline("Please type number and press Enter  to continute:  ")
    while (!exists("number")|number=="") {cat("You must give a number to continute!!!\n")
      number <- readline("Please type comps number value and press Enter  to continute: ")}
    number <- as.numeric(number)

    ##################construct final pls model###################
    cat(paste("Construct PLS model with all peaks using",number,"comps ...","\n"))
    pls2 <- plsreg1(int.scale,Y,comps = number)
    pls.temp <- plsreg2(int.scale,int.dummy, comps = number)
    vip <- pls.temp$VIP
    vip <- apply(vip, 1, mean)
  }

  tags <- MetFlowData[["tags"]]
  tags <- data.frame(tags, vip)
  MetFlowData[["tags"]] <- tags
  return(MetFlowData)
}





