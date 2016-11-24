SXTpaste<-function(x,sep=" ") {
  y<-NULL
  for (i in 1:length(x)) {
    if (i==1) {y<-paste(y,x[i],sep="")}
    else {y<-paste(y,x[i],sep=sep)}
  }
  y
}