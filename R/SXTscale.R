SXTscale<-function(x,method=c("pareto","auto")) {
  if (method=="pareto") {x<-apply(x,2, function(x) {(x-mean(x))/sqrt(sd(x))})}
  if (method=="auto") {x<-apply(x,2, function(x) {(x-mean(x))/sd(x)})}
  return(x)
}