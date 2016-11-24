#' @title MRMdataMerge
#' @description Merge different MRM batch datasets into one dataset.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param ... The MRM data to be merged.
#' @param path Directory for reading data and outputing results.
#' @param merge.term.index Use which term to merge datasets?
#' @param tags.index Which columna are tags in your data?
#' @param tags.same.in.datasets Which columna are same in you datasets?And those information only remain one in merging.
#' @param tags.should.be.mean.in.datasets Which columna should be averaged in the mergin
#' @export

MRMdataMerge <- function(...,
                         path =NULL,
                         merge.method = "name",
                         #merge according which?
                         merge.term.index = c(1),
                         #tags.index
                         tags.index = c(1:17),
                         #the information which are same in different datastets
                         tags.same.in.datasets = c(1:4),
                         #the information which should be mean in different datastets
                         tags.should.be.mean.in.datasets = c(11:13),
                         #the information which should be remain in different datastets
                         tags.should.be.remain.in.datasets = NULL,
                         #parameters for mz RT merge
                         mz.tolerance = 15,
                         rt.tolerance =15
                         ) {
# browser()
if (is.null(path)) path <- getwd()

  data <- list(...)

tags <- list()
sample <- list()
for (i in 1:length(data)) {
  tags[[i]] <- data[[i]][,tags.index]
  sample[[i]] <- data[[i]][,-tags.index]
}

union.name <- tags[[1]][,1]

for (i in 2:length(tags)) {
  union.name <- c(union.name,tags[[i]][,merge.term.index])
}

union.name <- unique(union.name)

new.sample <- list()
new.tags <- list()

for (i in 1:length(tags)) {
  this.name <- tags[[i]][,merge.term.index]
  this.tags <- tags[[i]]
  index <- match(this.name, union.name)

  new.sample[[i]] <- matrix(NA, ncol = ncol(sample[[i]]), nrow = length(union.name))
  colnames(new.sample[[i]]) <- colnames(sample[[i]])
  new.sample[[i]][index,] <- as.matrix(sample[[i]])

  new.tags[[i]] <- matrix(NA, ncol = ncol(tags[[i]]), nrow = length(union.name))
  colnames(new.tags[[i]]) <- colnames(tags[[i]])
  new.tags[[i]][index,] <- as.matrix(tags[[i]])

}


##combine new.sample
sample <- new.sample[[1]]
for (i in 2:length(new.sample)) {
  sample <- cbind(sample,new.sample[[i]])
}

##combine new.tags
tags1 <- new.tags[[1]][,tags.same.in.datasets]
tags2 <- new.tags[[1]][,tags.should.be.mean.in.datasets]
tags3 <- new.tags[[1]][,tags.should.be.remain.in.datasets]

for (i in 2:length(new.tags)) {
  new.tags1 <- new.tags[[i]][,tags.same.in.datasets]
  new.tags2 <- new.tags[[i]][,tags.should.be.mean.in.datasets]
  new.tags3 <- new.tags[[i]][,tags.should.be.remain.in.datasets]

  ##combine tags1
  for (j in 1:ncol(tags1)) {
    temp <- cbind(tags1[,j], new.tags1[,j])
    temp <- apply(temp, 1, function(x) { if (all(is.na(x))) {NA}
      else if (all(!is.na(x))) {x[1]} else {x[!is.na(x)]}})
    tags1[,j] <- temp
  }
# browser()
  ##combine tags2
  for (k in 1:ncol(tags2)) {
    temp <- cbind(as.numeric(tags2[,k]), as.numeric(new.tags2[,k]))
    temp <- apply(temp, 1, function(x) {mean(x,na.rm = T)})
    temp[is.nan(temp)] <- NA
    tags2[,k] <- temp
  }

  ##combine tags3
  tags3 <- cbind(tags3, new.tags3)
}

tags <- cbind(tags1, tags2, tags3)
data <- cbind(tags, sample)
write.csv(data,"new.data.csv",row.names = F)
return(data)
}