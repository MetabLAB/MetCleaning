#' @title SXTMinifrac
#' @description Remove feature or samples according to the MV or zero ratio.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param ... Datases you want to remove features and samples.
#' @param filter.item The item you want to filter. MV (missing values) or zero. Default is MV.
#' @param filter.rule: For feature filtering, which rule you want to use, intersect or union. Default is intersect.
#' @param minifrac.variable The cutoff for variable. Default is 0.5.
#' @param minifrac.observation: The cutoff for observation. Default is 0.5.
#' @return Return a SXTMinifracData.
#' @examples
#' \dontrun{
#' ## Generate
#' data1 <- matrix(1:20, ncol = 4)
#' data2 <- matrix(21:40, ncol = 4)
#' ## Give MV in to data1 and data2
#' data1[sample(1:20,5)] <- NA
#' data2[sample(1:20,5)] <- NA
#' ## Run SXTMinifrac
#'SXTMinifracData <- SXTMinifrac(data1, data2,
#'                               filter.item = "MV",
#'                               filter.rule = "intersect",
#'                               minifrac.variable = 0.4,
#'                               minifrac.observation = 0.4)
#' attributes(SXTMinifracData)
#' print(SXTMinifracData)
#' }

SXTMinifrac <- function(...,
                        filter.item = "mv",
                        filter.rule = "intersect",
                        minifrac.variable = 0.5,
                        minifrac.observation = 0.5) {
  # browser()
  data <- list(...)
  var.index <- list()
  obs.index <- list()

  for (i in 1:length(data)) {
    temp <- data[[i]]
    if (filter.item == "mv") {
      var.per <- apply(temp, 1, function(x) {
        sum(!is.na(x)) / ncol(temp)
      })
      obs.per <-
        apply(temp, 2, function(x) {
          sum(!is.na(x)) / nrow(temp)
        })
    }
    if (filter.item == "zero") {
      var.per <- apply(temp, 1, function(x) {
        sum(x != 0) / ncol(temp)
      })
      obs.per <-
        apply(temp, 2, function(x) {
          sum(x != 0) / nrow(temp)
        })
    }

    var.index[[i]] <- which(var.per >= minifrac.variable)
    obs.index[[i]] <- which(obs.per >= minifrac.observation)
  }

  if (filter.rule == "intersect") {
    var.index <- intersect.for.this.foo(var.index)
  }

  if (filter.rule == "union") {
    var.index <- union.for.this.foo(var.index)
  }

  for (i in 1:length(data)) {
    data[[i]] <- data[[i]][var.index, obs.index[[i]]]
  }

  return.data <- list(
    data = data,
    var.index = var.index,
    obs.index = obs.index,
    filter.item = filter.item,
    filter.rule = filter.rule,
    minifrac.variable = minifrac.variable,
    minifrac.observation = minifrac.observation
  )
  class(return.data) <- "SXTMinifracData"
  return(return.data)
}



## small functions for this function
intersect.for.this.foo <- function(x) {
  # browser()
  index <- x[[1]]
  for (i in 1:length(x)) {
    index <- intersect(index, x[[i]])
  }
  return(index)
}

union.for.this.foo <- function (x) {
  index <- x[[1]]
  for (i in 1:length(x)) {
    index <- union(index, x[[i]])
  }
  return(index)
}
