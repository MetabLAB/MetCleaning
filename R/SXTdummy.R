SXTdummy <- function (Y)
{
  dummy <- matrix(0, nrow = length(Y), ncol = length(table(Y)))
  for (i in 1:length(Y)) {
    for (j in 1:ncol(dummy)) {
      if (Y[i] == names(table(Y))[j])
        dummy[i, j] = 1
    }
  }
  return(dummy)
}