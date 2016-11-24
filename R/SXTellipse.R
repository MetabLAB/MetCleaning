#' @title SXTellipse
#' @description Give a SXTellipseData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param parameters All the parameters can see help of ellipse.
#' @return Return a SXTellipseData
#' @seealso \code{\link[ellipse]{ellipse}}

SXTellipse <- function (x,
                        scale = c(1, 1),
                        centre = c(0, 0),
                        level = 0.95,
                        t = sqrt(qchisq(level, 2)),
                        which = c(1, 2),
                        npoints = 100,
                        ...) {
  # browser()
  names <- c("x", "y")
  if (is.matrix(x)) {
    xind <- which[1]
    yind <- which[2]
    r <- x[xind, yind]
    if (missing(scale)) {
      scale <- sqrt(c(x[xind, xind], x[yind, yind]))
      if (scale[1] > 0)
        r <- r / scale[1]
      if (scale[2] > 0)
        r <- r / scale[2]
    }
    if (!is.null(dimnames(x)[[1]]))
      names <- dimnames(x)[[1]][c(xind, yind)]
  }
  else
    r <- x
  r <-
    min(max(r, -1), 1)  # clamp to -1..1, in case of rounding errors
  d <- acos(r)
  a <- seq(0, 2 * pi, len = npoints)
  ellipse.data = matrix(
    c(t * scale[1] * cos(a + d / 2) + centre[1], t * scale[2] *
        cos(a - d / 2) + centre[2]),
    npoints,
    2,
    dimnames = list(NULL,
                    names)
  )
  SXTellipseData <- list(
    ellipse.data = ellipse.data,
    t = t,
    a = a,
    d = d,
    scale = scale,
    centre = centre
  )
  class(SXTellipseData) = "SXTellipseData"
  return(SXTellipseData)
}