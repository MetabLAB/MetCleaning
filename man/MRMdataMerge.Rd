\name{MRMdataMerge}
\alias{MRMdataMerge}
\alias{MRMdataMerge}
\docType{package}
\title{
MRMdataMerge
}

\description{
Merge different MRM batch datasets into one dataset.
}

\usage{
MRMdataMerge(...,
             path =NULL,
             merge.method = "name",
             merge.term.index = c(1),
             tags.index = c(1:17),
             tags.same.in.datasets = c(1:4),
             tags.should.be.mean.in.datasets = c(11:13))
}
\arguments{
\tabular{ll}{
...: \tab The MRM data to be merged.\cr
path: \tab Directory for reading data and outputing results.\cr
merge.term.index: \tab Use which term to merge datasets?\cr
tags.index: \tab Which columna are tags in your data?\cr
tags.same.in.datasets: \tab Which columna are same in you datasets?And those information only remain one in merging.\cr
tags.should.be.mean.in.datasets: \tab Which columna should be averaged in the merging?\cr
}
}

\details{
\tabular{ll}{
Package: \tab MRMdataMerge\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-07-05\cr
License: \tab Nothing\cr
Depends: \tab R (> = 3.1.1) \cr
}
}

\value{
\tabular{ll}{
data: \tab It return a merged data.\cr
}
}

\author{
 Xiaotao Shen, Zhengjiang Zhu.

Maintainer: Xiaotao Shen <shenxt@sioc.ac.cn>
}

\references{
None.
}

\keyword{ package }
\seealso{
 \code{\link{none}}
}

\examples{
No.
}
