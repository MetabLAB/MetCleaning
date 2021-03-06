% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetabolitePlot.R
\name{MetabolitePlot}
\alias{MetabolitePlot}
\title{MetabolitePlot}
\usage{
MetabolitePlot(MetFlowData.before = MetFlowData1,
  MetFlowData.after = MetFlowData2, path = NULL, figure = "jpeg")
}
\arguments{
\item{MetFlowData.before}{MetFlowData before normalization or integration.}

\item{MetFlowData.after}{MetFlowData after normalization or integration.}

\item{path}{Work directory.}

\item{figure}{Figure type you want to draw. jpeg ot pdf, default is jpeg.}
}
\value{
Return the metabolite plot before and after processing.
}
\description{
Give scatter plot for each feature.
}
\author{
Xiaotao Shen
\email{shenxt@sioc.ac.cn}
}

