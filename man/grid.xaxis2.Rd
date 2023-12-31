% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grobFuncs.R
\name{grid.xaxis2}
\alias{grid.xaxis2}
\title{grid.xaxis2 function}
\usage{
grid.xaxis2(
  at = NULL,
  breaks = 5,
  labels = NULL,
  tick.len = 0.5,
  label.space = 0.5,
  side = c("bottom", "top"),
  rot = 0,
  ticks.gp = NULL,
  text.gp = NULL,
  hjust = NULL,
  vjust = NULL,
  draw = TRUE,
  name = NULL,
  gp = NULL,
  vp = NULL
)
}
\arguments{
\item{at}{numeric vector specifying the locations of tick marks on the axis.}

\item{breaks}{integer specifying the number of breaks on the axis.}

\item{labels}{character vector providing custom labels for the tick marks.}

\item{tick.len}{numeric specifying the length of the tick marks.}

\item{label.space}{numeric specifying the space between the tick marks and labels.}

\item{side}{character specifying the side to place the axis ("bottom" or "top").}

\item{rot}{numeric specifying the rotation angle for the tick mark labels.}

\item{ticks.gp}{gpar() seetings for the axis ticks.}

\item{text.gp}{gpar() seetings for the axis label.}

\item{hjust}{label hjust.}

\item{vjust}{label vjust.}

\item{draw}{if draw == "FALSE" , then return a gTree object.}

\item{name}{grob name.}

\item{gp}{gpar() seetings.}

\item{vp}{viewport seetings.}
}
\description{
grid.xaxis2 function
}
\author{
Jun Zhang

This function creates a x-axis in a graphical grid.
}
