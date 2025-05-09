% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/safuncs.R
\name{silencer}
\alias{silencer}
\title{Silence Code Output}
\usage{
silencer(x)
}
\arguments{
\item{x}{Code which output is to be directed to the sink.}
}
\value{
Code output without the directed outputs, e.g. output from \code{cat()}.
}
\description{
Hide output from R console by redirecting output using \code{sink(tempfile())} and subsequently \code{sink()}.
}
