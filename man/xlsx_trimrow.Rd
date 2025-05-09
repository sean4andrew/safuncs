% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/safuncs.R
\name{xlsx_trimrow}
\alias{xlsx_trimrow}
\title{Trim Rows Based on Non-NA Values}
\usage{
xlsx_trimrow(x, coli = 1)
}
\arguments{
\item{x}{A dataframe.}

\item{coli}{A number indicating the index of the column to base the trimming on.}
}
\value{
Returns a dataframe object without the "extra" NA values on the selected rows.
}
\description{
Remove rows after the last non-NA value in a selected column. Select column based on the \code{coli} argument.
}
