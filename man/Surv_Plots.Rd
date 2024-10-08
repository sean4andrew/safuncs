% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/safuncs.R
\name{Surv_Plots}
\alias{Surv_Plots}
\title{Generate Survival Plots}
\usage{
Surv_Plots(
  surv_db,
  plot_prefix = "plot_prefix",
  xlim = NULL,
  ylim = c(0, 1),
  xlab = "Days Post Challenge",
  lambda = NULL,
  phi = NULL,
  dailybin = TRUE,
  plot = "both",
  colours = NULL,
  theme = "ggplot",
  skip_error = FALSE,
  db = FALSE
)
}
\arguments{
\item{surv_db}{A survival dataframe as described in \bold{Details}.}

\item{plot_prefix}{A string specifying the prefix for the filename of the saved plots.}

\item{xlim}{A vector specifying the plots x-axis lower and upper limits, respectively.}

\item{ylim}{A vector specifying the Survival Plot y-axis lower and upper limits, respectively.}

\item{xlab}{A string specifying the plot x-axis label.}

\item{lambda}{Smoothing value for the hazard curve. Higher lambda produces greater smoothing. Defaults to NULL where \code{bshazard()} uses the provided survival data to estimate lambda; NULL specification is recommended for large sample size situations which usually occurs on our full-scale studies with many mortalities and tank-replication. At low sample sizes, the lambda estimate can be unreliable. Choosing a lambda of 10 (or anywhere between 1-100) probably produces the most accurate hazard curve for these situations. In place of choosing lambda, choosing \code{phi} is recommended; see below.}

\item{phi}{Dispersion parameter for the count model used in hazard curve estimation. Defaults to NULL where \code{bshazard()} uses the provided survival data to estimate phi; NULL specification is recommended for large sample size situations. At low sample sizes, the phi estimate can be unreliable. Choosing a phi value of 1 for low sample sizes is recommended. This value of 1 (or close) seems to be that estimated in past Tenaci data (QCATC997; phi ~ 0.8-1.4) where there are large sample sizes with tank-replication. The phi value of 1 indicates the set of counts (deaths) over time have a Poisson distribution, following the different hazard rates along the curve and are not overdispersed (phi > 1).}

\item{dailybin}{Whether to set time bins at daily (1 TTE) intervals. Refer to the \code{bshazard()} documentation for an understanding on the role of bins to hazard curve estimation. Please set to TRUE at low sample sizes and set to FALSE at large sample sizes with tank-replication. Defaults to TRUE.}

\item{plot}{Which plot to output. Use "surv" for the Kaplan-Meier Survival Curve, "haz" for the Hazard Curve, or "both" for both. Defaults to "both".}

\item{colours}{Vector of color codes for the different treatment groups in the plot. Defaults to ggplot2 default palette.}

\item{theme}{Character string specifying the graphics theme for the plots. Theme "ggplot2" and "prism" currently available. Defaults to "ggplot2".}

\item{skip_error}{Whether to skip the plotting of treatment groups with errors in hazard curve estimation from \code{bshazard()}. Defaults to FALSE.}

\item{db}{Whether to retun the hazard database. Defaults to FALSE.}
}
\value{
If \code{plot == "surv"}, returns a ggplot2 object reflecting the Kaplan-Meier Survival Curve.
If \code{plot == "haz"}, returns a ggplot2 object reflecting the Hazard Curve.
If \code{plot == "both"}, returns both ggplot2 objects in a list.
}
\description{
Produces a Kaplan-Meier Survival Plot and/or Hazard Plot from the given survival data. Each plot contains multiple curves for the different treatment groups. Saves plots automatically in the working directory.
}
\details{
The survival dataset should be a dataframe containing at least 4 different columns:
\itemize{
\item "Trt.ID" = Labels for treatment groups in the study.
\item "Tank.ID" = Labels for tanks in the study (each tank must have a unique label).
\item "TTE" = Time to Event. Event depends on "Status".
\item "Status" = Value indicating what happened at TTE. 1 for dead fish, 0 for survivors or those sampled and removed.
}

Each row should represent one fish. For an example dataframe, execute \code{data(surv_db_ex)} and view.

For details on the statistical methodology used by \code{bshazard()}, refer to: \url{https://www.researchgate.net/publication/287338889_bshazard_A_Flexible_Tool_for_Nonparametric_Smoothing_of_the_Hazard_Function}.

General concept: h(t) the hazard function is considered in an count model with the number of deaths as the response variable. I.e, death_count(t) = h(t) * P(t) where P(t) is the number alive as a function of time and h(t) is modeled over time using basis splines. The basis spline curvatures is assumed to have a normal distribution with mean 0 (a random effect). Based on this assumption, the author found that the variance of curvatures (i.e. smoothness) can be estimated as a function of over-dispersion (phi) of the death counts. Specifically, variance of curvatures = phi / lambda(smoothness parameter). Phi and lambda can be estimated from the data or specified by the user. Specification can be helpful in low sample size cases where estimates can be unreliable.
}
\examples{
data(surv_db_ex)
Surv_Plots(surv_db = surv_db_ex,
           plot_prefix = "QCATC777",
           xlim = c(0, 50),
           ylim = c(0, 1),
           xlab = "TTE",
           phi = 1,
           plot = "both")

}
