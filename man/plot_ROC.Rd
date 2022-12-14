% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rocFunctions.R
\name{plot_ROC}
\alias{plot_ROC}
\title{Plots the ROC of the indicated type (TOC, IOC, EOC) for the unforced choice model.}
\usage{
plot_ROC(
  data = NULL,
  ratesTable = NULL,
  type = "TOC",
  model = "IO",
  nSize = 4,
  mu = 1,
  sigmaT = 1,
  rho_0 = 0,
  rho_1 = 0,
  criterion = NULL,
  nSteps = 41,
  breaks = NULL,
  method.Ens = NULL,
  AUC = FALSE,
  AUCtext = "",
  title = NULL,
  subtitle = NULL,
  color = NULL,
  colorData = NULL,
  alpha = 0.3,
  xlabText = NULL,
  ylabText = NULL,
  linetype = "solid",
  theme = TRUE,
  ...
)
}
\arguments{
\item{data}{: Optional data.frame with columns FAR, HR, RRTA, and RRTP of empirical points to add to the plot.}

\item{ratesTable}{: Optional data.frame with the structure of the output of Compute_Rates. If provided the curve and shading in the figure will used this table.}

\item{type}{: String indicating the type of ROC for which the AUC will be computed. Default "TOC"}

\item{model}{: String indicating the model for which the Hit Rate will be computed. Default "IO"}

\item{nSize}{: Number of simultaneously presented stimuli}

\item{mu}{: Mean of the target marginal distributions}

\item{sigmaT}{: Standard deviation of the target marginal distributions}

\item{rho_0}{: Value of the common correlation among non-target stimuli}

\item{rho_1}{: Value of the common correlation among non-target stimuli and the target}

\item{criterion}{: Optional. Values of the criterion at which the rates will be computed.
If missing, breaks or nSteps will be used to determine the criteria.}

\item{nSteps}{: Used if neither criterion nor breaks is provided. criterion values are created as a vector
of length nSteps}

\item{breaks}{: Optional. Numeric values between 0 and 1. criterion  values are created as qnorm(breaks)}

\item{method.Ens}{: Optional. Used to specify the method.Ens for the Ensemble model.}

\item{AUC}{: Logical. Indicates whether to print or not AUCtext. Default is FALSE.}

\item{AUCtext}{: Optional character. If not ptovided, it will be paste0('AUC = ', Compute_AUC(.)) with the call to Compute_AUC using the parameters in the plot_ROC call.}

\item{title}{: Optional character. If not provided will be a text indicating the type of ROC plotted}

\item{subtitle}{: Optional character. If not provided it will show the model and parameter values.}

\item{color}{: Optional character to set shading color.}

\item{colorData}{: Optional character to set data points color.}

\item{alpha}{: Optional numeric value for alpha transparency.}

\item{xlabText}{: Optional character. If not provided it will be 'FAR + FIR[TA]' or 'False Alarm Rate' depending on type of ROC}

\item{ylabText}{: Optional character. If not provided it will be 'HR + FIR[TP]' or 'Hit Rate' depending on type of ROC}

\item{linetype}{: Optional character with the linetype for the ROC curve. See ?ggplot2::geom_line}

\item{theme}{: Logical. If TRUE (Default) theme_minimal will be used for the plot. See ?ggplot2::theme_minimal}

\item{...}{: Additional parameters to be passed to the specific model function}
}
\value{
Object of class ggplot. Plot of the ROC.
}
\description{
Plots the ROC of the indicated type (TOC, IOC, EOC) for the unforced choice model.
}
\examples{
plot_ROC(mu = 1, sigmaT = 1, nSize = 2, AUC = TRUE, type = "TOC")

}
