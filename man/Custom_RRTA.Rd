% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CustomModel.R
\name{Custom_RRTA}
\alias{Custom_RRTA}
\title{Computes the probability of rejecting the set when no target is shown (Rejection Rate) for
 Custom Models}
\usage{
Custom_RRTA(criterion = 0.25, rho_0 = 0, nSize = 4, cMatrix, ...)
}
\arguments{
\item{criterion}{: Value of the criterion for selection}

\item{rho_0}{: Value of the common correlation among non-target stimuli}

\item{nSize}{: Number of simultaneously presented stimuli. Must equal ncol(cMatrix).}

\item{cMatrix}{: Matrix that defines the custom model}

\item{...}{: Captures unused parameter when called from generic function unforcedChoice_RRTA()}
}
\value{
numeric value of the Hit Rate
}
\description{
Computes the probability of rejecting the set when no target is shown (Rejection Rate) for
 Custom Models
}
\examples{
unforcedChoice:::Custom_RRTA(nSize = 2, criterion = .25,
                                       cMatrix = diag(2))
          unforcedChoice:::Custom_RRTA(nSize = 2, criterion = .50,
                                       cMatrix = diag(2))
          unforcedChoice:::Custom_RRTA(nSize = 2, criterion = -Inf,
                                       cMatrix = diag(2))

}
