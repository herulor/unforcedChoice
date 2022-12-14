% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IntModel.R
\name{Int_RRTA}
\alias{Int_RRTA}
\title{Computes the probability of rejecting the set when no target is shown (Rejection Rate) for the
 Integration Model.}
\usage{
Int_RRTA(criterion = 0.25, rho_0 = 0, nSize = 4, ...)
}
\arguments{
\item{criterion}{: Value of the criterion for selection}

\item{rho_0}{: Value of the common correlation among non-target stimuli}

\item{nSize}{: Number of simultaneously presented stimuli}

\item{...}{: Captures unused parameter when called from generic function unforcedChoice_RRTA()}
}
\value{
numeric value of the Hit Rate
}
\description{
Computes the probability of rejecting the set when no target is shown (Rejection Rate) for the
 Integration Model.
}
\examples{
unforcedChoice:::Int_RRTP(nSize = 2, criterion = .25)
          unforcedChoice:::Int_RRTP(nSize = 2, criterion = .50)
          unforcedChoice:::Int_RRTP(nSize = 2, criterion = -Inf)

}
