# Functions for the Independent Observations Model ----
# Víctor H Cervantes
# victorhc@illinois.edu
# vhcervantesb@unal.edu.co


## IO_RRTA() ----
#'  Computes the probability of rejecting the set when no target is shown (Rejection Rate) for the
#'  Independent Observations Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @importFrom stats pnorm
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::IO_RRTP(nSize = 2, criterion = .25)
#'           unforcedChoice:::IO_RRTP(nSize = 2, criterion = .50)
#'           unforcedChoice:::IO_RRTP(nSize = 2, criterion = -Inf)
#'
IO_RRTA <- function (criterion = .25,
                     nSize = 4, ...) {
    RR <- (pnorm(criterion))^(nSize)
    return(RR)
}


## IO_RRTP() ----
#'  Computes the probability of rejecting the set when a target is shown (Rejection Rate) for the
#'  Independent Observations Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTP()
#'
#' @importFrom stats pnorm
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::IO_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::IO_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::IO_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
IO_RRTP <- function (criterion = .25, mu = 1, sigmaT = 1,
                     nSize = 4, ...) {
    RR <- pnorm((criterion - mu) / sigmaT) * (pnorm(criterion))^(nSize - 1)
    return(RR)
}


## IO_HR() ----
#'  Computes the probability of choosing the target stimulus when a target is shown (Hit Rate) for the
#'  Independent Observations Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_HR()
#'
#' @importFrom stats dnorm pnorm integrate
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::IO_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::IO_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::IO_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
IO_HR <- function (criterion = .25, mu = 1, sigmaT = 1,
                   nSize = 4, ...) {
    if (criterion == Inf) {
        return(0)
    }

    likChoice <- function (x) {
        dnorm((x - mu) / sigmaT) * (pnorm(x))^(nSize - 1)
    }
    HR <- integrate(f = likChoice, lower = criterion, upper = Inf)$value
    return(HR)
}


## IO_FAR() ----
#'  Computes the probability of choosing the one non-target stimulus when no target is shown (False Alarm Rate) for the
#'  Independent Observations Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_FAR()
#'
#' @importFrom stats dnorm pnorm integrate
#'
#' @return numeric value of the False Alarm Rate
#'
#' @examples unforcedChoice:::IO_FAR(nSize = 2, criterion = .25)
#'           unforcedChoice:::IO_FAR(nSize = 2, criterion = .50)
#'           unforcedChoice:::IO_FAR(nSize = 2, criterion = -Inf)
#
IO_FAR <- function (criterion = .25,
                    nSize = 4, ...) {
    if (criterion == Inf) {
        return(0)
    }
    likChoice <- function (x) {
        dnorm(x) * (pnorm(x))^(nSize - 1)
    }
    FAR <- integrate(f = likChoice, lower = criterion, upper = Inf)$value
    return(FAR)
}
