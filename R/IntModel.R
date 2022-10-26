# Functions for the Integration Model ----
# Víctor H Cervantes
# victorhc@illinois.edu
# vhcervantesb@unal.edu.co


## Int_RRTA() ----
#'  Computes the probability of rejecting the set when no target is shown (Rejection Rate) for the
#'  Integration Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @importFrom stats pnorm
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Int_RRTP(nSize = 2, criterion = .25)
#'           unforcedChoice:::Int_RRTP(nSize = 2, criterion = .50)
#'           unforcedChoice:::Int_RRTP(nSize = 2, criterion = -Inf)
#'
Int_RRTA <- function (criterion = .25, rho_0 = 0,
                      nSize = 4, ...) {

    if (criterion == -Inf) {
        return(0)
    } else if (criterion == Inf) {
        return(1)
    }

    RR <- pnorm(
        criterion /
            sqrt(nSize * (1 +
                              rho_0 * (nSize - 1)
            )
            )
    )
    return(RR)
}


## Int_RRTP() ----
#'  Computes the probability of rejecting the set when a target is shown (Rejection Rate) for the
#'  Integration Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTP()
#'
#' @importFrom stats pnorm
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Int_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::Int_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::Int_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
Int_RRTP <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                      nSize = 4, ...) {

    if (criterion == -Inf) {
        return(0)
    } else if (criterion == Inf) {
        return(1)
    }

    RR <- pnorm(
        (
            criterion - mu
        ) /
            (sqrt(
                (sigmaT^2) +
                    (
                        (nSize - 1) *
                            (
                                1 + (2 * sigmaT * rho_1) + ((nSize - 2) * rho_0)
                            )
                    )
            ))
    )
    return(RR)
}


## Int_HR() ----
#'  Computes the probability of choosing the target stimulus when a target is shown (Hit Rate) for the
#'  Integration Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_HR()
#'
#' @importFrom stats dnorm pnorm integrate
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Int_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::Int_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::Int_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
Int_HR <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                    nSize = 4, ...) {

    if (rho_0 == 1) {
        if (abs(rho_1) == 1) {
            rho_1 <- sign(rho_1) * .99999
        }
        rho_0 <- .99999
    }

    # Auxiliary integrand function
    likChoice <- function (x) {
        dnorm(x) *
            pnorm(
                (
                    (
                        sqrt(
                            (sigmaT^2) - (2 * sigmaT * rho_1) + rho_0
                        ) *
                            (mu - criterion)
                    ) + (
                        ((sigmaT^2) - ((nSize - 2) * (rho_0 - (sigmaT * rho_1))) - 1) * x
                    )
                ) /
                    sqrt(
                        (
                            (
                                (sigmaT^2) + (
                                    (nSize - 1) * (1 + 2 * sigmaT * rho_1 + (rho_0 * (nSize - 2)))
                                )
                            ) * (
                                (sigmaT^2) - (2 * sigmaT * rho_1) + rho_0
                            )
                        ) - (
                            (sigmaT^2) - ((nSize - 2) * (rho_0 - sigmaT * rho_1)) -1
                        )^2
                    )
            ) *
            (pnorm(
                (
                    mu /
                        sqrt(1 - rho_0)
                ) +
                    (
                        x *
                            sqrt( ((sigmaT^2) - (2 * sigmaT * rho_1) + rho_0) /
                                      (1 - rho_0)
                            )
                    )
            )
            )^(nSize - 1)
    }

    HR <- integrate(f = likChoice, lower = -Inf, upper = Inf)$value
    return(HR)
}


## Int_FAR() ----
#'  Computes the probability of choosing the one non-target stimulus when no target is shown (False Alarm Rate) for the
#'  Integration Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_FAR()
#'
#' @importFrom stats pnorm
#'
#' @return numeric value of the False Alarm Rate
#'
#' @examples unforcedChoice:::Int_FAR(nSize = 2, criterion = .25)
#'           unforcedChoice:::Int_FAR(nSize = 2, criterion = .50)
#'           unforcedChoice:::Int_FAR(nSize = 2, criterion = -Inf)
#
Int_FAR <- function (criterion = .25, rho_0 = 0,
                     nSize = 4, ...) {

    FAR <- (1  / nSize) * pnorm(
        -criterion /
            sqrt(
                nSize * (1 +
                             (
                                 (nSize - 1) * rho_0
                             )
                )
            )
    )
    return(FAR)
}
