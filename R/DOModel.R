# Functions for the Dependent Observations Model ----
# Víctor H Cervantes
# victorhc@illinois.edu
# vhcervantesb@unal.edu.co


## DO_RRTA() ----
#'  Computes the probability of rejecting the set when no target is shown (Rejection Rate) for the
#'  Dependent Observations Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @importFrom stats dnorm pnorm integrate
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::DO_RRTP(nSize = 2, criterion = .25)
#'           unforcedChoice:::DO_RRTP(nSize = 2, criterion = .50)
#'           unforcedChoice:::DO_RRTP(nSize = 2, criterion = -Inf)
#'
DO_RRTA <- function (criterion = .25, rho_0 = 0,
                     nSize = 4, ...) {

    # Parameter verification
    if (rho_0 == 0) {
        message("Independence assumed between all stimuli.\nIndependent Observations Model")
        return(IO_RRTA(criterion = criterion, nSize = nSize))
    }

    if (criterion == -Inf) {
        return(0)
    } else if (criterion == Inf) {
        return(1)
    }

    if (rho_0 == 1) {
        rho_0 <- .99999
    }

    # Auxiliary integrand function
    likChoice <- function (x) {
        dnorm(x) *
            (pnorm(
                (
                    criterion /
                        sqrt(1 - rho_0)
                ) +
                    (
                        x *
                            sqrt( rho_0 /
                                      (1 - rho_0)
                            )
                    )
            )
            )^(nSize)
    }

    RR <- integrate(f = likChoice, lower = -Inf, upper = Inf)$value
    return(RR)
}


## DO_RRTP() ----
#'  Computes the probability of rejecting the set when a target is shown (Rejection Rate) for the
#'  Dependent Observations Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTP()
#'
#' @importFrom stats dnorm pnorm integrate
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::DO_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::DO_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::DO_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
DO_RRTP <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                     nSize = 4, ...) {

    # Parameter verification
    if (rho_0 < (rho_1^2)) {
        rho_1 <- sign(rho_1) * sqrt(rho_0) * .9999
        warning("rho_0 must be at least as large as rho_1^2. Changing rho_1 to ", rho_1)
    }
    if (rho_0 == 0 && rho_1 == 0) {
        message("Independence assumed between all stimuli.\nIndependent Observations Model")
        return(IO_RRTP(criterion = criterion, mu = mu, sigmaT = sigmaT, nSize = nSize))
    }

    if (criterion == -Inf) {
        return(0)
    } else if (criterion == Inf) {
        return(1)
    }

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
                    (sqrt(rho_0) * (criterion - mu)) /
                        (sigmaT * sqrt(rho_0 - (rho_1^2)))
                ) +
                    (
                        x * rho_1 /
                            sqrt(
                                rho_0 - (rho_1^2)
                            )
                    )
            ) *
            (pnorm(
                (
                    criterion /
                        sqrt(1 - rho_0)
                ) +
                    (
                        x *
                            sqrt( rho_0 /
                                      (1 - rho_0)
                            )
                    )
            )
            )^(nSize - 1)
    }

    RR <- integrate(f = likChoice, lower = -Inf, upper = Inf)$value
    return(RR)
}


## DO_HR() ----
#'  Computes the probability of choosing the target stimulus when a target is shown (Hit Rate) for the
#'  Dependent Observations Model.
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
#' @examples unforcedChoice:::DO_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::DO_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::DO_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
DO_HR <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                   nSize = 4, ...) {

    # Parameter verification
    if (rho_0 < (rho_1^2)) {
        rho_1 <- sign(rho_1) * sqrt(rho_0) * .9999
        warning("rho_0 must be at least as large as rho_1^2. Changing rho_1 to ", rho_1)
    }
    if (rho_0 == 0 && rho_1 == 0) {
        message("Independence assumed between all stimuli.\nIndependent Observations Model")
        return(IO_HR(criterion = criterion, mu = mu, sigmaT = sigmaT, nSize = nSize))
    }

    if (criterion == Inf) {
        return(0)
    }

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
                    (sqrt((sigmaT^2) - (2 * sigmaT * rho_1) + rho_0) * (mu - criterion)) /
                        (sigmaT * sqrt(rho_0 - (rho_1^2)))
                ) +
                    (
                        x * (sigmaT - rho_1) /
                            sqrt(
                                rho_0 - (rho_1^2)
                            )
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


## DO_FAR() ----
#'  Computes the probability of choosing the one non-target stimulus when no target is shown (False Alarm Rate) for the
#'  Dependent Observations Model.
#'
#' @param criterion : Value of the criterion for selection
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_FAR()
#'
#' @importFrom stats dnorm pnorm integrate
#'
#' @return numeric value of the False Alarm Rate
#'
#' @examples unforcedChoice:::DO_FAR(nSize = 2, criterion = .25)
#'           unforcedChoice:::DO_FAR(nSize = 2, criterion = .50)
#'           unforcedChoice:::DO_FAR(nSize = 2, criterion = -Inf)
#
DO_FAR <- function (criterion = .25, rho_0 = 0,
                    nSize = 4, ...) {

    # Parameter verification
    if (rho_0 == 0) {
        message("Independence assumed between all stimuli.\nIndependent Observations Model")
        return(IO_FAR(criterion = criterion, nSize = nSize))
    }

    if (criterion == -Inf) {
        return(1 / nSize)
    } else if (criterion == Inf) {
        return(0)
    }


    # Auxiliary integrand function
    likChoice <- function (x) {
        dnorm(x) *
            pnorm(
                (
                    (- criterion) /
                        sqrt(rho_0)
                ) +
                    (
                        x * sqrt(1 - rho_0) /
                            sqrt(rho_0)
                    )
            ) *
            (
                pnorm(x)
            )^(nSize - 1)
    }

    FAR <- integrate(f = likChoice, lower = -Inf, upper = Inf)$value
    return(FAR)
}
