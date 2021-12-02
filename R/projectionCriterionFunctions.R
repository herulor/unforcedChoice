#' filler_ID_TA_projectionCrit
#'  Computes the probability of choosing one non-target stimulus when no target is shown in the lineup for the
#'  projection criterion models (Independent Observations, Wixted, and Dependent Observations, Melisa)
#'  This provides the False Alarm Rate.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu_0 : Mean of the target absent marginal distributions
#' @param sigma_0 : Standard deviation of the target absent marginal distributions
#' @param rho : Value of the common correlation among stimuli latent
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Integration method to pass to cubature::cubintegrate
#' @param ... : Additional (optional) parameters to pass to cubature::cubintegrate
#'
#' @return numeric value of the probability or list from cubature::cubintegrate if integration fails (non-zero returnCode)
#' @export
#' @importFrom stats dnorm pnorm
#' @importFrom cubature cubintegrate
#'
#' @examples filler_ID_TA_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = .25)
#'           filler_ID_TA_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = .50)
#'           filler_ID_TA_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = -Inf)
#'
filler_ID_TA_projectionCrit <- function (criterion = .25, mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 4, method = "cuhre", ...) {

    if (criterion == -Inf) {
        #        return(1)
        return(1 / nSize)
    } else if (criterion == Inf) {
        return(0)
    }

    filler_TA_integrand <- function (x) {
        (1 / sigma_0) *
            #nSize * (1 / sigma_0) *
            (dnorm(
                (x[1] - mu_0) / sigma_0
            )) *
            dnorm(x[2]) *
            (
                (pnorm(
                    (
                        (
                            (x[1] - mu_0) / sigma_0
                        ) * (
                            sqrt(1 - rho)
                        )
                    ) - (
                        (sqrt(rho) * x[2])
                    )
                )
                )^(nSize - 1)
            )
    }

    lowerLimits <- c(criterion, -Inf)
    upperLimits <- c(Inf, Inf)

    probability <- cubature::cubintegrate(f = filler_TA_integrand, lower = lowerLimits, upper = upperLimits,
                                          method = method, ...)

    if (probability[["returnCode"]] == 0) {
        probability <- probability[["integral"]]
    }

    return(probability)

}


#' target_ID_TP_projectionCrit
#'  Computes the probability of choosing the target stimulus when target is shown in the lineup for the
#'  projection criterion models (Independent Observations, Wixted, and Dependent Observations, Melisa)
#'  This provides the Hit Rate.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu_0 : Mean of the target absent marginal distributions
#' @param mu_1 : Mean of the target present marginal distributions
#' @param sigma_0 : Standard deviation of the target absent marginal distributions
#' @param sigma_1 : Standard deviation of the target present marginal distributions
#' @param rho : Value of the common correlation among stimuli latent
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Integration method to pass to cubature::cubintegrate
#' @param ... : Additional (optional) parameters to pass to cubature::cubintegrate
#'
#' @return numeric value of the probability or list from cubature::cubintegrate if integration fails (non-zero returnCode)
#' @export
#' @importFrom stats dnorm pnorm
#' @importFrom cubature cubintegrate
#'
#' @examples target_ID_TP_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = .25)
#'           target_ID_TP_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = .50)
#'           target_ID_TP_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = -Inf)
#'
target_ID_TP_projectionCrit <- function (criterion = .25, mu_0 = 0, mu_1 = 1, sigma_0 = 1, sigma_1 = 1, rho = 0,
                                         nSize = 4, method = "cuhre", ...) {

    if (criterion == Inf) {
        return(0)
    }

    if (sigma_1 != sigma_0) {
        stop("Heteroscedasticity is not yet implemented")
    } else {

        dPrime <- (mu_1 - mu_0) / sigma_0

        target_TP_integrand <- function (x) {
            (1 / sigma_1) *
                (dnorm(
                    (x[1] - mu_1) / sigma_1
                )
                ) *
                dnorm(x[2]) *
                ((
                    pnorm(
                        (
                            dPrime *
                                (
                                    sqrt( (rho^2) / (1 - rho) )
                                )
                        ) + (
                            (
                                (x[1] - mu_0) / sigma_0
                            ) *
                                sqrt(1 - rho)
                        ) - (
                            sqrt(rho) * x[2]
                        )
                    )
                )^(nSize - 1)
                )
        }

        lowerLimits <- c(criterion, -Inf)
        upperLimits <- c(Inf, Inf)

        probability <- cubature::cubintegrate(f = target_TP_integrand, lower = lowerLimits, upper = upperLimits,
                                              method = method, ...)

        if (probability[["returnCode"]] == 0) {
            probability <- probability[["integral"]]
        }

        return(probability)
    }
}


#' rejection_TP_projectionCrit
#'  Computes the probability of rejecting the lineup when the target is present for the
#'  projection criterion models (Independent Observations, Wixted, and Dependent Observations, Melisa)
#'
#' @param criterion : Value of the criterion for selection
#' @param mu_0 : Mean of the target absent marginal distributions
#' @param mu_1 : Mean of the target present marginal distributions
#' @param sigma_0 : Standard deviation of the target absent marginal distributions
#' @param sigma_1 : Standard deviation of the target present marginal distributions
#' @param rho : Value of the common correlation among stimuli latent
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Integration method to pass to cubature::cubintegrate
#' @param ... : Additional (optional) parameters to pass to cubature::cubintegrate
#'
#' @return numeric value of the probability or list from cubature::cubintegrate if integration fails (non-zero returnCode)
#' @export
#' @importFrom stats dnorm pnorm
#' @importFrom cubature cubintegrate
#'
#' @examples rejection_TP_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = .25)
#'           rejection_TP_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = .50)
#'           rejection_TP_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = -Inf)
#'
rejection_TP_projectionCrit <- function (criterion = .25, mu_0 = 0, mu_1 = 1, sigma_0 = 1, sigma_1 = 1, rho = 0,
                                         nSize = 4, method = "cuhre", ...) {

    if (criterion == -Inf) {
        return(0)
    } else if (criterion == Inf) {
        return(1)
    }

    if (sigma_1 != sigma_0) {
        stop("Heteroscedasticity is not yet implemented")
    } else {

        rejection_TP_integrand <- function (x) {
            dnorm(x) *
                pnorm(
                    (
                        (criterion - mu_1) /
                            (sigma_0 * sqrt(1 - rho))
                    ) +
                        (
                            x *
                                sqrt( rho /
                                          (1 - rho)
                                )
                        )
                ) *
                (pnorm(
                    (
                        (criterion - mu_0) /
                            (sigma_0 * sqrt(1 - rho))
                    ) +
                        (
                            x *
                                sqrt( rho /
                                          (1 - rho)
                                )
                        )
                )
                )^(nSize - 1)
        }

        lowerLimits <- -Inf
        upperLimits <-  Inf

        probability <- cubature::cubintegrate(f = rejection_TP_integrand, lower = lowerLimits, upper = upperLimits,
                                              method = method, ...)

        if (probability[["returnCode"]] == 0) {
            probability <- probability[["integral"]]
        }

        return(probability)
    }
}



#' rejection_TA_projectionCrit
#'  Computes the probability of rejecting the lineup when no target is present for the
#'  projection criterion models (Independent Observations, Wixted, and Dependent Observations, Melisa)
#'
#' @param criterion : Value of the criterion for selection
#' @param mu_0 : Mean of the target absent marginal distributions
#' @param sigma_0 : Standard deviation of the target absent marginal distributions
#' @param rho : Value of the common correlation among stimuli latent
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Integration method to pass to cubature::cubintegrate
#' @param ... : Additional (optional) parameters to pass to cubature::cubintegrate
#'
#' @return numeric value of the probability or list from cubature::cubintegrate if integration fails (non-zero returnCode)
#' @export
#' @importFrom stats dnorm pnorm
#' @importFrom cubature cubintegrate
#'
#' @examples rejection_TA_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = .25)
#'           rejection_TA_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = .50)
#'           rejection_TA_projectionCrit(mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2, criterion = -Inf)
#'
rejection_TA_projectionCrit <- function (criterion = .25, mu_0 = 0, sigma_0 = 1, rho = 0,
                                         nSize = 4, method = "cuhre", ...) {

    if (criterion == -Inf) {
        return(0)
    } else if (criterion == Inf) {
        return(1)
    }

    rejection_TA_integrand <- function (x) {
        dnorm(x) *
            (pnorm(
                (
                    (criterion - mu_0) /
                        (sigma_0 * sqrt(1 - rho))
                ) +
                    (
                        x *
                            sqrt( rho /
                                      (1 - rho)
                            )
                    )
            )
            )^(nSize)
    }

    lowerLimits <- -Inf
    upperLimits <-  Inf

    probability <- cubature::cubintegrate(f = rejection_TA_integrand, lower = lowerLimits, upper = upperLimits,
                                          method = method, ...)

    if (probability[["returnCode"]] == 0) {
        probability <- probability[["integral"]]
    }

    return(probability)
}
