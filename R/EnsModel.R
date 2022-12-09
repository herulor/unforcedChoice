# Functions for the Ensemble Model ----
# Víctor H Cervantes
# victorhc@illinois.edu
# vhcervantesb@unal.edu.co


## Ens_RR() ----
#'  Computes the probability of rejecting the set when (Rejection Rate) for the
#'  Ensemble Model
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param method.Ens : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#' @param verbose : Logical. If TRUE print the nested integral to be computed when method.Ens == 'base'
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @importFrom stats integrate dnorm
#' @importFrom mvtnorm pmvnorm
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Ens_RR(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::Ens_RR(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::Ens_RR(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
Ens_RR <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                    nSize = 4, method.Ens = "base", verbose = FALSE, ...) {

    # Parameter verification
    if (criterion == 0) {
        return(0)
    }

    if (rho_0 < 0) {
        rho_0 <- 0
        warning("rho_0 must be non-negative.")
    }

    if (rho_0 < (rho_1^2)) {
        rho_1 <- sign(rho_1) * sqrt(rho_0) * .9999
        warning("rho_0 must be at least as large as rho_1^2. Changing rho_1 to ", rho_1)
    }

    if (rho_0 == 1) {
        if (abs(rho_1) == 1) {
            rho_1 <- sign(rho_1) * .99999
        }
        rho_0 <- .99999
    }

    if (nSize > 5 && method.Ens == "base") {
        message("Consider switching method.Ens to 'mvnorm'. 'base' method.Ens may be too slow.")
    }

    if (method.Ens == "base") {

        sigmaTE <- (nSize - 1) * (
            ((nSize - 1) * (
                (sigmaT^2) - (2 * sigmaT * rho_1) + rho_0
            )) + (
                1 - rho_0
            )
        )

        if (nSize > 2) {
            vectorZ <- paste0("z", seq(nSize - 2))
        } else {
            vectorZ <- NA
        }

        lowerLimits <- upperLimits <- character(nSize - 1)

        upperLimits[1] <- as.character((
            (nSize * criterion) - ((nSize - 1) * mu)
        ) / (
            sqrt(sigmaTE)
        )
        )

        lowerLimits[1] <- as.character(
            -(nSize - 1) * (
                (nSize * criterion) + mu
            ) / (
                sqrt(sigmaTE)
            )
        )

        if (nSize > 2) {
            upperLimits[2] <- paste0(as.character(
                sqrt( (nSize - 1) / (nSize - 2) )
            ), " * (",
            as.character(
                ( ((nSize * criterion) + mu) / (nSize * sqrt(1 - rho_0) ) )
            ),
            " + ",
            paste0(
                "z0 * ",
                as.character(
                    sqrt(sigmaTE) /
                        (
                            nSize * (nSize - 1) * sqrt(1 - rho_0)
                        )
                ),")"
            )
            )

            lowerLimits[2] <- paste0(
                as.character(
                    -(nSize - 2)
                ),
                " * ",
                upperLimits[2]
            )
        }

        if (nSize > 3) {

            for (kk in 3:(nSize - 1)) {
                upperLimits[kk] <- paste0(as.character(
                    sqrt( (nSize - kk + 1) / (nSize - kk) )
                ), " * (",
                as.character(
                    ( ((nSize * criterion) + mu) / (nSize * sqrt(1 - rho_0) ) )
                ),
                " + ",
                paste0(
                    "z0 * ",
                    as.character(
                        sqrt(sigmaTE) /
                            (
                                nSize * (nSize - 1) * sqrt(1 - rho_0)
                            )
                    )
                ),
                " + ",
                paste("(", paste0(vectorZ[seq(kk - 2)], " / ",
                                  as.character(
                                      sqrt( (nSize - seq(kk - 2) - 1) * (nSize - seq(kk - 2)) )
                                  ),
                                  " )"),
                      collapse = " + "),
                " )"
                )

                lowerLimits[kk] <- paste0(
                    as.character(
                        -(nSize - kk)
                    ),
                    " * ",
                    upperLimits[kk]
                )

            }

        }

        if (any(is.na(vectorZ))) {
            vectorZ <- "z0"
        } else {
            vectorZ <- c("z0", vectorZ)
        }

        integral <- "integrate("

        for (ii in seq_along(vectorZ)){
            integral <- paste0(
                integral,
                "lower = ", lowerLimits[ii], ", ",
                "upper = ", upperLimits[ii], ", ",
                "f = function (y) { sapply(y, function (", vectorZ[ii], ") { dnorm(", vectorZ[ii], ")",
                ifelse(ii == length(vectorZ), "", " * integrate(")
            )
        }

        for (ii in seq_along(vectorZ)){
            integral <- paste0(
                integral,
                ifelse(ii == 1, " } ) } )", "$value } ) } )")
            )
        }


        if (verbose) {
            print(integral)
        }
        probability <- eval(parse(text = integral))

        if (probability[["message"]] == "OK") {
            error <- probability[["abs.error"]]
            probability <- probability[["value"]]
            attr(probability, "error") <- error
            attr(probability, "msg") <- "Normal Completion"
        }

    } else if (method.Ens == "mvnorm") {

        aFC <- (sigmaT^2) - (2 * sigmaT * rho_1) + rho_0
        vectorMu <- mu * c((nSize - 1), rep(-1, nSize - 1))

        covMatrix <- cbind(c(
            (nSize - 1) * (
                (nSize - 1) * aFC + (1 - rho_0)
            ),
            -rep(
                (nSize - 1) * aFC + (1 - rho_0),
                (nSize - 1))
        ),
        rbind(
            -rep(
                (nSize - 1) * aFC + (1 - rho_0),
                (nSize - 1)),
            (
                (nSize^2) * (1 - rho_0) * diag(nrow = nSize - 1)
            ) - (
                ((nSize + 1) * (1 - rho_0))  - aFC
            )
        )
        )

        probability <- pmvnorm(lower = rep(-Inf, nSize),
                                        upper = rep(nSize * criterion, nSize),
                                        mean = vectorMu,
                                        sigma = covMatrix)

    } else if (method.Ens == "truncated") {
        stop("Approximation not implemented to compute rejection rates.")
    } else {
        stop("Unsupported method.Ens. method.Ens should be either 'base'  or 'mvnorm'.")
    }

    return(probability)
}


## Ens_RRTA() ----
#'  Computes the probability of rejecting the set when no target is shown (Rejection Rate) for the
#'  Ensemble Model
#'
#' @param criterion : Value of the criterion for selection
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param method.Ens : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @importFrom mvtnorm pmvnorm
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Ens_RRTA(nSize = 2, criterion = .25)
#'           unforcedChoice:::Ens_RRTA(nSize = 2, criterion = .50)
#'           unforcedChoice:::Ens_RRTA(nSize = 2, criterion = -Inf)
#'
Ens_RRTA <- function (criterion = .25, rho_0 = 0,
                      nSize = 4, method.Ens = "base", ...) {

    if (nSize == 1) {
        stop("Ensemble model does not apply to set size 1.")
    } else if (nSize > 1) {
        probability <- Ens_RR(criterion = criterion, mu = 0, sigmaT = 1, rho_0 = rho_0, rho_1 = rho_0,
                              nSize = nSize, method.Ens = method.Ens)
    }

    return(probability)
}


## Ens_RRTP() ----
#'  Computes the probability of rejecting the set when a target is shown (Rejection Rate) for the
#'  Ensemble Model
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param method.Ens : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTP()
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Ens_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::Ens_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::Ens_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
Ens_RRTP <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                      nSize = 4, method.Ens = "base", ...) {

    if (nSize == 1) {
        stop("Ensemble model does not apply to set size 1.")
    } else if (nSize > 1) {
        probability <- Ens_RR(criterion = criterion, mu = mu, sigmaT = sigmaT,
                              rho_0 = rho_0, rho_1 = rho_1,
                              nSize = nSize, method.Ens = method.Ens)
    }

    return(probability)
}


## Ens_ID() ----
#'  Computes the probability of selecting the arbitrary first stimulus for the
#'  Ensemble Model
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param method.Ens : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#' @param verbose : Logical. If TRUE print the nested integral to be computed when method.Ens == 'base'
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @importFrom stats integrate dnorm pnorm
#' @importFrom mvtnorm pmvnorm
#'
#' @return numeric value of the Hit Rate
#' @export
#'
#' @examples unforcedChoice:::Ens_RR(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::Ens_RR(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::Ens_RR(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
Ens_ID <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                    nSize = 4, method.Ens = "base", verbose = FALSE, ...) {

    # Parameter verification
    if (criterion == Inf) {
        return(0)
    }

    if (rho_0 < 0) {
        rho_0 <- 0
        warning("rho_0 must be non-negative.")
    }

    if (rho_0 < (rho_1^2)) {
        rho_1 <- sign(rho_1) * sqrt(rho_0) * .9999
        warning("rho_0 must be at least as large as rho_1^2. Changing rho_1 to ", rho_1)
    }

    if (rho_0 == 1) {
        if (abs(rho_1) == 1) {
            rho_1 <- sign(rho_1) * .99999
        }
        rho_0 <- .99999
    }

    if (nSize > 5 && method.Ens == "base") {
        message("Consider switching method.Ens to 'mvnorm'. 'base' method.Ens may be too slow.")
    }

    if (sigmaT )

    sigmaTE <- (nSize - 1) * (
        ((nSize - 1) * (
            (sigmaT^2) - (2 * sigmaT * rho_1) + rho_0
        )) + (
            1 - rho_0
        )
    )

    if (method.Ens == "base") {

        if (nSize > 2) {
            vectorZ <- paste0("z", seq(nSize - 2))
        } else {
            vectorZ <- NA
        }

        lowerLimits <- upperLimits <- character(nSize - 1)

        lowerLimits[1] <- as.character((
            (nSize * criterion) - ((nSize - 1) * mu)
        ) / (
            sqrt(sigmaTE)
        )
        )

        upperLimits[1] <- as.character(
            Inf
        )

        if (nSize > 2) {
            lowerLimits[2] <- paste0(
                "-((z0 * ", as.character(sqrt(sigmaTE)), ") + ",
                as.character(
                    (nSize - 1) * mu
                ), ") /",
                as.character(
                    sqrt((nSize - 1) * (nSize - 2) * (1 - rho_0))
                )
            )

            upperLimits[2] <- paste0(
                as.character(
                    -(nSize - 2)
                ),
                " * ",
                lowerLimits[2]
            )
        }

        if (nSize > 3) {

            for (kk in 3:(nSize - 1)) {
                lowerLimits[kk] <- paste0("-", as.character(
                    sqrt( (nSize - kk + 1) / (nSize - kk) )
                ), " * (",
                as.character(
                    (mu / sqrt(1 - rho_0))
                ),
                " + ",
                paste0(
                    "z0 * ",
                    as.character(
                        sqrt(sigmaTE) /
                            (
                                (nSize - 1) * sqrt(1 - rho_0)
                            )
                    )
                ),
                " - ",
                paste("(", paste0(vectorZ[seq(kk - 2)], " / ",
                                  as.character(
                                      sqrt( (nSize - seq(kk - 2) - 1) * (nSize - seq(kk - 2)) )
                                  ),
                                  " )"),
                      collapse = " + "),
                " )"
                )

                upperLimits[kk] <- paste0(
                    as.character(
                        -(nSize - kk)
                    ),
                    " * ",
                    lowerLimits[kk]
                )

            }

        }

        if (any(is.na(vectorZ))) {
            vectorZ <- "z0"
        } else {
            vectorZ <- c("z0", vectorZ)
        }

        integral <- "integrate("

        for (ii in seq_along(vectorZ)){
            integral <- paste0(
                integral,
                "lower = ", lowerLimits[ii], ", ",
                "upper = ", upperLimits[ii], ", ",
                "f = function (y) { sapply(y, function (", vectorZ[ii], ") { dnorm(", vectorZ[ii], ")",
                ifelse(ii == length(vectorZ), "", " * integrate(")
            )
        }

        for (ii in seq_along(vectorZ)){
            integral <- paste0(
                integral,
                ifelse(ii == 1, " } ) } )", "$value } ) } )")
            )
        }


        if (verbose) {
            print(integral)
        }
        probability <- eval(parse(text = integral))

        if (probability[["message"]] == "OK") {
            error <- probability[["abs.error"]]
            probability <- probability[["value"]]
            attr(probability, "error") <- error
            attr(probability, "msg") <- "Normal Completion"
        }

    } else if (method.Ens == "mvnorm") {

        aFC <- (sigmaT^2) - (2 * sigmaT * rho_1) + rho_0
        vectorMu <- mu * c((nSize - 1), rep(1, nSize - 1))

        covMatrix <- cbind(c(
            sigmaTE,
            rep(
                sigmaTE / (nSize - 1),
                (nSize - 1))
        ),
        rbind(
            rep(
                sigmaTE / (nSize - 1),
                (nSize - 1)),
            (
                (1 - rho_0) * diag(nrow = nSize - 1)
            ) + (
                aFC
            )
        )
        )

        probability <- pmvnorm(lower = c(nSize * criterion, rep(0, nSize - 1)),
                                        upper = rep(Inf, nSize),
                                        mean = vectorMu,
                                        sigma = covMatrix)

    } else if (method.Ens == "truncated") {
        if (rho_0 != rho_1) {
            stop("Approximation is only defined for equal rho_0 and rho_1")
        }
        if (sigmaT != 1) {
            stop("Approximation only implemented for equal variances")
        }

        rho <- rho_0
        sigmaT <- sigmaT * sqrt(1 - rho)

        ComputezConst <- function (x, mu, sigmaT, rho) {
            dnorm(scale(x, mu, sigmaT)) / pnorm(scale(x, mu, sigmaT))
        }

        muTrunc <- function (x, mu, sigmaT, rho) {
            mu - sigmaT * ComputezConst(x, mu, sigmaT)
        }

        sigmaTTrunc <- function (x, mu, sigmaT, rho) {
            zConst <- ComputezConst(x, mu, sigmaT)

            sigmaT * sqrt(1 - (zConst * scale(x, mu, sigmaT)) - (zConst^2))
        }

        probability <- integrate(lower = -Inf, upper = Inf,
                                        f = function (x) {
                                            ifelse(is.nan(sigmaTTrunc(x, mu, sigmaT, rho)) | is.na(sigmaTTrunc(x, mu, sigmaT, rho)),
                                                   0,
                                                   dnorm(x, mean = mu, sd = sigmaT) * (
                                                       (pnorm(x, mean = 0, sd = sqrt(1 - rho)))^(nSize - 1)
                                                   ) * (
                                                       1 - pnorm(criterion,
                                                                        mean =
                                                                            (1 - (1 / nSize)) *
                                                                            (x - muTrunc(x, mu = 0, sigmaT = sqrt(1 - rho))),
                                                                        sd = (sigmaTTrunc(x, mu = 0, sigmaT = sqrt(1 - rho)) * sqrt(nSize - 1)) / nSize
                                                       )
                                                   )
                                            )
                                        }
        )

        if (probability[["message"]] == "OK") {
            error <- probability[["abs.error"]]
            probability <- probability[["value"]]
            attr(probability, "error") <- error
            attr(probability, "msg") <- "Normal Completion"
        }
    } else {
        stop("Unsupported method.Ens. method.Ens should be either 'base', 'mvnorm', or 'truncated'.")
    }

    return(probability)
}


## Ens_HR() ----
#'  Computes the probability of choosing the target stimulus when a target is shown (Hit Rate) for the
#'  Ensemble Model
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param method.Ens : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_HR()
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Ens_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice:::Ens_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice:::Ens_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
Ens_HR <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                    nSize = 4, method.Ens = "base", ...) {

    if (nSize == 1) {
        stop("Ensemble model does not apply to set size 1.")
    } else if (nSize > 1) {
        probability <- Ens_ID(criterion = criterion, mu = mu, sigmaT = sigmaT,
                              rho_0 = rho_0, rho_1 = rho_1,
                              nSize = nSize, method.Ens = method.Ens)
    }

    return(probability)
}


## Ens_FAR() ----
#'  Computes the probability of choosing the one non-target stimulus when no target is shown (False Alarm Rate) for the
#'  Ensemble Model
#'
#' @param criterion : Value of the criterion for selection
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param method.Ens : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_FAR()
#'
#' @return numeric value of the False Alarm Rate
#' @export
#'
#' @examples unforcedChoice:::Ens_FAR(nSize = 2, criterion = .25)
#'           unforcedChoice:::Ens_FAR(nSize = 2, criterion = .50)
#'           unforcedChoice:::Ens_FAR(nSize = 2, criterion = -Inf)
#
Ens_FAR <- function (criterion = .25, rho_0 = 0,
                     nSize = 4, method.Ens = "base", ...) {

    if (nSize == 1) {
        stop("Ensemble model does not apply to set size 1.")
    } else if (nSize > 1) {
        probability <- Ens_ID(criterion = criterion, mu = 0, sigmaT = 1, rho_0 = rho_0, rho_1 = rho_0,
                              nSize = nSize, method.Ens = method.Ens)
    }

    return(probability)
}
