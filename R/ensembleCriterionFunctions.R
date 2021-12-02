#'  rejection_ensembleCrit
#'  Computes the probability of rejecting the lineup
#'  ensemble criterion model.
#'  The rate is computed for the Target Absent or the Target Present case by setting
#'  the parameter mu equal to mu_0 or mu_1, the respective means of the
#'  Target Absent and Target Present distributions.
#'
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the marginal distributions
#' @param sigma : The common Standard deviation of the marginal distributions (Heteroscedasticity is not yet implemented)
#' @param rho : Value of the common correlation among stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#'
#' @return numeric value of the probability
#' @export
#' @importFrom mvtnorm pmvnorm
#'
#' @examples rejection_ensembleCrit(criterion = 1, mu = 0, sigma = 1, rho = 0, nSize = 2)
#'           rejection_ensembleCrit(criterion = 1, mu = 1, sigma = 1, rho = 0, nSize = 2)
#'
rejection_ensembleCrit <- function (criterion = .25, mu = 0, sigma = 1, rho = 0,
                                    nSize = 4, method = "base") {

    if (criterion == 0) {
        return(0)
    } else if (criterion == Inf) {
        return(1)
    }

    if (rho == 1) {
        if ((nSize * mu) < (nSize * criterion)) {
            return(1)
        } else {
            return(0)
        }
    }

    if (method == "base") {

        vectorZ <- paste0("z", seq(nSize - 1))

        lowerLimits <- upperLimits <- character(nSize - 1)

        upperLimits[1] <- as.character((
            (nSize * criterion) - ((nSize - 1) * mu)
        ) / (
            sigma * sqrt(
                nSize * (nSize - 1) * (1 - rho)
            )
        ))

        lowerLimits[1] <- as.character(
            -(nSize - 1) * (
                (nSize * criterion) + mu
            ) / (
                sigma * sqrt(
                    nSize * (nSize - 1) * (1 - rho)
                )
            ))

        if (nSize > 2) {

            for (kk in 2:(nSize - 1)) {
                upperLimits[kk] <- paste0(as.character(
                    sqrt( (nSize - kk + 1) / (nSize - kk) )
                ), " * (",
                as.character(
                    ( (nSize * criterion + mu ) / (nSize * sigma * sqrt(1 - rho) ) )
                ),
                " + ",
                paste("(", paste0(vectorZ[seq(kk - 1)], " / ",
                                  as.character(
                                      sqrt( (nSize - seq(kk - 1) + 1) * (nSize - seq(kk - 1)) )
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


        probability <- eval(parse(text = integral))


    } else if (method == "mvnorm") {

        vectorMu <- mu * c((nSize - 1), rep(-1, nSize - 1))

        covMatrix <- (nSize^2) * (sigma^2) * (1 - rho) *
            (diag(nSize) - (1 / nSize))

        probability <- mvtnorm::pmvnorm(lower = rep(-Inf, nSize),
                                        upper = rep(nSize * criterion, nSize),
                                        mean = vectorMu,
                                        sigma = covMatrix)

    } else {
        stop("Unsupported method. Method should be either 'base'  or 'mvnorm'. See ?rejection_ensembleCrit.")
    }

    return(probability)
}


#'  ID_ensembleCrit
#'  Computes the hit and false alarm rates for the lineup
#'  ensemble criterion model.
#'  The rate is computed for the False Alarm Rate or the Hit Rate by setting
#'  the parameter mu equal to mu_0 or mu_1, the respective means of the
#'  Target Absent and Target Present distributions.
#'
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the marginal distribution
#' @param mu_0 : Mean of the noise distributions
#' @param sigma : The common Standard deviation of the marginal distributions (Heteroscedasticity is not yet implemented)
#' @param rho : Value of the common correlation among stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster). Argument 'truncated' can also be used
#'                  to coompute Wixted et al. (2018) approximation (faster than 'base' but quite biased for small nSize).
#'
#' @return numeric value of the probability
#' @export
#' @importFrom mvtnorm pmvnorm
#'
#' @examples ID_ensembleCrit(criterion = 1, mu = 0, sigma = 1, rho = 0, nSize = 2)
#'           ID_ensembleCrit(criterion = 1, mu = 1, sigma = 1, rho = 0, nSize = 2)
#'
ID_ensembleCrit <- function (criterion = .25, mu = 0, mu_0 = 0, sigma = 1, rho = 0,
                             nSize = 4, method = "base") {

    # if (criterion == 0 & mu == 0) {
    #     return(1 / nSize)
    # } else if (criterion == Inf) {
    #     return(0)
    # }

    if (rho == 1) {
        if (mu >= criterion) {
            return(1)
        } else {
            return(0)
        }
    }

    if (method == "base") {

        vectorZ <- paste0("z", seq(nSize - 1))

        lowerLimits <- upperLimits <- character(nSize - 1)

        lowerLimits[1] <- nSize * criterion
        upperLimits[1] <- Inf

        if (nSize > 2) {

            lowerLimits[2] <- paste0("-z1 / ",
                                     as.character(
                                         sigma * sqrt(
                                             (nSize - 1) * (nSize - 2) * (1 - rho)
                                         )
                                     ))


            upperLimits[2] <- paste0("(",
                                     as.character((nSize - 2)),
                                     " * z1) / ",
                                     as.character(
                                         sigma * sqrt(
                                             (nSize - 1) * (nSize - 2) * (1 - rho)
                                         )
                                     ))
        }

        if (nSize > 3) {

            for (kk in 3:(nSize - 1)) {
                lowerLimits[kk] <- paste0(as.character(
                    sqrt( (nSize - kk + 1) / (nSize - kk) )
                ), " * ( (-z1 / ",
                as.character(
                    (nSize - 1)
                ),
                ") + ",
                paste("(", paste0(vectorZ[seq(from = 2, to = kk - 1)], " / ",
                                  as.character(
                                      sqrt( (nSize - seq(from = 2, to = kk - 1) + 1) * (nSize - seq(from = 2, to = kk - 1)) )
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

        integral <- "integrate("

        muVector    <- c((nSize - 1) * mu, rep(0, nSize - 2))
        sigmaVector <- c(sigma * sqrt(nSize * (nSize - 1) * (1 - rho)), rep(1, nSize - 2))

        for (ii in seq_along(vectorZ)){
            integral <- paste0(
                integral,
                "lower = ", lowerLimits[ii], ", ",
                "upper = ", upperLimits[ii], ", ",
                "f = function (y) { sapply(y, function (", vectorZ[ii],
                ") { dnorm(", vectorZ[ii], ", ",
                "mean = ",  muVector[ii], ", ",
                "sd = ",    sigmaVector[ii],
                ")",
                ifelse(ii == length(vectorZ), "", " * integrate(")
            )
        }

        for (ii in seq_along(vectorZ)){
            integral <- paste0(
                integral,
                ifelse(ii == 1, " } ) } )", "$value } ) } )")
            )
        }

        probability <- eval(parse(text = integral))

        if (probability[["message"]] == "OK") {
            error <- probability[["abs.error"]]
            probability <- probability[["value"]]
            attr(probability, "error") <- error
            attr(probability, "msg") <- "Normal Completion"
        }

    } else if (method == "mvnorm") {

        vectorMu <- mu * c((nSize - 1), rep(1, nSize - 1))

        covMatrix <- (sigma^2) * (1 - rho) *
            rbind(matrix(c(nSize * (nSize - 1), rep(nSize, nSize - 1)), nrow = 1),
                  cbind(matrix(rep(nSize, nSize - 1), ncol = 1),
                        (diag(nSize - 1) + 1)))

        probability <- mvtnorm::pmvnorm(lower = c(nSize * criterion, rep(0, nSize - 1)),
                                        upper = rep(Inf, nSize),
                                        mean = vectorMu,
                                        sigma = covMatrix)

    } else if (method == "truncated") {


        sigma <- sigma * sqrt(1 - rho)

        ComputezConst <- function (x, mu, sigma, rho) {
            dnorm(scale(x, mu, sigma)) / pnorm(scale(x, mu, sigma))
        }

        muTrunc <- function (x, mu, sigma, rho) {
            mu - sigma * ComputezConst(x, mu, sigma)
        }

        sigmaTrunc <- function (x, mu, sigma, rho) {
            zConst <- ComputezConst(x, mu, sigma)

            sigma * sqrt(1 - (zConst * scale(x, mu, sigma)) - (zConst^2))
        }

        probability <- integrate(lower = -Inf, upper = Inf,
                                 f = function (x) {
                                    ifelse(is.nan(sigmaTrunc(x, mu, sigma, rho)) | is.na(sigmaTrunc(x, mu, sigma, rho)),
                                           0,
                                           dnorm(x, mean = mu, sd = sigma) * (
                                               (pnorm(x, mean = mu_0, sd = sigma))^(nSize - 1)
                                           ) * (
                                               1 - pnorm(criterion,
                                                     mean =
                                                         (1 - (1 / nSize)) *
                                                         (x - muTrunc(x, mu = mu_0, sigma = sigma)),
                                                     sd = (sigmaTrunc(x, mu = mu_0, sigma = sigma) * sqrt(nSize - 1)) / nSize
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
        stop("Unsupported method. Method should be either 'base', 'mvnorm', or 'truncated'. See ?ID_ensembleCrit.")
    }

    return(probability)
}


#'  rejection_TA_ensembleCrit
#'  Computes the probability of rejecting the lineup when the target is absent for the
#'  ensemble criterion model.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu_0 : Mean of the marginal distributions of the non target stimuli
#' @param sigma_0 : The common Standard deviation of the marginal distributions (Heteroscedasticity is not yet implemented)
#' @param rho : Value of the common correlation among stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#'
#' @return numeric value of the probability
#' @export
#' @importFrom mvtnorm pmvnorm
#'
#' @examples rejection_TA_ensembleCrit(criterion = 1, mu_0 = 0, sigma_0 = 1, rho = 0, nSize = 2)
#'
rejection_TA_ensembleCrit <- function (criterion = .25, mu_0 = 0, sigma_0 = 1, rho = 0,
                                       nSize = 4, method = "base") {

    if (nSize == 1) {
        probability <- rejection_TA_projectionCrit(criterion = criterion, mu_0 = mu_0, sigma_0 = sigma_0, rho = rho,
                                                   nSize = nSize, method = method)
    } else if (nSize > 1) {
        probability <- rejection_ensembleCrit(criterion = criterion, mu = mu_0, sigma = sigma_0, rho = rho, nSize = nSize, method = method)
    }

    return(probability)
}


#'  rejection_TP_ensembleCrit
#'  Computes the probability of rejecting the lineup when the target is present for the
#'  ensemble criterion model.
#'
#' @param criterion : Value of the criterion for selection
#' @param mu_1 : Mean of the marginal distributions of the non target stimuli
#' @param sigma_0 : The common Standard deviation of the marginal distributions (Heteroscedasticity is not yet implemented)
#' @param rho : Value of the common correlation among stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#'
#' @return numeric value of the probability
#' @export
#' @importFrom mvtnorm pmvnorm
#'
#' @examples rejection_TP_ensembleCrit(criterion = 1, mu_1 = 1, sigma_0 = 1, rho = 0, nSize = 2)
#'
rejection_TP_ensembleCrit <- function (criterion = .25, mu_1 = 0, sigma_0 = 1, rho = 0,
                                       nSize = 4, method = "base") {

    if (nSize == 1) {
        probability <- rejection_TP_projectionCrit(criterion = criterion, mu_0 = mu_1, mu_1 = mu_1,
                                                   sigma_0 = sigma_0, sigma_1 = sigma_0, rho = rho,
                                                   nSize = nSize, method = method)
    } else if (nSize > 1) {
        probability <- rejection_ensembleCrit(criterion = criterion, mu = mu_1, sigma = sigma_0, rho = rho, nSize = nSize, method = method)
    }

    return(probability)
}


#'  target_ID_TP_ensembleCrit
#'  Computes the hit rate for the ensemble criterion model.
#'
#'
#' @param criterion : Value of the criterion for selection
#' @param mu_1 : Mean of the marginal distribution of the target stimulus
#' @param sigma_0 : The common Standard deviation of the marginal distributions (Heteroscedasticity is not yet implemented)
#' @param rho : Value of the common correlation among stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#'
#' @return numeric value of the probability
#' @export
#' @importFrom mvtnorm pmvnorm
#'
#' @examples target_ID_TP_ensembleCrit(criterion = 1, mu_1 = 1, sigma_0 = 1, rho = 0, nSize = 2)
#'
target_ID_TP_ensembleCrit <- function (criterion = .25, mu_1 = 0, sigma_0 = 1, rho = 0,
                                       nSize = 4, method = "base") {

    if (nSize == 1) {
        probability <- target_ID_TP_projectionCrit(criterion = criterion, mu_0 = mu_1, mu_1 = mu_1,
                                                   sigma_0 = sigma_0, sigma_1 = sigma_0, rho = rho, nSize = nSize, method = method)
    } else if (nSize > 1) {
        probability <- ID_ensembleCrit(criterion = criterion, mu = mu_1, sigma = sigma_0, rho = rho, nSize = nSize, method = method)
    }

    return(probability)
}


#'  filler_ID_TA_ensembleCrit
#'  Computes the false alarm rate for the ensemble criterion model.
#'
#'
#' @param criterion : Value of the criterion for selection
#' @param mu_0 : Mean of the marginal distribution of the target stimulus
#' @param sigma_0 : The common Standard deviation of the marginal distributions (Heteroscedasticity is not yet implemented)
#' @param rho : Value of the common correlation among stimuli
#' @param nSize : Number of simultaneously presented stimuli
#' @param method : Whether to integrate using base-R integrate (if set to "base" default. More precise, slower) or
#'                  mvtnorm::pmvnorm (if set to "mvnorm". Less precise, faster)
#'
#' @return numeric value of the probability
#' @export
#' @importFrom mvtnorm pmvnorm
#'
#' @examples filler_ID_TA_ensembleCrit(criterion = 1, mu_0 = 1, sigma_0 = 1, rho = 0, nSize = 2)
#'
filler_ID_TA_ensembleCrit <- function (criterion = .25, mu_0 = 0, sigma_0 = 1, rho = 0,
                                       nSize = 4, method = "base") {

    if (nSize == 1) {
        probability <- filler_ID_TA_projectionCrit(criterion = criterion, mu_0 = mu_0, sigma_0 = sigma_0, rho = rho,
                                                   nSize = nSize, method = method)
    } else if (nSize > 1) {
        probability <- ID_ensembleCrit(criterion = criterion, mu = mu_0, sigma = sigma_0, rho = rho, nSize = nSize, method = method)
    }

    return(probability)
}
