# Functions to estimate model parameters from data for several unforced choice SDT models ----
# Víctor H Cervantes
# victorhc@illinois.edu
# vhcervantesb@unal.edu.co


## parsInData() ----
#' Plugs the values of pars in the data data.frame
#'
#' @param pars : Numeric vector. Should contain the parameter values to assign to the unique values of
#'               criterion, mu, sigmaT, rho_0, and rho_1, excluding those that appear on 'fixed'.
#' @param fixed : Optional character vector. Name of parameter that will fixed at values in data data.frame. Can include any of
#'                 'criterion', 'mu', 'sigmaT', 'rho_0', 'rho_1'
#' @param data : data.frame. Must have columns:
#'                 'model': character. Can be any of 'IO', 'DO', 'Ens', 'Int', 'Custom'
#'                 'nTP': numeric. Number of target present trials for the row.
#'                 'nTA': numeric. Number of target absent trials for the row.
#'                 'nSize': numeric. Number of stimuli presented on each trial for the row.
#'                 'criterion': If numeric and fixed includes 'criterion', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a criterion to be estimated equal across rows with the same value.
#'                 'mu': If numeric and fixed includes 'mu', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a mu to be estimated equal across rows with the same value.
#'                 'sigmaT': If numeric and fixed includes 'sigmaT', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a sigmaT to be estimated equal across rows with the same value.
#'                 'rho_0': If numeric and fixed includes 'rho_0', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_0 to be estimated equal across rows with the same value.
#'                 'rho_1': If numeric and fixed includes 'rho_1', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_1 to be estimated equal across rows with the same value.
#'                 'HR' : Numeric. Either a value between 0 and 1 given the Hit Rate or an integer giving the number of Hits for the row.
#'                 'FAR' : Numeric. Either a value between 0 and 1 given the False Alarm Rate or an integer giving the number of Misses for the row.
#'                 'RRTP' : Numeric. Either a value between 0 and 1 given the Rejection Rate on target present trials or an integer giving the number of Incorrect Rejections for the row.
#'                 'RRTA' : Numeric. Either a value between 0 and 1 given the Rejection Rate on target absent or an integer giving the number of Correct Rejections for the row.
#'
#' @return data.frame with pars replacing certain columns of data. See usage.
#'
parsInData <- function (pars, fixed = NULL, data) {

    parameterNames <- c("criterion", "mu", "sigmaT", "rho_0", "rho_1")
    columns <- c("model", "nTP", "nTA", "nSize", parameterNames, "HR", "FAR", "RRTP", "RRTA")
    if (!all(columns %in% names(data))) {
        stop("data must be a data.frame with names '", paste(columns, collapse = ", "), "'")
    }

    if (!is.null(fixed) && all(parameterNames %in% fixed)) {
        message("All parameters fixed. Nothing done.")
    } else {

        if (!("criterion" %in% fixed)) {
            ncriterion <- table(data[['criterion']])

            if (length(pars) < length(ncriterion)) {
                stop("More values need replacing than available in pars. Check pars and fixed.")
            }

            criterion  <- pars[1:length(ncriterion)]
            valuescriterion <- rep(criterion, times = ncriterion)

            data[['criterion']] <- valuescriterion
            pars <- pars[-seq(length(ncriterion))]
        }

        if (!("mu" %in% fixed)) {
            nMeans <- table(data[['mu']])

            if (length(pars) < length(nMeans)) {
                stop("More values need replacing than available in pars. Check pars and fixed.")
            }

            means  <- pars[1:length(nMeans)]
            valuesMeans <- rep(means, times = nMeans)

            data[['mu']] <- valuesMeans
            pars <- pars[-seq(length(nMeans))]
        }

        if (!("sigmaT" %in% fixed)) {
            nSigmas <- table(data[['sigmaT']])

            if (length(pars) < length(nSigmas)) {
                stop("More values need replacing than available in pars. Check pars and fixed.")
            }

            sigmas  <- pars[1:length(nSigmas)]
            valuesSigmas <- rep(sigmas, times = nSigmas)

            data[['sigmaT']] <- valuesSigmas
            pars <- pars[-seq(length(nSigmas))]
        }

        if (!("rho_0" %in% fixed)) {
            nRho0s <- table(data[['rho_0']])

            if (length(pars) < length(nRho0s)) {
                stop("More values need replacing than available in pars. Check pars and fixed.")
            }

            rho0s  <- pars[1:length(nRho0s)]
            valuesRho0s <- rep(rho0s, times = nRho0s)

            data[['rho_0']] <- valuesRho0s
            pars <- pars[-seq(length(nRho0s))]
        }

        if (("r1=r0" %in% fixed) || ("r1=sr0" %in% fixed)) {

            if ("rho_1" %in% fixed) {
                message("rho_1 is set to be a function of rho_0. Ignoring values in data")
            }

            if (("r1=r0" %in% fixed) & ("r1=sr0" %in% fixed)) {
                message("Only one 'r1=r0' or 'r1=sr0' can be included in fixed'. Choose one.")
            }

            if ("r1=r0" %in% fixed) {
                data[['rho_1']] <- data[['rho_0']]
            }

            if ("r1=sr0" %in% fixed) {
                data[['rho_1']] <- data[['sigmaT']] * data[['rho_0']]
            }

        } else {
            if (!("rho_1" %in% fixed)) {
                nRho1s <- table(data[['rho_1']])

                if (length(pars) < length(nRho1s)) {
                    stop("More values need replacing than available in pars. Check pars and fixed.")
                }

                rho1s  <- pars[1:length(nRho1s)]
                valuesRho1s <- rep(rho1s, times = nRho1s)

                data[['rho_1']] <- valuesRho1s
                pars <- pars[-seq(length(nRho1s))]
            }
        }
    }
    if (length(pars) > 0) {
        warning("pars contains more values than values need replacement. Check that pars and fixed are correct")
    }

    return(data)
}


## predictData() ----
#' Predicts the HR, FAR, RRTP, and RRTA values the data data.frame (replaces those values is data if existing)
#'
#' @param data : data.frame. Must have columns:
#'                 'model': character. Can be any of 'IO', 'DO', 'Ens', 'Int', 'Custom'
#'                 'nTP': numeric. Number of target present trials for the row.
#'                 'nTA': numeric. Number of target absent trials for the row.
#'                 'nSize': numeric. Number of stimuli presented on each trial for the row.
#'                 'criterion': If numeric and fixed includes 'criterion', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a criterion to be estimated equal across rows with the same value.
#'                 'mu': If numeric and fixed includes 'mu', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a mu to be estimated equal across rows with the same value.
#'                 'sigmaT': If numeric and fixed includes 'sigmaT', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a sigmaT to be estimated equal across rows with the same value.
#'                 'rho_0': If numeric and fixed includes 'rho_0', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_0 to be estimated equal across rows with the same value.
#'                 'rho_1': If numeric and fixed includes 'rho_1', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_1 to be estimated equal across rows with the same value.
#' @param ... : Additional arguments passed to Compute_Rates
#'
#' @return data.frame with pars replacing certain columns of data. See usage.
#'
predictData <- function (data, ...) {

    parameterNames <- c("criterion", "mu", "sigmaT", "rho_0", "rho_1")
    columns <- c("model", "nTP", "nTA", "nSize", parameterNames)
    if (!all(columns %in% names(data))) {
        stop("data must be a data.frame with names '", paste(columns, collapse = ", "), "'")
    }

    data[, c("HR", "FAR", "RRTP", "RRTA")] <-
        matrix(unlist(apply(data[, c("model",
                                     "nSize",
                                     "mu",
                                     "sigmaT",
                                     "rho_0",
                                     "rho_1",
                                     "criterion")], 1,
                            function(x) suppressMessages(Compute_Rates(model    = as.character(x[1]),
                                                                       nSize    = as.numeric(x[2]),
                                                                       mu       = as.numeric(x[3]),
                                                                       sigmaT   = as.numeric(x[4]),
                                                                       rho_0    = as.numeric(x[5]),
                                                                       rho_1    = as.numeric(x[6]),
                                                                       criterion = as.numeric(x[7]),
                                                                       ... = ...))[, -1]
        )
        ), ncol = 4, byrow = TRUE)

    data[, "FIRTP"] <- 1 - data[, "HR"]  - data[, "RRTP"]
    data[, "FIRTA"] <- 1 - data[, "FAR"] - data[, "RRTA"]

    data[, c("HR", "FAR", "RRTP", "RRTA", "FIRTP", "FIRTA")][data[, c("HR", "FAR", "RRTP", "RRTA", "FIRTP", "FIRTA")] < 0] <- 0
    data[, c("HR", "FAR", "RRTP", "RRTA", "FIRTP", "FIRTA")][data[, c("HR", "FAR", "RRTP", "RRTA", "FIRTP", "FIRTA")] > 1] <- 1

    return(data)
}


## unforcedChoice_wls() ----
#' Computes the weighed least squares objective for the unforced choice SDT models.
#'
#' @param pars : Numeric vector with the values of the non-fixed parameters
#' @param data : data.frame. Must have columns:
#'                 'model': character. Can be any of 'IO', 'DO', 'Ens', 'Int', 'Custom'
#'                 'nTP': numeric. Number of target present trials for the row.
#'                 'nTA': numeric. Number of target absent trials for the row.
#'                 'nSize': numeric. Number of stimuli presented on each trial for the row.
#'                 'criterion': If numeric and fixed includes 'criterion', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a criterion to be estimated equal across rows with the same value.
#'                 'mu': If numeric and fixed includes 'mu', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a mu to be estimated equal across rows with the same value.
#'                 'sigmaT': If numeric and fixed includes 'sigmaT', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a sigmaT to be estimated equal across rows with the same value.
#'                 'rho_0': If numeric and fixed includes 'rho_0', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_0 to be estimated equal across rows with the same value.
#'                 'rho_1': If numeric and fixed includes 'rho_1', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_1 to be estimated equal across rows with the same value.
#'                 'HR' : Numeric. Either a value between 0 and 1 given the Hit Rate or an integer giving the number of Hits for the row.
#'                 'FAR' : Numeric. Either a value between 0 and 1 given the False Alarm Rate or an integer giving the number of Misses for the row.
#'                 'RRTP' : Numeric. Either a value between 0 and 1 given the Rejection Rate on target present trials or an integer giving the number of Incorrect Rejections for the row.
#'                 'RRTA' : Numeric. Either a value between 0 and 1 given the Rejection Rate on target absent or an integer giving the number of Correct Rejections for the row.
#' @param fixed : Optional character vector. Name of parameter that will fixed at values in data data.frame. Can include any of
#'                 'criterion', 'mu', 'sigmaT', 'rho_0', 'rho_1'
#' @param rates : Logical. Indicates whether columns HR, FAR, RRTP, and RRTA of data ought to be treated as proportions/rates (default) or as counts.
#' @param ... : additional parameters to be passed to the model specific functions.
#'
#' @return numeric value of the wls objective computed on data for pars
#'
unforcedChoice_wls <- function (pars, data, fixed = NULL, rates = TRUE, ...) {

    if (any(data[['sigmaT']] <= 0)) {
        stop("sigmaT must be positive")
    }

    if (any(data[['sigmaT']] <= 0)) {
        stop("sigmaT must be positive")
    }

    if (any((data[['rho_0']] - (data[['rho_1']])^2) < 0)) {
        stop("rho_0 must be at least as large as rho_1^2")
    }

    if (!rates) {
        data[, 'HR']   <- data[, 'HR'] / data[, 'nTP']
        data[, 'FAR']  <- data[, 'FAR'] / data[, 'nTA']
        data[, 'RRTP'] <- data[, 'RRTP'] / data[, 'nTP']
        data[, 'RRTA'] <- data[, 'RRTA'] / data[, 'nTA']
    }

    data[, "FIRTP"] <- 1 - data[, "HR"]  - data[, "RRTP"]
    data[, "FIRTA"] <- 1 - data[, "FAR"] - data[, "RRTA"]

    dataPred <- predictData(parsInData(pars = pars, fixed = fixed, data = data), ...)

    varPred <- dataPred
    varPred <-
        dataPred[, c("HR", "FAR", "RRTP", "RRTA", "FIRTP", "FIRTA")] *
        (1 - dataPred[, c("HR", "FAR", "RRTP", "RRTA", "FIRTP", "FIRTA")]
        ) /
        cbind(
            data[, c("nTP", "nTA")],
            data[, c("nTP", "nTA")],
            data[, c("nTP", "nTA")]
        )

    num <- (data[, c("HR", "FAR", "RRTP", "RRTA", "FIRTP", "FIRTA")] -
                dataPred[, c("HR", "FAR", "RRTP", "RRTA", "FIRTP", "FIRTA")])^2

    varPred[varPred < .Machine$double.eps] <- .00001

    obj <- num / varPred

    return(sum(obj))
}


## unforcedChoice_lik() ----
#' Computes the likelihood for the unforced choice SDT models.
#'
#' @param pars : Numeric vector with the values of the non-fixed parameters
#' @param data : data.frame. Must have columns:
#'                 'model': character. Can be any of 'IO', 'DO', 'Ens', 'Int', 'Custom'
#'                 'nTP': numeric. Number of target present trials for the row.
#'                 'nTA': numeric. Number of target absent trials for the row.
#'                 'nSize': numeric. Number of stimuli presented on each trial for the row.
#'                 'criterion': If numeric and fixed includes 'criterion', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a criterion to be estimated equal across rows with the same value.
#'                 'mu': If numeric and fixed includes 'mu', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a mu to be estimated equal across rows with the same value.
#'                 'sigmaT': If numeric and fixed includes 'sigmaT', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a sigmaT to be estimated equal across rows with the same value.
#'                 'rho_0': If numeric and fixed includes 'rho_0', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_0 to be estimated equal across rows with the same value.
#'                 'rho_1': If numeric and fixed includes 'rho_1', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_1 to be estimated equal across rows with the same value.
#'                 'HR' : Numeric. Either a value between 0 and 1 given the Hit Rate or an integer giving the number of Hits for the row.
#'                 'FAR' : Numeric. Either a value between 0 and 1 given the False Alarm Rate or an integer giving the number of Misses for the row.
#'                 'RRTP' : Numeric. Either a value between 0 and 1 given the Rejection Rate on target present trials or an integer giving the number of Incorrect Rejections for the row.
#'                 'RRTA' : Numeric. Either a value between 0 and 1 given the Rejection Rate on target absent or an integer giving the number of Correct Rejections for the row.
#' @param fixed : Optional character vector. Name of parameter that will fixed at values in data data.frame. Can include any of
#'                 'criterion', 'mu', 'sigmaT', 'rho_0', 'rho_1'
#' @param rates : Logical. Indicates whether columns HR, FAR, RRTP, and RRTA of data ought to be treated as proportions/rates (default) or as counts.
#' @param ... : additional parameters to be passed to the model specific functions.
#'
#' @importFrom stats dmultinom
#'
#' @return numeric value of -log-likelihood of data computed on pars
#'
unforcedChoice_lik <- function (pars, data, fixed = NULL, rates = TRUE, ...) {

    if (any(data[['sigmaT']] <= 0)) {
        stop("sigmaT must be positive")
    }

    if (any(data[['sigmaT']] <= 0)) {
        stop("sigmaT must be positive")
    }

    if (any((data[['rho_0']] - (data[['rho_1']])^2) < 0)) {
        stop("rho_0 must be at least as large as rho_1^2")
    }

    if (rates) {
        data[, 'HR']   <- data[, 'HR'] * data[, 'nTP']
        data[, 'FAR']  <- data[, 'FAR'] * data[, 'nTA']
        data[, 'RRTP'] <- data[, 'RRTP'] * data[, 'nTP']
        data[, 'RRTA'] <- data[, 'RRTA'] * data[, 'nTA']
    }

    data[, "FIRTP"] <- data[, "nTP"] - data[, "HR"]  - data[, "RRTP"]
    data[, "FIRTA"] <- data[, "nTA"] - data[, "FAR"] - data[, "RRTA"]

    dataPred <- predictData(parsInData(pars = pars, fixed = fixed, data = data), ...)

    obj <- numeric(nrow(data))

    for (ii in seq(obj)) {
        obj[ii] <- dmultinom(x    = as.numeric(data[ii, c("HR", "RRTP", "FIRTP")]),
                             prob = as.numeric(dataPred[ii, c("HR", "RRTP", "FIRTP")]),
                             log = TRUE
        ) +
            dmultinom(x    = as.numeric(data[ii, c("FAR", "RRTA", "FIRTA")]),
                      prob = as.numeric(dataPred[ii, c("FAR", "RRTA", "FIRTA")]),
                      log = TRUE
            )
    }

    obj[obj == -Inf] <- -5000

    return(-sum(obj))
}


## unforcedChoice_Estimation() ----
#' Estimates the parameters for unforced choice SDT models.
#'
#' @param data : data.frame. Must have columns:
#'                 'model': character. Can be any of 'IO', 'DO', 'Ens', 'Int', 'Custom'
#'                 'nTP': numeric. Number of target present trials for the row.
#'                 'nTA': numeric. Number of target absent trials for the row.
#'                 'nSize': numeric. Number of stimuli presented on each trial for the row.
#'                 'criterion': If numeric and fixed includes 'criterion', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a criterion to be estimated equal across rows with the same value.
#'                 'mu': If numeric and fixed includes 'mu', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a mu to be estimated equal across rows with the same value.
#'                 'sigmaT': If numeric and fixed includes 'sigmaT', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a sigmaT to be estimated equal across rows with the same value.
#'                 'rho_0': If numeric and fixed includes 'rho_0', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_0 to be estimated equal across rows with the same value.
#'                 'rho_1': If numeric and fixed includes 'rho_1', values that will be used for computations for the row.
#'                              Otherwise, it will be treated as a label of a rho_1 to be estimated equal across rows with the same value.
#'                 'HR' : Numeric. Either a value between 0 and 1 given the Hit Rate or an integer giving the number of Hits for the row.
#'                 'FAR' : Numeric. Either a value between 0 and 1 given the False Alarm Rate or an integer giving the number of Misses for the row.
#'                 'RRTP' : Numeric. Either a value between 0 and 1 given the Rejection Rate on target present trials or an integer giving the number of Incorrect Rejections for the row.
#'                 'RRTA' : Numeric. Either a value between 0 and 1 given the Rejection Rate on target absent or an integer giving the number of Correct Rejections for the row.
#' @param fixed : Optional character vector. Name of parameter that will fixed at values in data data.frame. Can include any of
#'                 'criterion', 'mu', 'sigmaT', 'rho_0', 'rho_1'.
#' @param rates : Logical. Indicates whether columns HR, FAR, RRTP, and RRTA of data ought to be treated as proportions/rates (default) or as counts.
#' @param objective : Character. Whether to estimate by maximum likelihood, 'maxLik' (default) or by weighed least squares 'wls'
#' @param nInit : Numeric. Number of different initial parameters to try for the optimization routine. Defaults to 5. Should be at least 1
#' @param startingPars: Numeric vector. Starting values to use for (at least one of) the optimization runs.
#' @param control : List. control parameters to be passed to optim. Default is list(trace = 3)). See ?optim for details.
#' @param ... : additional parameters to be passed to the model specific functions.
#'
#' @importFrom stats optim qnorm runif
#' @importFrom magrittr "%>%"
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by inner_join mutate arrange desc group_map select
#' @importFrom tidyselect starts_with
#'
#' @export
#'
#' @return output
#'
#' @examples # Not run:
#'           # (exampleIO.ml    <- unforcedChoice_Estimation(data = toyIO,
#'           #                                               fixed = c("rho_0", "rho_1"),
#'           #                                               objective = "maxLik"))
#'           # (exampleIO.c2    <- unforcedChoice_Estimation(data = toyIO,
#'           #                                               fixed = c("rho_0", "rho_1"),
#'           #                                               objective = "wls"))
#'           # (exampleDO.ml    <- unforcedChoice_Estimation(data = toyDO,
#'           #                                               fixed = c("r1=r0"),
#'           #                                               objective = "maxLik"))
#'           # (exampleDO.c2    <- unforcedChoice_Estimation(data = toyDO,
#'           #                                               fixed = c("r1=r0"),
#'           #                                               objective = "wls"))
#'           # (exampleInt.ml   <- unforcedChoice_Estimation(data = toyInt,
#'           #                                               fixed = c("r1=r0"),
#'           #                                               objective = "maxLik"))
#'           # (exampleInt.ml.r <- unforcedChoice_Estimation(data = toyInt,
#'           #                                               fixed = c("r1=r0", "rho_0"),
#'           #                                               objective = "maxLik"))
#'           # (exampleInt.ml.s <- unforcedChoice_Estimation(data = toyInt,
#'           #                                               fixed = c("r1=r0", "sigmaT"),
#'           #                                               objective = "maxLik"))
#'           # (exampleInt.c2   <- unforcedChoice_Estimation(data = toyInt,
#'           #                                               fixed = c("r1=r0"),
#'           #                                               objective = "wls"))
#'           # (exampleInt.c2.r <- unforcedChoice_Estimation(data = toyInt,
#'           #                                               fixed = c("r1=r0", "rho_0"),
#'           #                                               objective = "wls"))
#'           # (exampleInt.c2.s <- unforcedChoice_Estimation(data = toyInt,
#'           #                                               fixed = c("r1=r0", "sigmaT"),
#'           #                                               objective = "wls"))
#'           # (exampleEns.ml   <- unforcedChoice_Estimation(data = toyEns,
#'           #                                               fixed = c("r1=r0"),
#'           #                                               method.Ens = "mvnorm",
#'           #                                               objective = "maxLik"))
#'           # (exampleEns.ml.r <- unforcedChoice_Estimation(data = toyEns,
#'           #                                               fixed = c("r1=r0", "rho_0"),
#'           #                                               method.Ens = "mvnorm",
#'           #                                               objective = "maxLik"))
#'           # (exampleEns.ml.s <- unforcedChoice_Estimation(data = toyEns,
#'           #                                               fixed = c("r1=r0", "sigmaT"),
#'           #                                               method.Ens = "mvnorm",
#'           #                                               objective = "maxLik"))
#'           # (exampleEns.c2   <- unforcedChoice_Estimation(data = toyEns,
#'           #                                               fixed = c("r1=r0"),
#'           #                                               method.Ens = "mvnorm", objective = "wls"))
#'           # (exampleEns.c2.r <- unforcedChoice_Estimation(data = toyEns,
#'           #                                               fixed = c("r1=r0", "rho_0"),
#'           #                                               method.Ens = "mvnorm", objective = "wls"))
#'           # (exampleEns.c2.s <- unforcedChoice_Estimation(data = toyEns,
#'           #                                               fixed = c("r1=r0", "sigmaT"),
#'           #                                               method.Ens = "mvnorm", objective = "wls"))
#'
unforcedChoice_Estimation <- function (data, fixed = NULL, rates = TRUE, objective = "maxLik",
                                       nInit = 5, startingPars = NULL, control = list(trace = 3), ...) {

    parameterNames <- c("criterion", "mu", "sigmaT", "rho_0", "rho_1")
    columns <- c("model", "nTP", "nTA", "nSize", parameterNames, "HR", "FAR", "RRTP", "RRTA")
    if (!all(columns %in% names(data))) {
        stop("data must be a data.frame with names '", paste(columns, collapse = ", "), "'")
    }

    if (any(data[, c("HR", "FAR", "RRTP", "RRTA")] < 0)) {
        stop("'HR', 'FAR', 'RRTP', and 'RRTA' must all be non-negative")
    }

    if (rates) {
        if (any(data[["HR"]] + data[["RRTP"]] > 1)) {
            stop("HR + RRTP cannot add up to more than 1")
        }
        if (any(data[["FAR"]] + data[["RRTA"]] > 1)) {
            stop("FAR + RRTA cannot add up to more than 1")
        }
    } else {
        if (any(data[["HR"]] + data[["RRTP"]] - data[['nTP']] < 0)) {
            stop("HR + RRTP cannot add up to more than nTP")
        }
        if (any(data[["FAR"]] + data[["RRTA"]] - data[['nTA']] < 0)) {
            stop("FAR + RRTA cannot add up to more than nTA")
        }
    }

    if (objective %in% c("maxLik", "wls")) {
        objFun <- switch(objective,
                         "maxLik" = unforcedChoice_lik,
                         "wls"    = unforcedChoice_wls
        )
    } else {
        stop("objective must be wither maxLik or wls")
    }


    initialPars <- lowerLimits <- upperLimits <- numeric()


    if (!is.null(fixed) && all(parameterNames %in% fixed)) {
        message("All parameters fixed. Nothing done.")
    } else {
        if (!("criterion" %in% fixed)) {

            if (any(by(data[['model']], data[['criterion']], length) > 1)) {
                warning("criterion should not be constrained to be equal for different models. Treating criterion values as different parameters for different models.")
                data[['criterion']] <- paste(data[['criterion']], data[['model']], sep = ".")
            }

            criterion <- qnorm((data[['RRTA']])^(1/data[['nSize']]))
            criterion[criterion == -Inf] <- -5
            criterion[criterion ==  Inf] <-  5
            criterion <- ifelse(data[['model']] == 'Ens' & criterion < 0, 0, criterion)
            lowerLimits <- c(lowerLimits, as.numeric(ifelse(data[['model']] == 'Ens', 0, -7)))

            criterion <- by(criterion, data[['criterion']], mean, simplify = TRUE)
            lowerLimits <- by(lowerLimits, data[['criterion']], mean, simplify = TRUE)

            names(criterion) <- paste("criterion", unique(data[['criterion']]), sep = ".")
            initialPars <- c(initialPars, criterion)

            upperLimits <- c(upperLimits, rep(7, length(criterion)))
        }

        if (!("mu" %in% fixed)) {
            means <-  qnorm((data[['RRTA']])^(1/data[['nSize']])) -
                qnorm(
                    data[['RRTP']] *
                        (data[['RRTA']])^(
                            data[['nSize']] /
                                (data[['nSize']] - 1)
                        )
                )
            means[is.nan(means)] <- 0

            means <- by(means, data[['mu']], mean, simplify = TRUE)
            names(means) <- paste("mu", unique(data[['mu']]), sep = ".")
            initialPars <- c(initialPars, means)

            lowerLimits <- c(lowerLimits, rep(0, length(means)))
            upperLimits <- c(upperLimits, rep(5, length(means)))
        }

        if (!("sigmaT" %in% fixed)) {
            sigmas <- rep(1, length(unique(data[['sigmaT']])))
            names(sigmas) <- paste("sigmaT", unique(data[['sigmaT']]), sep = ".")
            initialPars <- c(initialPars, sigmas)

            lowerLimits <- c(lowerLimits, rep(1, length(sigmas)))
            upperLimits <- c(upperLimits, rep(10, length(sigmas)))
        }

        if (!("rho_0" %in% fixed)) {
            rho0s <- ifelse(data[['model']] == 'IO', 0, .2)
            rho0s <- by(rho0s, data[['rho_0']], mean, simplify = TRUE)
            names(rho0s) <- paste("rho_0", unique(data[['rho_0']]), sep = ".")
            initialPars <- c(initialPars, rho0s)

            lowerLimits <- c(lowerLimits, rep(0, length(rho0s)))
            upperLimits <- c(upperLimits, rep(1, length(rho0s)))
        }
        if (("r1=r0" %in% fixed) || ("r1=sr0" %in% fixed)) {

            if ("rho_1" %in% fixed) {
                message("rho_1 is set to be a function of rho_0. Ignoring values in data")
            }

            if (("r1=r0" %in% fixed) & ("r1=sr0" %in% fixed)) {
                message("Only one 'r1=r0' or 'r1=sr0' can be included in fixed'. Choose one.")
            }

        } else {
            if (!("rho_1" %in% fixed)) {
                rho1s <- ifelse(data[['model']] == 'IO', 0, .15)
                rho1s <- by(rho1s, data[['rho_1']], mean, simplify = TRUE)
                names(rho1s) <- paste("rho_1", unique(data[['rho_1']]), sep = ".")
                initialPars <- c(initialPars, rho1s)

                lowerLimits <- c(lowerLimits, rep(0, length(rho1s)))
                upperLimits <- c(upperLimits, rep(1, length(rho1s)))
            }

        }


    }

    if (nInit > 1) {
        initialPars <- cbind(initialPars,
                             matrix(runif((nInit - 1) * length(upperLimits),
                                          min = lowerLimits,
                                          max = upperLimits
                             ), ncol = nInit - 1)
        )
    } else {
        initialPars <- matrix(runif(length(upperLimits),
                                    min = lowerLimits,
                                    max = upperLimits
        ), ncol = 1)
    }

    if (!is.null(startingPars)) {
        initialPars <- cbind(startingPars, initialPars)
    }

    optimized <- list()

    findOptim <- function (par, fn, data, fixed, rates, method,
                           control, lower, upper, ...) {
        optim(par = par,
              fn  = fn,
              data = data,
              fixed = fixed,
              rates = rates,
              method = method,
              control = control,
              lower = lower,
              upper = upper,
              ... = ...)
    }

    for (ii in seq(nInit)) {
        optimized[[ii]] <- findOptim(par = initialPars[, ii],
                                     fn = objFun,
                                     data = data,
                                     fixed = fixed,
                                     rates = rates,
                                     method = "L-BFGS-B",
                                     control = control,
                                     lower = lowerLimits,
                                     upper = upperLimits,
                                     ... = ...)
    }

    optimized <- optimized[[
        which.min(sapply(optimized, function (x) x[['value']]))
    ]]

    dataPars <- parsInData(pars = optimized$par, data = data, fixed = fixed)
    prediction <- predictData(dataPars, ... = ...)
    dataPred <- inner_join(x = dataPars,
                           y = prediction,
                           by = c("model", "nTP", "nTA", "nSize", "criterion",
                                  "mu", "sigmaT", "rho_0", "rho_1"),
                           suffix = c(".data", ".prediction"))

    dataLong <- pivot_longer(dataPred,
                             cols = starts_with(c("HR", "FAR", "RRTP", "RRTA")),
                             names_to = c(".value", "set"),
                             names_pattern = "(HR|FAR|RRTP|RRTA).(data|prediction)") %>%
        select(.data$model, .data$nTP, .data$nTA, .data$nSize, .data$criterion,
               .data$mu, .data$sigmaT, .data$rho_0, .data$rho_1,
               .data$set, .data$HR, .data$FAR, .data$RRTP, .data$RRTA)

    groupedData <- group_by(dataLong,
                            .data$model, .data$nSize, .data$mu, .data$sigmaT, .data$rho_0, .data$rho_1)

    ## TOCs
    groupedData %>% arrange(desc(.data$set)) %>%
        group_map(~ plot_ROC(type   = 'TOC',
                             model  = .y[['model']],
                             nSize  = .y[['nSize']],
                             mu     = .y[['mu']],
                             sigmaT = .y[['sigmaT']],
                             rho_0  = .y[['rho_0']],
                             rho_1  = .y[['rho_1']],
                             method.Ens = 'mvnorm'
        ) +
            geom_point(data = .x, aes(shape = .data$set, color = .data$set)) +
            scale_color_manual(values = c("#13294b", "#ff552e")) +
            theme(legend.title=element_blank())
        ) ->
        dataTOCs

    ## IOCs
    groupedData %>% arrange(desc(.data$set)) %>%
        group_map(~ plot_ROC(type   = 'IOC',
                             model  = .y[['model']],
                             nSize  = .y[['nSize']],
                             mu     = .y[['mu']],
                             sigmaT = .y[['sigmaT']],
                             rho_0  = .y[['rho_0']],
                             rho_1  = .y[['rho_1']],
                             method.Ens = 'mvnorm'
        ) +
            geom_point(data = mutate(.x, FAR = 1 - RRTA), aes(shape = .data$set, color = .data$set)) +
            scale_color_manual(values = c("#13294b", "#ff552e")) +
            theme(legend.title=element_blank())
        ) ->
        dataIOCs

    ## EOCs
    groupedData %>% arrange(desc(.data$set)) %>%
        group_map(~ plot_ROC(type   = 'EOC',
                             model  = .y[['model']],
                             nSize  = .y[['nSize']],
                             mu     = .y[['mu']],
                             sigmaT = .y[['sigmaT']],
                             rho_0  = .y[['rho_0']],
                             rho_1  = .y[['rho_1']],
                             method.Ens = 'mvnorm'
        ) +
            geom_point(data = mutate(.x, HR = 1 - RRTP, FAR = 1 - RRTA), aes(shape = .data$set, color = .data$set)) +
            scale_color_manual(values = c("#13294b", "#ff552e")) +
            theme(legend.title=element_blank())
        ) ->
        dataEOCs



    output <- list(estimates = optimized,
                   data = dataPred,
                   TOCs = dataTOCs,
                   IOCs = dataIOCs,
                   EOCs = dataEOCs
    )

    return(output)
}
