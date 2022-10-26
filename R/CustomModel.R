# Functions for the Custom Model ----
# Víctor H Cervantes
# victorhc@illinois.edu
# vhcervantesb@unal.edu.co


## Custom_RR() ----
#'  Computes the probability of rejecting the set when (Rejection Rate) for
#'  Custom Models
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli. Must equal ncol(cMatrix).
#' @param cMatrix : Matrix that defines the custom model
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @importFrom mvtnorm pmvnorm
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Custom_RR(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = .25, cMatrix = diag(2))
#'           unforcedChoice:::Custom_RR(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = .50, cMatrix = diag(2))
#'           unforcedChoice:::Custom_RR(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = -Inf, cMatrix = diag(2))
#'
Custom_RR <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                       nSize = 4, cMatrix, ...) {

    # Parameter verification
    if (criterion == 0) {
        return(0)
    }

    if (rho_0 == 1) {
        if (abs(rho_1) == 1) {
            rho_1 <- sign(rho_1) * .99999
        }
        rho_0 <- .99999
    }

    if (nSize != ncol(cMatrix)) {
        stop("nSize must equal ncol(cMatrix)")
    }

    if (length(criterion) == 1 && nrow(cMatrix) > 1) {
        message("length(criterion) equals 1 and nrow(cMatrix) is greater than 1.\n rep(criterion, nrow(cMatrix)) will be used as criterion.")
        criterion <- rep(criterion, nrow(cMatrix))
    } else if (length(criterion) != nrow(cMatrix)) {
        stop("length(criterion) should equal 1 or nrow(cMatrix)")
    }

    vectorMu <- as.numeric(cMatrix %*% c(mu, rep(0, nSize - 1)))

    covMatrix <- cbind(c(
        sigmaT^2,
        rep(
            sigmaT * rho_1,
            (nSize - 1))
    ),
    rbind(
        rep(
            sigmaT * rho_1,
            (nSize - 1)
        ),
        (
            (1 - rho_0) * diag(nrow = nSize - 1)
        ) + (
            rho_0
        )
    )
    )

    covMatrix <- cMatrix %*% covMatrix %*% t(cMatrix)

    probability <- pmvnorm(lower = rep(-Inf, nSize),
                                    upper = criterion,
                                    mean = vectorMu,
                                    sigma = covMatrix)

    return(probability)
}


## Custom_RRTA() ----
#'  Computes the probability of rejecting the set when no target is shown (Rejection Rate) for
#'  Custom Models
#'
#' @param criterion : Value of the criterion for selection
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param nSize : Number of simultaneously presented stimuli. Must equal ncol(cMatrix).
#' @param cMatrix : Matrix that defines the custom model
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Custom_RRTA(nSize = 2, criterion = .25,
#'                                        cMatrix = diag(2))
#'           unforcedChoice:::Custom_RRTA(nSize = 2, criterion = .50,
#'                                        cMatrix = diag(2))
#'           unforcedChoice:::Custom_RRTA(nSize = 2, criterion = -Inf,
#'                                        cMatrix = diag(2))
#'
Custom_RRTA <- function (criterion = .25, rho_0 = 0,
                         nSize = 4, cMatrix, ...) {

    if (nSize == 1) {
        stop("Customemble model does not apply to set size 1.")
    } else if (nSize > 1) {
        probability <- Custom_RR(criterion = criterion, mu = 0, sigmaT = 1, rho_0 = rho_0, rho_1 = rho_0,
                                 nSize = nSize, cMatrix = cMatrix)
    }

    return(probability)
}


## Custom_RRTP() ----
#'  Computes the probability of rejecting the set when a target is shown (Rejection Rate) for
#'  Custom Models
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli. Must equal ncol(cMatrix).
#' @param cMatrix : Matrix that defines the custom model
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTP()
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Custom_RRTP(mu = 1, sigmaT = 1, nSize = 2,
#'                                        criterion = .25, cMatrix = diag(2))
#'           unforcedChoice:::Custom_RRTP(mu = 1, sigmaT = 1, nSize = 2,
#'                                        criterion = .50, cMatrix = diag(2))
#'           unforcedChoice:::Custom_RRTP(mu = 1, sigmaT = 1, nSize = 2,
#'                                        criterion = -Inf, cMatrix = diag(2))
#'
Custom_RRTP <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                         nSize = 4, cMatrix, ...) {

    if (nSize == 1) {
        stop("Customemble model does not apply to set size 1.")
    } else if (nSize > 1) {
        probability <- Custom_RR(criterion = criterion, mu = mu, sigmaT = sigmaT,
                                 rho_0 = rho_0, rho_1 = rho_1,
                                 nSize = nSize, cMatrix = cMatrix)
    }

    return(probability)
}


## Custom_ID() ----
#'  Computes the probability of selecting the arbitrary first stimulus for
#'  Custom Models
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli. Must equal ncol(cMatrix).
#' @param cMatrix : Matrix that defines the custom model
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_RRTA()
#'
#' @importFrom mvtnorm pmvnorm
#'
#' @return numeric value of the Hit Rate
#' @export
#'
#' @examples unforcedChoice:::Custom_ID(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = .25, cMatrix = matrix(c(1, 0), nrow = 1))
#'           unforcedChoice:::Custom_ID(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = .50, cMatrix = matrix(c(1, 0), nrow = 1))
#'           unforcedChoice:::Custom_ID(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = -Inf, cMatrix = matrix(c(1, 0), nrow = 1))
#'
Custom_ID <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                       nSize = 4, cMatrix, ...) {
    # Parameter verification
    if (criterion == 0) {
        return(0)
    }

    if (rho_0 == 1) {
        if (abs(rho_1) == 1) {
            rho_1 <- sign(rho_1) * .99999
        }
        rho_0 <- .99999
    }

    if (nSize != ncol(cMatrix)) {
        stop("nSize must equal ncol(cMatrix)")
    }

    if (length(criterion) == 1 && nrow(cMatrix) > 1) {
        message("length(criterion) equals 1 and nrow(cMatrix) is greater than 1.\n rep(criterion, nrow(cMatrix)) will be used as criterion.")
        criterion <- rep(criterion, nrow(cMatrix))
    } else if (length(criterion) != nrow(cMatrix)) {
        stop("length(criterion) should equal 1 or nrow(cMatrix)")
    }

    cMatrix <- rbind(cMatrix,
                     cbind(1, -diag(nSize - 1)))

    vectorMu <- as.numeric(cMatrix %*% c(mu, rep(0, nSize - 1)))

    covMatrix <- cbind(c(
        sigmaT^2,
        rep(
            sigmaT * rho_1,
            (nSize - 1))
    ),
    rbind(
        rep(
            sigmaT * rho_1,
            (nSize - 1)
        ),
        (
            (1 - rho_0) * diag(nrow = nSize - 1)
        ) + (
            rho_0
        )
    )
    )

    covMatrix <- cMatrix %*% covMatrix %*% t(cMatrix)

    probability <- pmvnorm(lower = c(criterion, rep(0, nSize - 1)),
                                    upper = rep(Inf, nSize),
                                    mean = vectorMu,
                                    sigma = covMatrix)


    return(probability)
}


## Custom_HR() ----
#'  Computes the probability of choosing the target stimulus when a target is shown (Hit Rate) for
#'  Custom Models
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli. Must equal ncol(cMatrix).
#' @param cMatrix : Matrix that defines the custom model
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_HR()
#'
#' @return numeric value of the Hit Rate
#'
#' @examples unforcedChoice:::Custom_HR(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = .25, cMatrix = matrix(c(1, 0), nrow = 1))
#'           unforcedChoice:::Custom_HR(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = .50, cMatrix = matrix(c(1, 0), nrow = 1))
#'           unforcedChoice:::Custom_HR(mu = 1, sigmaT = 1, nSize = 2,
#'                                      criterion = -Inf, cMatrix = matrix(c(1, 0), nrow = 1))
#'
Custom_HR <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                       nSize = 4, cMatrix, ...) {

    if (nSize == 1) {
        stop("Customemble model does not apply to set size 1.")
    } else if (nSize > 1) {
        probability <- Custom_ID(criterion = criterion, mu = mu, sigmaT = sigmaT,
                                 rho_0 = rho_0, rho_1 = rho_1,
                                 nSize = nSize, cMatrix = cMatrix)
    }

    return(probability)
}


## Custom_FAR() ----
#'  Computes the probability of choosing the one non-target stimulus when no target is shown (False Alarm Rate) for
#'  Custom Models
#'
#' @param criterion : Value of the criterion for selection
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param nSize : Number of simultaneously presented stimuli. Must equal ncol(cMatrix).
#' @param cMatrix : Matrix that defines the custom model
#' @param ... : Captures unused parameter when called from generic function unforcedChoice_FAR()
#'
#' @return numeric value of the False Alarm Rate
#' @export
#'
#' @examples unforcedChoice:::Custom_FAR(nSize = 2, criterion = .25,
#'                                       cMatrix = matrix(c(1, 0), nrow = 1))
#'           unforcedChoice:::Custom_FAR(nSize = 2, criterion = .50,
#'                                       cMatrix = matrix(c(1, 0), nrow = 1))
#'           unforcedChoice:::Custom_FAR(nSize = 2, criterion = -Inf,
#'                                       cMatrix = matrix(c(1, 0), nrow = 1))
#
Custom_FAR <- function (criterion = .25, rho_0 = 0,
                        nSize = 4, cMatrix, ...) {

    if (nSize == 1) {
        stop("Customemble model does not apply to set size 1.")
    } else if (nSize > 1) {
        probability <- Custom_ID(criterion = criterion, mu = 0, sigmaT = 1, rho_0 = rho_0, rho_1 = rho_0,
                                 nSize = nSize, cMatrix = cMatrix)
    }

    return(probability)
}
