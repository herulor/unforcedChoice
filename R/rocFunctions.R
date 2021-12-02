#' plotROC
#'
#' @param data G
#' @param ratesTable G
#' @param type G
#' @param nSize G
#' @param mu_0 G
#' @param sigma_0 G
#' @param mu_1 G
#' @param sigma_1 G
#' @param rho G
#' @param title G
#' @param subtitle G
#' @param color G
#' @param alpha G
#' @param nSteps G
#' @param method G
#' @param fillerFAR G
#' @param AUC G
#' @param add G
#' @param ... G G
#'
#' @importFrom rlang expr
#' @importFrom stats qnorm integrate
#' @import ggplot2
#'
#' @return G
#' @export
#'
#' @examples plotROC(AUC = FALSE)
plotROC <- function(data = NULL, ratesTable = NULL, type = "ensemble",
                    nSize = 4, mu_0 = 0, sigma_0 = 1, mu_1 = 1, sigma_1 = 1, rho = 0,
                    title = NULL, subtitle = NULL, color = NULL, alpha = .3,
                    nSteps = 41, method = NULL, fillerFAR = TRUE, AUC = TRUE, add = FALSE,
                    ...) {


    if (sigma_0 != sigma_1) {
        warning("Heteroscedasticity not (yet) implemented. Taking sigma_1 equal to sigma_0")
        sigma_1 <- sigma_0
    }

    if (is.null(ratesTable)) {
        ratesTable <- ComputeRates(type = type, nSize = nSize,
                                   mu_0 = mu_0, mu_1 = mu_1, sigma_0 = sigma_0, sigma_1 = sigma_1, rho = rho,
                                   nSteps = nSteps, method = method, fillerFAR = fillerFAR, ...)
    }

    if (AUC) {
        if (rho < 1) {
            algorithm <- "simpson"
        } else {
            algorithm <- "trapezoid"
        }
        AUC <- ComputeAUC(lower = 0, upper = 1,
                          alg = algorithm, method = method, fillerFAR = fillerFAR,
                          type = type, nSize = nSize,
                          mu_0 = mu_0, mu_1 = mu_1, sigma_0 = sigma_0, sigma_1 = sigma_1, rho = rho,
                          ...)
    }

    if (is.null(title)) {
        if (nSize > 1) {
            title <- "pROC Curve"
        } else if (nSize == 1) {
            title <- "ROC Curve"
        }
    }

    if (is.null(subtitle)) {
        if (type == "ensemble") {
            model <- "Ensemble"

            if (rho == 1) {
                ratesTable <- ratesTable[c(1,  nrow(ratesTable) %/% 2, nrow(ratesTable)), ]
            }

        } else if (type == "projection") {
            if (rho == 0) {
                model <- "Independent Observations"
            } else if (rho > 0) {
                model <- "Dependent Observations"
            }
        }

        n <- mu <- sigma <- NULL
        subtitle <- expr(paste(!!model, " Model. ",
                               n, " = ", !!nSize, ", ",
                               mu[0], " = ", !!mu_0, ", ",
                               mu[T], " = ", !!mu_1, ", ",
                               sigma, " = ", !!sigma_0, ", ",
                               rho,   " = ", !!rho
        ))
    }

    if (is.null(color)) {
        color <- "#034521"
    }

    ratesPoly <- rbind(ratesTable, ratesTable[1, ])
    ratesPoly[nrow(ratesPoly), "HR"] <- 0

    if (type == "ensemble") {
        if (rho == 1) {
            ratesTable <- ratesTable[c(1, nrow(ratesTable), nrow(ratesTable)) %/% 2, ]
        }
    }
    square   <- data.frame(FAR = c(0, 1, 1, 0), HR = c(0, 0, 1, 1), FAR2 = c(1, 1, 0, 0), HR2 = c(0, 1, 1, 0))
    diagonal <- data.frame(FAR = 0, HR = 0, FAR2 = 1, HR2 = 1)

    plotROC <- ggplot(ratesTable, aes(x = .data$FAR, y = .data$HR)) +
        geom_polygon(data = ratesPoly, fill = color, alpha = alpha, linetype = 0) +
        geom_line(colour = color) +
        xlim(0, 1.01) + ylim(0, 1.01) +
        xlab("False Alarm Rate") + ylab("Hit Rate") +
        labs(title = title, subtitle = subtitle) +
        coord_fixed(ratio = 1) +
        geom_segment(data = square, aes(xend = .data$FAR2, yend = .data$HR2),   color = "grey55") +
        geom_segment(data = diagonal, aes(xend = .data$FAR2, yend = .data$HR2), color = "grey45") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))

    if (AUC) {
        if (nSize > 1) {
            AUCText <- paste0("pAUC = ", round(AUC, 3))
        } else if (nSize == 1) {
            AUCText <- paste0("AUC = ", round(AUC, 3))
        }
        plotROC <- plotROC +
            annotate("text", x = .8, y = .1, label = AUCText, size = 3.5)
    }


    if (!is.null(data)) {
        plotROC <- plotROC + geom_point(data = data)
    }

    return(plotROC)
}



ComputeRates <- function(type = "ensemble", nSize = 4, mu_0 = 0, sigma_0 = 1, mu_1 = 1, sigma_1 = 1, rho = 0,
                         fillerFAR = TRUE,
                         nSteps = 41, breaks = NULL, method = NULL, ...) {

    if (nSteps < 3) {
        warning("nSteps should be at least 3")
        nSteps <- 3
    }
    if (type == "projection") {
        if (is.null(breaks)) {
            criteria <- qnorm(seq(from = 0, to = 1, length = nSteps), mean = mu_1 - mu_0, sd = sqrt(2 * nSize) * sigma_0)
        } else {
            criteria <- qnorm(breaks)
        }

        if (is.null(method)) {
            method <- "cuhre"
        }


        HR <- sapply(criteria, function (y) target_ID_TP_projectionCrit(criterion = y,
                                                                        mu_0 = mu_0, sigma_0 = sigma_0,
                                                                        mu_1 = mu_1, sigma_1 = sigma_1,
                                                                        rho = rho, nSize = nSize,
                                                                        method = method,
                                                                        ...))


        FAR <- sapply(criteria, function (y) filler_ID_TA_projectionCrit(criterion = y,
                                                                         mu_0 = mu_0, sigma_0 = sigma_0,
                                                                         rho = rho, nSize = nSize,
                                                                         method = method,
                                                                         ...))

        if (fillerFAR) {
            scaleK <- nSize
        } else {
            scaleK <- 1
        }

        rates <- data.frame(criteria = criteria, HR = HR, FAR = scaleK * FAR)

    } else if (type == "ensemble") {
        if (is.null(breaks)) {
            criteria <- qnorm(seq(from = .5, to = 1, length = nSteps))
        } else {
            criteria <- qnorm(breaks)
        }

        if (is.null(method)) {
            method <- "mvnorm"
        }

        HR <- sapply(criteria, function (y) target_ID_TP_ensembleCrit(criterion = y,
                                                                   mu_1 = mu_1, sigma_0 = sigma_1,
                                                                   rho = rho, nSize = nSize,
                                                                   method = method,
                                                                   ...))


        FAR <- sapply(criteria, function (y) filler_ID_TA_ensembleCrit(criterion = y,
                                                                    mu_0 = mu_0, sigma_0 = sigma_0,
                                                                    rho = rho, nSize = nSize,
                                                                    method = method,
                                                                    ...))

        rates <- data.frame(criteria = criteria, HR = HR, FAR = FAR)

    } else {
        stop("type must be either 'projection' or 'ensemble'.")
    }

    return(rates)

}





ComputeAUC <- function (type = "ensemble", nSize = 4, mu_0 = 0, sigma_0 = 1, mu_1 = 1, sigma_1 = 1, rho = 0,
                        lower = 0, upper = 1, alg = "simpson", fillerFAR = TRUE,
                        maxIter = 100, tol = 5e-4, method = NULL, ...) {

    if (type == "ensemble") {
        if (lower < .5) {
            lower <- .5
        }
        if (is.null(method)) {
        method <- "mvnorm"
        }
    } else if (type == "projection") {
        if (is.null(method)) {
        method <- "cuhre"
        }
    } else {
        stop("type must be either 'ensemble' or 'projection'")
    }

    trapezoid <- function(a, b, n = 100, type = "ensemble") {
        # numerical integral from a to b
        # using the trapezoid rule with n subdivisions
        # assume a < b and n is a positive integer
        # Modified from
        # https://www.r-bloggers.com/2010/12/one-dimensional-integrals/
        # Guangchuang Yu

        if (n < 2) {
            stop("n needs to be 2 or larger")
        }

        h <- (b - a) / n
        x <- seq(a, b, by = h)

        rates <- ComputeRates(type = type, nSize = nSize,
                              mu_0 = mu_0, sigma_0 = sigma_0, mu_1 = mu_1, sigma_1 = sigma_1, rho = rho,
                              fillerFAR = fillerFAR, breaks = x,
                              method = method)


        y  <- rates[['HR']]
        hi <- -diff(rates[['FAR']])

        s <- sum(.5 * hi * (y[-1] + y[-length(y)]))
        return(s)
    }

    simpson <- function(a, b, n = 100, type = "ensemble") {
        # numerical integral using Simpson's rule
        # assume a < b and n is an even positive integer
        # Modified from
        # https://www.r-bloggers.com/2010/12/one-dimensional-integrals/
        # Guangchuang Yu
        #
        # Unequal intervals Simpson's rule from
        # A.  K.  Singh  and  G.  R.  Thorpe,  Simpson’s  1/3-rule  of  integration  for  unequal  divisions of integration domain,
        # J. Concrete Applicable Maths., 1(3) (2003), 247-252.


        if ((n %% 2) != 0) {
            warning("n should be even. Using n + 1")
            n <- n + 1
        }

        h <- (b - a) / n
        x <- seq(a, b, by = h)

        rates <- ComputeRates(type = type, nSize = nSize,
                              mu_0 = mu_0, sigma_0 = sigma_0, mu_1 = mu_1, sigma_1 = sigma_1, rho = rho,
                              fillerFAR = fillerFAR, breaks = x,
                              method = method)


        y  <- rates[['HR']]
        hi <- -diff(rates[['FAR']])

        hiOdd  <- hi[seq_along(hi) %% 2 == 1]
        hiEven <- hi[seq_along(hi) %% 2 == 0]

        yOdd  <- y[seq_along(y) %% 2 == 1]
        yEven <- y[seq_along(y) %% 2 == 0]

        seqA <- ifelse(hiOdd == 0, 0, (hiOdd + hiEven) * (2 * hiOdd - hiEven) / hiOdd)
        seqB <- ifelse(hiOdd * hiEven == 0, 0, ((hiOdd + hiEven)^3) / (hiOdd * hiEven))
        seqC <- ifelse(hiEven == 0, 0, (hiOdd + hiEven) * (2 * hiEven - hiOdd) / hiEven)

        seqA <- seqA * yOdd[-length(yOdd)]
        seqB <- seqB * yEven
        seqC <- seqC * yOdd[-1]

        s <- sum(seqA + seqB + seqC) / 6

        return(s)
    }

    quadrature_internal <- function(S.old, a, m, b, tol, level, type = "ensemble", alg = "simpson", maxIter = maxIter) {
        # numerical integration using adaptive quadrature
        #
        # Modified from
        # https://www.r-bloggers.com/2010/12/one-dimensional-integrals/
        # Guangchuang Yu
        if (level > maxIter) {
            cat ("recursion limit reached: singularity likely\n")
            return (NULL)
        }

        if (alg == "simpson") {
            S.left  <- simpson(a = a, b = m, n = 2, type = type)
            S.right <- simpson(a = m, b = b, n = 2, type = type)
        } else if (alg == "trapezoid") {
            S.left  <- trapezoid(a = a, b = m, n = 2, type = type)
            S.right <- trapezoid(a = m, b = b, n = 2, type = type)
        }
        S.new <- S.left + S.right
        error <- abs(S.new - S.old)

        if (error > tol) {
            S.left  <- quadrature_internal(S.old = S.left, a = a, m = (a + m) / 2, b = m,
                                           tol = tol / 2, level = level + 1,
                                           type = type, alg = alg, maxIter = maxIter)
            S.right <- quadrature_internal(S.old = S.right, a = m, m = (m + b) / 2, b = b,
                                           tol = tol / 2, level = level + 1,
                                           type = type, alg = alg, maxIter = maxIter)
            S.new   <- S.left + S.right
        }

        attr(S.new, "tol") <- tol
        return(S.new)
    }

    level <- 1
    S.old <- .5
    S.new <- quadrature_internal(S.old = S.old, a = lower, m = (lower + upper) / 2, b = upper,
                                 tol = tol, level = level + 1, type = type, alg = alg, maxIter = maxIter)

    return(S.new)
}
