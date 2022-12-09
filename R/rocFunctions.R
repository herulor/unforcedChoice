# Functions to plot the ROCs of the models ----
# TOC: Target-Operating Characteristic curve
# IOC: Identification-Operating Characteristic curve
# EOC: Ensemble-Operating Characteristic curve

# Víctor H Cervantes
# victorhc@illinois.edu
# vhcervantesb@unal.edu.co


## Compute_Rates() ----
#' Computes the HR, FAR, RR_TA, and RR_TP of the model for a number of criterion values.
#'
#' @param model : String indicating the model for which the Hit Rate will be computed. Default "IO"
#' @param nSize : Number of simultaneously presented stimuli
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param criterion : Optional. Values of the criterion at which the rates will be computed.
#'                   If missing, breaks or nSteps will be used to determine the criteria.
#' @param breaks : Optional. Numeric values between 0 and 1. criterion values are created as qnorm(breaks)
#' @param nSteps : Used if neither criterion nor breaks is provided. criterion values are created as a vector
#'                 of length nSteps
#' @param method.Ens : Optional. Used to specify the method.Ens for the Ensemble model.
#' @param ...   : Additional parameters to be passed to the specific model function
#'
#' @importFrom stats qnorm
#'
#' @return data.frame with columns criterion, HR, FAR, RRTP, and RRTA
#' @export
#'
#' @examples Compute_Rates(model = "IO",  mu = 1, sigmaT = 1, nSize = 2)
#'           Compute_Rates(model = "DO",  mu = 1, sigmaT = 1, nSize = 2)
#'           Compute_Rates(model = "Int", mu = 1, sigmaT = 1, nSize = 2)
#'           Compute_Rates(model = "Ens", mu = 1, sigmaT = 1, nSize = 2)
#'
Compute_Rates <- function(model = "IO", nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                          criterion = NULL, breaks = NULL, nSteps = 41, method.Ens = NULL, ...) {

    if (nSteps < 3) {
        warning("nSteps should be at least 3. Computing for 3 steps.")
        nSteps <- 3
    }


    if (is.null(criterion)) {
        if (is.null(breaks)) {
            if (model %in% c("IO", "DO", "Int")) {
                criterion <- qnorm(seq(from = 0, to = 1, length = nSteps), mean = mu, sd = sqrt(2 * nSize))
            } else if (model == "Ens") {
                if (is.null(method.Ens)) {
                    method.Ens <- "mvnorm"
                } else {
                }
                criterion <- qnorm(seq(from = .5, to = 1, length = nSteps))
            } else if (model == "Custom") {
                criterion <- qnorm(seq(from = .5, to = 1, length = nSteps))
                warning(paste("Curstom model requested. It is preferable to also provide a suitable set of criterion breaks.",
                              "\nUsing", criterion))
            }
        } else {
            criterion <- qnorm(breaks)
        }
    }

    HR <- sapply(criterion, function (y) unforcedChoice_HR(criterion = y,
                                                          mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                                                          model = model, nSize = nSize, method.Ens = method.Ens,
                                                          ...))


    FAR <- sapply(criterion, function (y) unforcedChoice_FAR(criterion = y,
                                                            mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                                                            model = model, nSize = nSize, method.Ens = method.Ens,
                                                            ...))

    RRTP <- sapply(criterion, function (y) unforcedChoice_RRTP(criterion = y,
                                                              mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                                                              model = model, nSize = nSize, method.Ens = method.Ens,
                                                              ...))

    RRTA <- sapply(criterion, function (y) unforcedChoice_RRTA(criterion = y,
                                                              mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                                                              model = model, nSize = nSize, method.Ens = method.Ens,
                                                              ...))

    rates <- data.frame(criterion = criterion, HR = HR, FAR = FAR, RRTP = RRTP, RRTA = RRTA)

    return(rates)

}



## Compute_AUC() ----
#' Computes the area under the curve of the specified type ROC
#'
#' @param type : String indicating the type of ROC for which the AUC will be computed. Default "TOC"
#' @param model : String indicating the model for which the Hit Rate will be computed. Default "IO"
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param method.Ens : Optional. Used to specify the method.Ens for the Ensemble model.
#' @param lower : Numeric. Lower value of pnorm(criterion) for which the area will be computed. Defaults to 0
#' @param upper : Numeric. Upper value of pnorm(criterion) for which the area will be computed. Defaults to 1
#' @param alg : String indicating the algorithm to approximate the AUC. Default "simpson"
#' @param maxIter : Numeric indicating maximum number of iteration of the algorithm at each recursion level.
#' @param tol : Numeric indicating maximum tolerated difference between iterations of the algorithm
#' @param ...   : Additional parameters to be passed to the specific model function
#'
#' @return numeric value of the AUC with tol attribute
#' @export
#'
#' @examples Compute_AUC(mu = 1, sigmaT = 1, nSize = 2)
#'           Compute_AUC(mu = 1, sigmaT = 1, nSize = 2)
#'
Compute_AUC <- function (type = "TOC", model = "IO",
                         mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                         nSize = 4, method.Ens = NULL,
                         lower = 0, upper = 1, alg = "simpson",
                         maxIter = 100, tol = .001, ...) {

    if (model == "Ens") {
        if (lower < .5) {
            lower <- .5
        }
        if (is.null(method.Ens)) {
            method.Ens <- "mvnorm"
        }
        if (rho_0 >= .95 || abs(rho_1) >= .95) {
            alg <- "trapezoid"
        }
    }


    trapezoid <- function(a, b, n = 100, type = "TOC", model = "IO",
                          nSize, mu, sigmaT, rho_0, rho_1, method.Ens, ...) {
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

        rates <- Compute_Rates(model = model, nSize = nSize,
                               mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                               breaks = x,
                               method.Ens = method.Ens, ...)

        if (type != "EOC") {
            y  <- rates[['HR']]
        } else {
            y <- 1 - rates[['RRTP']]
        }

        if (type == "TOC") {
            hi <- -diff(rates[['FAR']])
        } else {
            hi <- -diff(1 - rates[['RRTA']])
        }

        s <- sum(.5 * hi * (y[-1] + y[-length(y)]))
        return(s)
    }

    simpson <- function(a, b, n = 100, type = "TOC", model = "IO",
                        nSize, mu, sigmaT, rho_0, rho_1, method.Ens, ...) {
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

        rates <- Compute_Rates(model = model, nSize = nSize,
                               mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                               breaks = x,
                               method.Ens = method.Ens, ...)

        if (type != "EOC") {
            y  <- rates[['HR']]
        } else {
            y <- 1 - rates[['RRTP']]
        }

        if (type == "TOC") {
            hi <- -diff(rates[['FAR']])
        } else {
            hi <- -diff(1 - rates[['RRTA']])
        }

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

    quadrature_internal <- function(S.old, a, m, b, tol, level, model = "IO", type = "TOC",
                                    nSize = nSize, method.Ens = NULL,
                                    mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                                    alg = "simpson", maxIter = maxIter, ...) {
        # numerical integration using adaptive quadrature
        #
        # Modified from
        # https://www.r-bloggers.com/2010/12/one-dimensional-integrals/
        # Guangchuang Yu
        if (level > maxIter) {
            cat("recursion limit reached: singularity likely\nReducing tol or increasig maxIter might reach a solution.")
            if (model == "Int" && alg == "trapezoid") {
                cat("For model == 'Int', try running alg = 'simpson'")
            }
            return(NULL)
        }

        if (alg == "simpson") {
            S.left  <- simpson(a = a, b = m, n = 2, model = model, type = type,
                               nSize = nSize, method.Ens = method.Ens,
                               mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1, ...)
            S.right <- simpson(a = m, b = b, n = 2, model = model, type = type,
                               nSize = nSize, method.Ens = method.Ens,
                               mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1, ...)
        } else if (alg == "trapezoid") {
            S.left  <- trapezoid(a = a, b = m, n = 2, model = model, type = type,
                                 nSize = nSize, method.Ens = method.Ens,
                                 mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1, ...)
            S.right <- trapezoid(a = m, b = b, n = 2, model = model, type = type,
                                 nSize = nSize, method.Ens = method.Ens,
                                 mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1, ...)
        }
        S.new <- S.left + S.right
        error <- abs(S.new - S.old)

        if (error > tol) {
            S.left  <- quadrature_internal(S.old = S.left, a = a, m = (a + m) / 2, b = m,
                                           tol = tol / 2, level = level + 1, model = model,
                                           type = type, alg = alg, maxIter = maxIter,
                                           nSize = nSize, method.Ens = method.Ens,
                                           mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1, ...)
            S.right <- quadrature_internal(S.old = S.right, a = m, m = (m + b) / 2, b = b,
                                           tol = tol / 2, level = level + 1, model = model,
                                           type = type, alg = alg, maxIter = maxIter,
                                           nSize = nSize, method.Ens = method.Ens,
                                           mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1, ...)
            S.new   <- S.left + S.right
        }

        attr(S.new, "tol") <- tol
        return(S.new)
    }

    level <- 1
    S.old <- 1 / nSize
    S.new <- quadrature_internal(S.old = S.old, a = lower, m = (lower + upper) / 2, b = upper,
                                 tol = tol, level = level + 1, type = type, alg = alg, maxIter = maxIter,
                                 nSize = nSize, model = model, method.Ens = method.Ens,
                                 mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1, ...)

    return(S.new)
}


## plot_ROC() ----
#' Plots the ROC of the indicated type (TOC, IOC, EOC) for the unforced choice model.
#'
#' @param model : String indicating the model for which the Hit Rate will be computed. Default "IO"
#' @param nSize : Number of simultaneously presented stimuli
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param criterion : Optional. Values of the criterion at which the rates will be computed.
#'                   If missing, breaks or nSteps will be used to determine the criteria.
#' @param breaks : Optional. Numeric values between 0 and 1. criterion  values are created as qnorm(breaks)
#' @param nSteps : Used if neither criterion nor breaks is provided. criterion values are created as a vector
#'                 of length nSteps
#' @param method.Ens : Optional. Used to specify the method.Ens for the Ensemble model.
#' @param ...   : Additional parameters to be passed to the specific model function
#' @param data : Optional data.frame with columns FAR, HR, RRTA, and RRTP of empirical points to add to the plot.
#' @param ratesTable : Optional data.frame with the structure of the output of Compute_Rates. If provided the curve and shading in the figure will used this table.
#' @param type : String indicating the type of ROC for which the AUC will be computed. Default "TOC"
#' @param AUC : Logical. Indicates whether to print or not AUCtext. Default is FALSE.
#' @param AUCtext : Optional character. If not ptovided, it will be paste0('AUC = ', Compute_AUC(.)) with the call to Compute_AUC using the parameters in the plot_ROC call.
#' @param title : Optional character. If not provided will be a text indicating the type of ROC plotted
#' @param subtitle : Optional character. If not provided it will show the model and parameter values.
#' @param color : Optional character to set shading color.
#' @param colorData : Optional character to set data points color.
#' @param alpha : Optional numeric value for alpha transparency.
#' @param xlabText : Optional character. If not provided it will be 'FAR + FIR[TA]' or 'False Alarm Rate' depending on type of ROC
#' @param ylabText : Optional character. If not provided it will be 'HR + FIR[TP]' or 'Hit Rate' depending on type of ROC
#' @param linetype : Optional character with the linetype for the ROC curve. See ?ggplot2::geom_line
#' @param theme : Logical. If TRUE (Default) theme_minimal will be used for the plot. See ?ggplot2::theme_minimal
#'
#' @return Object of class ggplot. Plot of the ROC.
#' @importFrom rlang expr
#' @import ggplot2
#' @export
#'
#' @examples plot_ROC(mu = 1, sigmaT = 1, nSize = 2, AUC = TRUE, type = "TOC")
#'
plot_ROC <- function(data = NULL, ratesTable = NULL, type = "TOC",
                     model = "IO", nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                     criterion = NULL, nSteps = 41, breaks = NULL, method.Ens = NULL,
                     AUC = FALSE, AUCtext = "",
                     title = NULL, subtitle = NULL, color = NULL, colorData = NULL, alpha = .3,
                     xlabText = NULL, ylabText = NULL,
                     linetype = "solid", theme = TRUE,
                     ...) {


    if (is.null(ratesTable)) {
        ratesTable <- Compute_Rates(model = model, nSize = nSize,
                                    mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                                    criterion = criterion, breaks = breaks, nSteps = nSteps,
                                    method.Ens = method.Ens, ...)
    }

    if (AUC & AUCtext == "") {
        if (rho_0 > .95 || abs(rho_1) > .95) {
            algorithm <- "trapezoid"
        } else {
            algorithm <- "simpson"
        }
        cat("Computing AUC... This may take a while")
        AUCtext <- Compute_AUC(lower = 0, upper = 1,
                               alg = algorithm, model = model, method.Ens = method.Ens,
                               mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                               type = type, nSize = nSize,
                               ...)
        AUCtext <- paste0("AUC = ", round(AUCtext, 3))
    }


    FIR <- TA <- TP <- NULL
    if (type == "TOC") {
        if (is.null(xlabText)) {
            xlabText <- "False Alarm Rate"
        }
    } else {
        if (is.null(xlabText)) {
            xlabText <-  expr(paste("FAR + ", FIR[TA]))
        }
        ratesTable[['FAR']] <- 1 - ratesTable[['RRTA']]
    }

    if (type != "EOC") {
        if (is.null(ylabText)) {
            ylabText <- "Hit Rate"
        }
    } else {
        if (is.null(ylabText)) {
            ylabText <-  expr(paste("HR + ", FIR[TP]))
        }
        ratesTable[['HR']] <- 1 - ratesTable[['RRTP']]
    }

    if (is.null(title)) {
        if (type == "TOC") {
            title <- "Target-Operating Characteristic Curve"
        } else if (type == "IOC") {
            title <- "Identification-Operating Characteristic Curve"
        } else if (type == "EOC") {
            title <- "Ensemble-Operating Characteristic Curve"
        }
    }

    if (is.null(subtitle)) {
        if (model == "Ens") {
            model <- "Ensemble"

            if (rho_0 == 1 || abs(rho_1) == 1) {
                ratesTable <- ratesTable[c(1,  nrow(ratesTable) %/% 2, nrow(ratesTable)), ]
            }

        } else if (model == "IO") {
            model <- "Independent Observations"
        } else if (model == "DO") {
            if (rho_0 == 0 & rho_1 == 0) {
                model <- "Independent Observations"
            } else {
                model <- "Dependent Observations"
            }

        } else if (model == "Int") {
            model <- "Integration"
        }

        n <- rho <- sigma <- NULL
        subtitle <- expr(paste(!!model, " Model. ",
                               n, " = ", !!nSize, ", ",
                               mu, " = ", !!round(mu, 3), ", ",
                               sigma[T], " = ", !!round(sigmaT, 3), ", ",
                               rho[0],   " = ", !!round(rho_0, 3), ", ",
                               rho[1],   " = ", !!round(rho_1, 3)
        ))
    }

    if (is.null(color)) {
        color <- "#FF552E"
    }

     if (is.null(colorData)) {
        colorData <- "#13294B"
    }

    ratesPoly <- rbind(ratesTable, ratesTable[1, ])
    ratesPoly[nrow(ratesPoly), "HR"] <- 0

    if (model == "Ens") {
        if (rho_0 == 1 || abs(rho_1) == 1) {
            ratesTable <- ratesTable[c(1, nrow(ratesTable), nrow(ratesTable)) %/% 2, ]
        }
    }
    square   <- data.frame(FAR = c(0, 1, 1, 0), HR = c(0, 0, 1, 1), FAR2 = c(1, 1, 0, 0), HR2 = c(0, 1, 1, 0))
    diagonal <- data.frame(FAR = 0, HR = 0, FAR2 = 1, HR2 = 1)

    plotROC <- ggplot(ratesTable, aes(x = .data$FAR, y = .data$HR)) +
        geom_polygon(data = ratesPoly, fill = color, alpha = alpha, linetype = 0) +
        geom_line(colour = color, linetype = linetype) +
        xlim(0, 1.01) + ylim(0, 1.01) +
        xlab(xlabText) + ylab(ylabText) +
        labs(title = title, subtitle = subtitle) +
        coord_fixed(ratio = 1) +
        geom_segment(data = square, aes(xend = .data$FAR2, yend = .data$HR2),   color = "grey55") +
        geom_segment(data = diagonal, aes(xend = .data$FAR2, yend = .data$HR2), color = "grey45") +
        theme(plot.title = element_text(hjust = 0.5))

    if (theme) {
        plotROC <- plotROC +
            theme_minimal()
    }

    if (AUC) {
        plotROC <- plotROC +
            annotate("text", x = .8, y = .1, label = AUCtext, size = 3.5)
    }


    if (!is.null(data)) {
        plotROC <- plotROC + geom_point(data = data, color = colorData, size = 1.5)
    }

    return(plotROC)
}
