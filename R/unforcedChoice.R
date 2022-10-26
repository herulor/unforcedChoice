# Unforced choice models ----
# Víctor H Cervantes
# victorhc@illinois.edu
# vhcervantesb@unal.edu.co


## unforcedChoice_HR() ----
#'  Computes the probability of choosing the target stimulus when a target is shown (Hit Rate) for the
#'  model in the argument. Default is Independent Observations Model (model = "IO")
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param model : String indicating the model for which the Hit Rate will be computed. Default "IO"
#' @param ...   : Additional parameters to be passed to the specific model function
#'
#' @return numeric value of the Hit Rate
#' @export
#' @importFrom mvtnorm pmvnorm
#'
#' @examples unforcedChoice_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice_HR(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
unforcedChoice_HR  <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                                nSize = 4, model = "IO", ...) {
    funHR <- switch(model,
                    "IO"  = IO_HR,
                    "DO"  = DO_HR,
                    "Int" = Int_HR,
                    "Ens" = Ens_HR,
                    "Custom" = Custom_HR)

    HR <- funHR(criterion = criterion, mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                nSize = nSize, ... = ...)
    return(HR)
}


## unforcedChoice_FAR() ----
#'  Computes the probability of choosing the target stimulus when no target is shown (False Alarm Rate) for the
#'  model in the argument. Default is Independent Observations Model (model = "IO")
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param model : String indicating the model for which the Hit Rate will be computed. Default "IO"
#' @param ...   : Additional parameters to be passed to the specific model function
#'
#' @return numeric value of the False Alarm Rate
#' @export
#'
#' @examples unforcedChoice_FAR(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice_FAR(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice_FAR(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
unforcedChoice_FAR <- function (criterion = .25, mu = 0, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                                nSize = 4, model = "IO", ...) {
    funHR <- switch(model,
                    "IO"  = IO_FAR,
                    "DO"  = DO_FAR,
                    "Int" = Int_FAR,
                    "Ens" = Ens_FAR,
                    "Custom" = Custom_FAR)

    FAR <- funHR(criterion = criterion, mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                 nSize = nSize, ... = ...)
    return(FAR)
}


## unforcedChoice_RRTP() ----
#'  Computes the probability of rejecting the set when no target is shown (Rejection Rate) for the
#'  model in the argument. Default is Independent Observations Model (model = "IO")
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param model : String indicating the model for which the Hit Rate will be computed. Default "IO"
#' @param ...   : Additional parameters to be passed to the specific model function
#'
#' @return numeric value of the Rejection Rate
#' @export
#'
#' @examples unforcedChoice_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice_RRTP(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
unforcedChoice_RRTP <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                                 nSize = 4, model = "IO", ...) {
    funHR <- switch(model,
                    "IO"  = IO_RRTP,
                    "DO"  = DO_RRTP,
                    "Int" = Int_RRTP,
                    "Ens" = Ens_RRTP,
                    "Custom" = Custom_RRTP)

    RR <- funHR(criterion = criterion, mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                nSize = nSize, ... = ...)
    return(RR)
}


## unforcedChoice_RRTA() ----
#'  Computes the probability of rejecting the set when no target is shown (Rejection Rate) for the
#'  model in the argument. Default is Independent Observations Model (model = "IO")
#'
#' @param criterion : Value of the criterion for selection
#' @param mu : Mean of the target marginal distributions
#' @param sigmaT : Standard deviation of the target marginal distributions
#' @param rho_0 : Value of the common correlation among non-target stimuli
#' @param rho_1 : Value of the common correlation among non-target stimuli and the target
#' @param nSize : Number of simultaneously presented stimuli
#' @param model : String indicating the model for which the Hit Rate will be computed. Default "IO"
#' @param ...   : Additional parameters to be passed to the specific model function
#'
#' @return numeric value of the Rejection Rate
#' @export
#'
#' @examples unforcedChoice_RRTA(mu = 1, sigmaT = 1, nSize = 2, criterion = .25)
#'           unforcedChoice_RRTA(mu = 1, sigmaT = 1, nSize = 2, criterion = .50)
#'           unforcedChoice_RRTA(mu = 1, sigmaT = 1, nSize = 2, criterion = -Inf)
#'
unforcedChoice_RRTA <- function (criterion = .25, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
                                 nSize = 4, model = "IO", ...) {
    funHR <- switch(model,
                    "IO"  = IO_RRTA,
                    "DO"  = DO_RRTA,
                    "Int" = Int_RRTA,
                    "Ens" = Ens_RRTA,
                    "Custom" = Custom_RRTA)

    RR <- funHR(criterion = criterion, mu = mu, sigmaT = sigmaT, rho_0 = rho_0, rho_1 = rho_1,
                nSize = nSize, ... = ...)
    return(RR)
}
