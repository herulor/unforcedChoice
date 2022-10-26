#' @name toyIO
#' @title Example data that perfectly follows the IO model
#' @description data.frame with rates of an IO model
#' @docType data
#' @usage data(toyIO)
#' @format data.frame
NULL

#' @name toyDO
#' @title Example data that perfectly follows the DO model
#' @description data.frame with rates of an DO model
#' @docType data
#' @usage data(toyDO)
#' @format data.frame
NULL

#' @name toyInt
#' @title Example data that perfectly follows the Integration model
#' @description data.frame with rates of an Integration model
#' @docType data
#' @usage data(toyInt)
#' @format data.frame
NULL

#' @name toyEns
#' @title Example data that perfectly follows the Ensemble model
#' @description data.frame with rates of an Ensemble model
#' @docType data
#' @usage data(toyEns)
#' @format data.frame
NULL

# toyIO  <- cbind(model = "IO", nTP = 100, nTA = 100,
#                nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0, rho_1 = 0,
#                Compute_Rates(model = "IO", criterion = c(-Inf, 0, .5, .8, 1, 1.2, 1.5, 1.7, 2, 2.5, Inf)))
# save(toyIO, file = "./data/toyIO.rda")
#
# toyDO  <- cbind(model = "DO", nTP = 100, nTA = 100,
#                nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0.2, rho_1 = 0.2,
#                Compute_Rates(model = "DO", criterion = c(-Inf, 0, .5, .8, 1, 1.2, 1.5, 1.7, 2, 2.5, Inf),
#                              nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0.2, rho_1 = 0.2))
# save(toyDO, file = "./data/toyDO.rda")
#
# toyInt  <- cbind(model = "Int", nTP = 100, nTA = 100,
#                nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0.2, rho_1 = 0.2,
#                Compute_Rates(model = "Int", criterion = c(-Inf, -3, -2, -1, 0, .5, 1, 2, 3, Inf),
#                              nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0.2, rho_1 = 0.2))
# save(toyInt, file = "./data/toyInt.rda")
#
# toyEns <- cbind(model = "Ens", nTP = 100, nTA = 100,
#                nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0.2, rho_1 = 0.2,
#                Compute_Rates(model = "Ens", criterion = c(0, .5, .8, 1, 1.2, 1.5, 1.7, 2, Inf),
#                              nSize = 4, mu = 1, sigmaT = 1, rho_0 = 0.2, rho_1 = 0.2, method.Ens = "mvnorm"))
# save(toyEns, file = "./data/toyEns.rda")
#
