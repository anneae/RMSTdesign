#' Create a Survdef Object for a Weibull Distribution
#'
#' @param shape the shape parameter for the Weibull distribution.
#' Parametrization is according to built-in R functions for the Weibull
#' distribution, see `?Weibull` for more information.
#' @param scale the scale parameter for the Weibull distribution.
#'
#' @return a list with components:
#' \item{S}{a vectorized function that takes time as input and returns the survival probability at that time}
#' \item{h}{a vectorized function that takes time as input and returns the hazard at that time}
#' @export
#'
#' @importFrom stats pweibull
#' @importFrom stats dweibull
#'
#' @examples survdefWeibull(shape = 1.05, scale = 8573)
survdefWeibull<-function(shape, scale){
    return(list(S = function(t) pweibull(t,shape = shape, scale = scale, lower.tail = F),
                h = function(t) dweibull(t, shape = shape, scale = scale)/pweibull(t, shape = shape, scale = scale, lower.tail = F)))
}
