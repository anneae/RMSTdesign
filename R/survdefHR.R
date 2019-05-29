#' Create a Survdef Object Based on a Hazard Ratio
#'
#' Creates a new object which stores user-specified survival distribution
#' information in the format needed for the main function, `RMSTpow`.
#' `survdefHR` is used when the user wishes to specify a survival distribution
#' that is defined by its relationship to another distribution via a constant
#' hazard ratio.
#'
#' @param survdefC the survival distribution for the reference/control group,
#' as a list in the form output by `survdef`.
#' @param HR the hazard ratio defining the relationship between the two distributions.
#'
#' @return a list with components:
#' \item{S}{a vectorized function that takes time as input and returns the survival probability at that time}
#' \item{h}{a vectorized function that takes time as input and returns the hazard at that time}
#' @export
#' @importFrom stats integrate
#' @examples con<-survdef(times = 3, surv = 0.5); survdefHR(con, 0.5)
survdefHR<-function(survdefC, HR){
    trthaz <- function(t) survdefC$h(t)*HR
    trtsur <- function(x) sapply(x, function(t) {
        if (t==0) 1
        else exp(-integrate(trthaz,lower=0,upper = t)$value)})
    return(list(S = trtsur, h = trthaz))
}
