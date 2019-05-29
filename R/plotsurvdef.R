#' Plots two survival distributions.
#'
#'
#' @param survdefC the survival distribution of the control group (will be
#' plotted as a solid line), as a list in the form output by `survdef`.
#' @param survdefT the survival distribution of the control group (will be
#'  plotted as a dashed line), as a list in the form output by `survdef`.
#' @param xupper the upper x axis limit for the plot.
#'
#' @export
#'
#' @examples RIC<-survdef(times = 18, surv = .6)
#' MAC<-survdef(times = 3, haz=c(4.375*RIC$h(1),RIC$h(1)))
#' plot.survdef(RIC, MAC, 24)
plotsurvdef<-function(survdefC, survdefT, xupper){
    xvec<-seq(from = 0, to = xupper, length.out = 100)
    plot(xvec, survdefC$S(xvec),ylim = c(0,1), type = 'l',
         xlab = 'Time', ylab = 'Survival')
    lines(xvec, survdefT$S(xvec), lty = 2)
}
