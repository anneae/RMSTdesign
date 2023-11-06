#' Sample Size and Power for the Test of the Difference in Restricted Mean Survival Time Under Biomarker Misclassification
#'
#' Determine the asymptotic power of the test of RMST under a given trial design when using
#' estimators that are adjusted for biomarker misclassification (i.e., the sensitivity and/or specificity of
#' the biomarker test is less than 1), or
#' calculate the samples size needed to achieve a desired power in this setting.
#'
#' @param survdef_t0_m0 the survival distribution of the control, marker negative group, as a list in the form output by `survdef`.
#' @param survdef_t1_m0 the survival distribution of the treatment, marker negative group, as a list in the form output by `survdef`.
#' @param survdef_t0_m1 the survival distribution of the control, marker positive group, as a list in the form output by `survdef`.
#' @param survdef_t1_m1 the survival distribution of the treatment, marker positive group, as a list in the form output by `survdef`.
#' @param k1 length of the accrual period. We assume subjects will accrue
#' uniformly over the interval `(0, k1)` and then be followed until trial time `k1+k2`.
#' @param k2 length of the follow-up period.
#' @param tau restriction time for the RMST.
#' @param prev the prevalence of the biomarker.
#' @param sens the sensitivity of the biomarker test.
#' @param spec the specificity of the biomarker test.
#' @param n total sample size for the clinical trial. Either `n` or `power` can
#' be specified, and the other value will be calculated. 1:1 randomization is assumed
#' as well as enrollment of equal numbers of patients who test positive for the biomarker and
#' patients who test negative for the biomarker, i.e., n/4 patients who test positive will be
#' randomized to the treatment group, n/4 patients who test negative for the biomarker will
#' be randomized to the control group, etc.
#' @param power the desired power.
#' @param alpha type I error level. Default is 0.025 if 'two.sided'=F and 0.05
#' if 'two.sided'=T.
#' @param two.sided whether a two-sided test is desired. Default is F.
#' If test = 'T in M+' or 'T in M-', this corresponds to a test of treatment over control;
#' if test = 'M in T+' or 'M in T-', it corresponds to a test of marker positive over marker negative;
#' if test = 'interaction' it corresponds to a test to show a larger treatment effect in the marker positive group.
#' If two.sided = T, the power for each one-sided test will be reported separately in the results;
#' the power of a two-sided test is the sum of two.
#' @param test the test for which you want to power the clinical trial This can be set to
#' 'T in M+' (treatment effect in marker positive group), 'T in M-' (treatment effect in marker negative group),
#'  'M in T+' (marker effect in treatment group), 'M in T-' (marker effect in control group) or
#'  'interaction'.
#'
#' @return a list with components
#' \item{n}{the user-specified n, or if n was left blank, the n needed to achieve the user-specified power.}
#' \item{powerRMST}{the user-specified power, or if power was left blank, the asymptotic power of the RMST test.
#' If `one-sided=T`, `powerRMST` is equivalent to `powerRMST1over0`.
#' If `one-sided=F`, `powerRMST` is equivalent to the sum of the power of a one-sided test in each direction, i.e.
#' `powerRMST1over0 + powerRMST0over1`.}
#' \item{powerRMST1over0}{the asymptotic power for a test of treatment over
#' control or marker positive over marker negative or a larger treatment effect
#' in the marker positive group. See the explanations of the `test` and `two.sided` parameters.}
#' \item{powerRMST0over1}{the asymptotic power for a test of control over
#' treatment or marker negative over marker positive or a larger treatment effect
#' in the marker negative group. See the explanations of the `test` and `two.sided` parameters.
#' If a one-sided test is specified, this is set to NA.}
#' \item{pKME}{the probability that you will be able to estimate RMST difference
#'  at time tau in the 4 groups (treatment/control and test positive/negative)
#'  using the standard Kaplan-Meier estimator. If the last observation
#'  in any group is censored, and the censoring time is less than tau, the
#'  Kaplan-Meier estimate is not defined through time tau, and the RMST difference
#'  cannot be estimated using the standard area under the Kaplan-Meier curve. A
#'  modified estimator must be used.}
#' @export
#'
#' @examples
#' surv_t0_m0 <- survdefWeibull(shape = .8, scale = 10)
#' surv_t1_m0 <- survdefWeibull(shape = .8, scale = exp(-.5)^(-1/.8)*10)
#' surv_t0_m1 <- survdefWeibull(shape = .8, scale = exp(.1)^(-1/.8)*10)
#' surv_t1_m1 <- survdefWeibull(shape = .8, scale = exp(-.5+.1+.3)^(-1/.8)*10)
#' RMSTpowadjusted(surv_t0_m0, surv_t1_m0, surv_t0_m1,surv_t1_m1,
#'                  k1=20, k2=5, tau = 10, prev=.3, sens=.8, spec=.8,
#'                  n=400, power=NA, alpha = .05, two.sided=T, test = 'T in M+')
RMSTpowadjusted <-function(survdef_t0_m0, survdef_t1_m0,survdef_t0_m1,survdef_t1_m1,
                            k1, k2, tau,
                            prev, sens, spec,
                            n=NA, power=NA,
                            alpha = NA, two.sided=F, test = 'T in M+'){

    if (is.na(alpha)) alpha <- ifelse(two.sided==T, 0.05, 0.025)
    if (is.na(power)+is.na(n)!=1) stop('One of n, power must be missing.')
    if (tau<=0) stop('Tau must be greater than zero.')
    if (!(test %in% c('T in M+','T in M-','M in T+','M in T-','interaction'))) stop('Test must be one of the following: T in M+, T in M-, M in T+, M in T-, interaction.')
    if (sens <0|sens>1|spec<0|spec>1) stop('Sensitivity and specificity must be between 0 and 1.')

    ppv<-prev*sens/(prev*sens+(1-prev)*(1-spec))
    npv<-(1-prev)*spec/(prev*(1-sens)+(1-prev)*spec)

    if (test == 'T in M+') RMST_truediff<-integrate(function(x) survdef_t1_m1$S(x)-survdef_t0_m1$S(x),
                                                    lower = 0, upper = tau)$value
    if (test == 'T in M-') RMST_truediff<-integrate(function(x) survdef_t1_m0$S(x)-survdef_t0_m0$S(x),
                                                    lower = 0, upper = tau)$value
    if (test == 'M in T+') RMST_truediff<-integrate(function(x) survdef_t1_m1$S(x)-survdef_t1_m0$S(x),
                                                    lower = 0, upper = tau)$value
    if (test == 'M in T-') RMST_truediff<-integrate(function(x) survdef_t0_m1$S(x)-survdef_t0_m0$S(x),
                                                    lower = 0, upper = tau)$value
    if (test == 'interaction') RMST_truediff<-integrate(function(x) survdef_t1_m1$S(x)-survdef_t0_m1$S(x)-(survdef_t1_m0$S(x)-survdef_t0_m0$S(x)),
                                                        lower = 0, upper = tau)$value

    # Need the variance of the naive estimators to calculate the variance of the adjusted estimators
    evar.t1m1.naive <- function(N) evar(survdefMix(survdef_t1_m1, survdef_t1_m0, ppv, 1-ppv), k1, k2, tau, k1+k2, N)
    evar.t1m0.naive <- function(N) evar(survdefMix(survdef_t1_m0, survdef_t1_m1, npv, 1-npv), k1, k2, tau, k1+k2, N)
    evar.t0m1.naive <- function(N) evar(survdefMix(survdef_t0_m1, survdef_t0_m0, ppv, 1-ppv), k1, k2, tau, k1+k2, N)
    evar.t0m0.naive <- function(N) evar(survdefMix(survdef_t0_m0, survdef_t0_m1, npv, 1-npv), k1, k2, tau, k1+k2, N)
    if (test == 'T in M+') myevar <- function(N) (npv^2*evar.t1m1.naive(N)+(1-ppv)^2*evar.t1m0.naive(N))/((npv*ppv-(1-npv)*(1-ppv))^2)+
        (npv^2*evar.t0m1.naive(N)+(1-ppv)^2*evar.t0m0.naive(N))/((npv*ppv-(1-npv)*(1-ppv))^2)
    if (test == 'T in M-') myevar <- function(N) (ppv^2*evar.t1m0.naive(N)+(1-npv)^2*evar.t1m1.naive(N))/((npv*ppv-(1-npv)*(1-ppv))^2)+
        (ppv^2*evar.t0m0.naive(N)+(1-npv)^2*evar.t0m1.naive(N))/((npv*ppv-(1-npv)*(1-ppv))^2)
    if (test == 'M in T+') myevar <- function(N) (evar.t1m1.naive(N)+evar.t1m0.naive(N))/((npv*ppv-(1-npv)*(1-ppv))^2)
    if (test == 'M in T-') myevar <- function(N) (evar.t0m1.naive(N)+evar.t0m0.naive(N))/((npv*ppv-(1-npv)*(1-ppv))^2)
    if (test == 'interaction') myevar <- function(N) (evar.t1m1.naive(N)+evar.t1m0.naive(N)+evar.t0m1.naive(N)+evar.t0m0.naive(N))/((npv*ppv-(1-npv)*(1-ppv))^2)

    if (is.na(n)){
        if (two.sided==F) {
            if (RMST_truediff<0) stop('True RMST in treatment/marker + arm is less than true RMST in treatment/marker - arm or interaction effect is negative; cannot design a trial to show a positive effect.')
            RMST_trueSE <- RMST_truediff/(qnorm(1-alpha)-qnorm(1-power))
            find_n<-function(N) sqrt(myevar(N))-RMST_trueSE
            if (find_n(10000)>0) stop('Trial size would be more than 10,000 patients per group; please select different design parameters.')
        }
        else {
            find_n<-function(N) 1-pnorm(qnorm(1-alpha/2)-RMST_truediff/sqrt(myevar(N)))+
                pnorm(qnorm(alpha/2)-RMST_truediff/sqrt(myevar(N)))-
                power
            if (find_n(10000)<0) stop('Trial size would be more than 10,000 patients per group; please select different design parameters.')
        }
        N<- uniroot(find_n, lower = 1, upper = 10000)$root
        n <- 4*ceiling(N) # total trial sample size
    }
    else{if (n%%4 !=0) n<- ifelse(floor(n)%%4==0,floor(n),floor(n)+floor(n)%%4)}
    if (two.sided==F) {
        powerRMST0over1 <- NA
        powerRMST1over0 <- powerRMST <- 1-pnorm(qnorm(1-alpha)-RMST_truediff/sqrt(myevar(n/4)))
    }
    else{
        powerRMST0over1 <- pnorm(qnorm(alpha/2)-RMST_truediff/sqrt(myevar(n/4)))
        powerRMST1over0 <- 1-pnorm(qnorm(1-alpha/2)-RMST_truediff/sqrt(myevar(n/4)))
        powerRMST<-powerRMST1over0+powerRMST0over1
    }
    pKME <- RMSTeval(survdefMix(survdef_t1_m1, survdef_t1_m0, ppv, 1-ppv), survdefMix(survdef_t1_m0, survdef_t1_m1, npv, 1-npv), k1, k2, tau, n/2)*
        RMSTeval(survdefMix(survdef_t0_m1, survdef_t0_m0, ppv, 1-ppv), survdefMix(survdef_t0_m0, survdef_t0_m1, npv, 1-npv), k1, k2, tau, n/2)

    to_ret<-list(test = test, n=n, powerRMST=powerRMST, powerRMST1over0=powerRMST1over0,powerRMST0over1=powerRMST0over1,pKME=pKME)
    return(to_ret)
}
