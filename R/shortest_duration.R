#' Find Shortest Duration Trial
#'
#' Find the trial with the shortest duration in calendar time with a specified
#' power and probability that RMST difference will be estimable with the
#' Kaplan-Meier estimator. Within all trials of that minimum duration, the
#' function finds the one with the smallest sample size.
#'
#' @param survdefC the survival distribution of the control group, as a list in the form output by survdef.
#' @param survdefT the survival distribution of the treatment group, as a list in the form output by survdef.
#' @param tau restriction time for the RMST.
#' @param power the desired power.
#' @param accrual_rate the planned accrual rate, per unit of time.
#' @param pKME The desired probability that the RMST difference will be
#' estimable using the Kaplan-Meier estimator. Default is 0.95.
#' @param alpha one-sided type I error level. Default is 0.025.
#' @param altdesign if T, the function will look for an alternative design that
#' is not the shortest in duration, but has duration equal to some multiple of
#' the shortest possible duration. Default is F. The sample size of the shortest
#'  duration trial can be much larger than a slightly longer trial, so we
#'  recommend considering an alternative design slightly longer than the
#'  shortest trial in addition to the shortest trial.
#' @param multiplier  the factor by which the duration of the shortest possible
#' trial is multiplied to acquire the duration of the alternative trial design.
#' Default is 1.1, meaning a trial that is 10% longer in duration than the
#' shortest possible trial. This argument is ignored if altdesign=F.
#'
#' @return a list with components:
#' \item{k1}{length of the accrual period. We assume subjects will accrue
#' uniformly over the interval (0, k1) and then be followed until trial time
#' k1+k2.}
#' \item{k2}{length of the follow-up period.  }
#' \item{duration}{trial duration in calendar time, k1+k2.  }
#' \item{powerRMST}{the asymptotic power of the RMST test.  }
#' \item{powerLR}{the asymptotic power of the log-rank test.  }
#' \item{pKME}{the probability that you will be able to estimate RMST
#' difference at time tau using the standard Kaplan-Meier estimator.
#' @export
#'
#' @examples
#' con<-survdef(times = 3, surv = 0.5)
#' trt<-survdef(haz = 0.67*con$h(1))
#' shortest_duration(con, trt, 3, .8, 552/4)
shortest_duration<-function(survdefC, survdefT, tau, power, accrual_rate, pKME=.95,
                            alpha = 0.025, altdesign = F, multiplier=1.1){
    RMST_truediff<-integrate(function(x) survdefT$S(x)-survdefC$S(x),
                             lower = 0, upper = tau)$value
    RMST_trueSE <- RMST_truediff/(qnorm(1-alpha)-qnorm(1-power))
    powtau<-RMSTpow(survdefC, survdefT,k1 = tau, k2 = 0,tau = tau, n = tau*accrual_rate)$powerRMST
    pkme_tau<-RMSTeval(survdefC, survdefT, tau, 0, tau, tau*accrual_rate)
    if (powtau>=power & pkme_tau>=pKME) return(tau)
    # The shortest trial will be the one that accrues the whole time. We'll find
    # that, then see how much we can shorten k1 while keeping trial duration constant.
    # First, find the trial that accrues the whole time that fits the first criteria
    if (powtau<power){
        find_k1<-function(k1) sqrt(evar(survdefT, k1, 0, tau, k1, k1*accrual_rate/2)+
                                       evar(survdefC, k1, 0, tau, k1, k1*accrual_rate/2))-RMST_trueSE
        k1_cand<-uniroot(find_k1, lower = tau, upper = tau*20)$root}
    else k1_cand<-tau
    # Check whether it fits the second criteria
    if (pkme_tau>=pKME) trial_length<-k1_cand
    else if (RMSTeval(survdefC, survdefT, k1_cand, 0, tau, k1_cand*accrual_rate)>=pKME) trial_length<-k1_cand
    else {
        find_k1<-function(k1) RMSTeval(survdefC, survdefT, k1, 0, tau, k1*accrual_rate)-pKME
        trial_length<-uniroot(find_k1, lower = k1_cand, upper = tau*20)$root
    } # trial_length is the shortest trial duration that fits both criteria.
    # Next find the shortest accrual that would keep power high
    if (altdesign==T) {
        if (multiplier<1) stop('multiplier must be >1.')
        trial_length<-multiplier*trial_length
    }
    k1_pow<-trial_length
    if(sqrt(evar(survdefT, trial_length, 0, tau, trial_length, trial_length*accrual_rate/2)+
            evar(survdefC, trial_length, 0, tau, trial_length, trial_length*accrual_rate/2))<RMST_trueSE)
    {
        find_k1<-function(k1) sqrt(evar(survdefT, k1, trial_length-k1, tau, trial_length, k1*accrual_rate/2)+
                                       evar(survdefC, k1, trial_length-k1, tau, trial_length, k1*accrual_rate/2))-RMST_trueSE
        k1_pow<-uniroot(find_k1, lower = 10/accrual_rate, upper = trial_length)$root
    }
    # Find the shortest accrual that would keep pKME high enough
    k1_pkme<-trial_length
    if (RMSTeval(survdefC, survdefT, trial_length, 0, tau, trial_length*accrual_rate)>pKME)
    {
        if (RMSTeval(survdefC, survdefT, 4/accrual_rate, trial_length-4/accrual_rate, tau, 4)>pKME) k1_pkme<-4/accrual_rate
        else {
            find_k1<-function(k1) RMSTeval(survdefC, survdefT, k1, trial_length-k1, tau, k1*accrual_rate)-pKME
            k1_pkme<-uniroot(find_k1, lower = 4/accrual_rate, upper = trial_length)$root
        }
    }
    # Set k1 to the longer of the two
    k1<-max(k1_pow, k1_pkme)
    n<-ceiling(k1*accrual_rate)
    if ((n %% 2) == 1) n<-n+1 # Get integer n, the smallest even number >= n
    k1<-n/accrual_rate
    k2<-max(0,trial_length-k1)
    return(c(list(k1=k1, k2 =k2, duration = k1+k2),
             RMSTpow(survdefC, survdefT, k1, k2, tau = tau, n=n)))
}
