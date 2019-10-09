#' Sample Size and Power for the Test of the Difference in Restricted Mean Survival Time
#'
#' Determine the asymptotic power of the test of RMST under a given trial design, or
#' calculate the samples size needed to achieve a desired power.
#'
#' @param survdefC the survival distribution of the control group, as a list in the form output by `survdef`.
#' @param survdefT the survival distribution of the treatment group, as a list in the form output by `survdef`.
#' @param k1 length of the accrual period. We assume subjects will accrue
#' uniformly over the interval `(0, k1)` and then be followed until trial time `k1+k2`.
#' @param k2 length of the follow-up period.
#' @param tau restriction time for the RMST.
#' @param n total sample size for both groups. 1:1 randomization is assumed.
#' Either `n` or `power` can be specified, and the other value will be calculated.
#' @param power the desired power.
#' @param plot if T, plots of the assumed survival distributions and power
#' as a function of sample size, accrual time ka1and follow-up time k2 will be produced.
#' Default is F.
#' The power of the RMST test, the log-rank test using all available followup and
#' the log-rank test using only followup to time tau after randomization will be
#' displayed. If two-sided=T, the power of the test for superiorty (treatment over
#'  control) and inferiority (control over treatment), are represented with solid
#'  and dashed lines, respectively.
#' @param sim if T, simulations will be conducted and empirical power and other
#'  summary statistics will be provided. Default is F. The hypothesis tests are
#'  carried out based on the normal approximation with the variance estimated
#'  according to the Greenwood plug-in/infinitesimal jackknife method.
#'  Specifying situations where survival doesn't go to zero in a reasonable
#'  amount of time (trial length times 1000) will lead to problems if the
#'  `sim=T` option is used.
#' @param M number of simualations. Default is 1000.
#' @param method modification to be used in simulations if the Kaplan-Meier
#' estimate is not defined at time `tau` in either group. Default is 'tau_star',
#' which changes the restriction time to the last censoring time,
#' if the last observation is censored at a time earlier than `tau`.
#' Other possible values are 'gill', 'efron', 'tau_star', 'risk1', 'risk2',
#' and 'risk5'. The riskX' options indicate estimating RMST difference at the
#' latest time at which at least X people are at risk in each group,
#' irrespective of the value of `tau`.
#' @param alpha type I error level. Default is 0.025 if 'two.sided'=F and 0.05
#' if 'two.sided'=T.
#' @param two.sided whether a two-sided test is desired. Default is F, meaning that
#' all reported power values correspond to a test of the superiority of treatment
#' over control. If set to T, the power for a test of superiority (treatment over control)
#'  and inferiority (control over treatment) will be reported separately in the results;
#' the power of a two-sided test is the sum of two.
#'
#' @return a list with components
#' \item{n}{the user-specified n, or if n was left blank, the n needed to achieve the user-specified power.}
#' \item{powerRMST}{the user-specified power, or if power was left blank, the asymptotic power of the RMST test.}
#' \item{powerRMSToverC}{the asymptotic power for a test of superiority of treatment over control.}
#' \item{powerRMSCoverT}{the asymptotic power for a test of superiority of control over treatment.
#'  If a one-sided test is specified, this is set to NA.}
#' \item{powerLRToverC}{the asymptotic power of the log-rank test of superiority
#'  of treatment over control.}
#' \item{powerLRCoverT}{the asymptotic power of the log-rank test of superiority
#'  of control over treatment. If a one-sided test is specified, this is set to NA.}
#' \item{powerLRtauToverC}{the asymptotic power of the log-rank test of superiority
#'  of treatment over control, using only data up to time tau after randomization.}
#' \item{powerLRtauCoverT}{the asymptotic power of the log-rank test of superiority
#'  of control over treatment, using only data up to time tau after randomization.
#'  If a one-sided test is specified, this is set to NA.}
#' \item{pKME}{the probability that you will be able to estimate RMST difference
#'  at time tau using the standard Kaplan-Meier estimator. If the last observation
#'  in either group is censored, and the censoring time is less than tau, the
#'  Kaplan-Meier estimate is not defined through time tau, and the RMST difference
#'  cannot be estimated using the standard area under the Kaplan-Meier curve. A
#'  modified estimator must be used.}
#' \item{simout}{a list returned if `sim = T`, with components:}
#'
#' \itemize{
#' \item{`emppowRMSTToverC`}{ empirical power of the RMST test for the superiority of
#'  treatment over control.}
#' \item{`emppowRMSTCoverT`}{ empirical power of the RMST test for the superiority of
#' control over treatment. If a one-sided test is specified, this is set to NA.}
#' \item{`emppowLRToverC`}{ empirical power of the log-rank test for the superiority of
#'  treatment over control.}
#' \item{`emppowLRCoverT`}{ empirical power of the log-rank testthe superiority of
#' control over treatment. If a one-sided test is specified, this is set to NA.}
#' \item{`emppowLRtauToverC`}{ empirical power of the log-rank test for the superiority of
#'  treatment over control, using only data up to time tau after randomization.}
#' \item{`emppowLRtauCoverT`}{ empirical power of the log-rank testthe superiority of
#' control over treatment, using only data up to time tau after randomization.
#' If a one-sided test is specified, this is set to NA.}
#' \item{`emppKME`}{ proportion of simulations where the standard KM estimator was used.}
#' \item{`meandiff`}{ mean estimated difference in RMST across the simulated datasets.}
#' \item{`SDdiff`}{ standard deviation of the estimated difference in RMST across the simulated datasets.}
#' \item{`meantrunc`}{ mean truncation time used in the simulated datasets (may be smaller than tau if method = 'tau_star' or 'riskX' options are used).}
#' \item{`SDtrunc`}{ standard deviation of the truncation time used in the simulated datasets.}
#' }
#' @export
#' @import survival
#' @import graphics
#' @import stats
#'
#' @examples
#' con<-survdef(times = 3, surv = 0.5)
#' trt<-survdef(haz = 0.67*con$h(1))
#' RMSTpow(con, trt, k1 = 0, k2 = 3, tau = 3, power = 0.8)
#' RMSTpow(con, trt, k1 = 0, k2 = 3, tau = 3, n = 552)
RMSTpow<-function(survdefC, survdefT, k1, k2, tau, n=NA, power=NA,
                  plot = F, sim = F, M = 1000, method = 'tau_star',
                  alpha = NA, two.sided=F){
    if (is.na(alpha)) alpha <- ifelse(two.sided==T, 0.05, 0.025)
    if (is.na(power)+is.na(n)!=1) stop('One of n, power must be missing.')
    if (tau<=0) stop('Tau must be greater than zero.')
    if (is.na(n)){
        RMST_truediff<-integrate(function(x) survdefT$S(x)-survdefC$S(x),
                                 lower = 0, upper = tau)$value
        if (two.sided==F) {
            if (RMST_truediff<0) stop('True RMST in treatment arm is less than true RMST in control arm; cannot design a trial to show treatment is superior.')
            RMST_trueSE <- RMST_truediff/(qnorm(1-alpha)-qnorm(1-power))
            find_n<-function(N) sqrt(evar(survdefT, k1, k2, tau, k1+k2, N)+
                                     evar(survdefC, k1, k2, tau, k1+k2, N))-RMST_trueSE
            if (find_n(10000)>0) stop('Trial size would be more than 20,000 patients; please select different design parameters.')
        }
        else {
            find_n<-function(N) 1-pnorm(qnorm(1-alpha/2)-RMST_truediff/sqrt(evar(survdefT, k1, k2, tau, k1+k2, N)+evar(survdefC, k1, k2, tau, k1+k2, N)))+
                pnorm(qnorm(alpha/2)-RMST_truediff/sqrt(evar(survdefT, k1, k2, tau, k1+k2, N)+evar(survdefC, k1, k2, tau, k1+k2, N)))-
                power
            if (find_n(10000)<0) stop('Trial size would be more than 20,000 patients; please select different design parameters.')
        }
        N<- uniroot(find_n, lower = 1, upper = 10000)$root
        n <- 2*ceiling(N)
    }
    else{if (n%%2 !=0) n<- ifelse(floor(n)%%2==0,floor(n),floor(n)-1)}
#    powerRMST <- powfn(survdefC, survdefT, k1, k2, tau, n, alpha)
#    powerLR <- LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha)
#    powerLRtau <- LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha, tau = tau)
    if (two.sided==F) {
        powerRMSTToverC <- powerRMST <- powfn(survdefC, survdefT, k1, k2, tau, n, alpha)
        powerRMSTCoverT <- NA
        powerLRToverC <- LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha)
        powerLRCoverT <- NA
        powerLRtauToverC <- LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha, tau = tau)
        powerLRtauCoverT <- NA
    }
    else{
        powerRMSTToverC<-powfn(survdefC, survdefT, k1, k2, tau, n, alpha/2)
        powerRMSTCoverT<-powfn(survdefT, survdefC, k1, k2, tau, n, alpha/2)
        powerRMST<-powerRMSTToverC+powerRMSTCoverT
        powerLRToverC <- LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha/2)
        powerLRCoverT <- LRpow_nonPH(survdefT, survdefC, k1,k2,n, alpha/2)
        powerLRtauToverC <- LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha/2, tau = tau)
        powerLRtauCoverT <- LRpow_nonPH(survdefT, survdefC, k1,k2,n, alpha/2, tau = tau)
    }
    pKME <- RMSTeval(survdefC, survdefT, k1, k2, tau, n)
    to_ret<-list(n=n, powerRMST=powerRMST, powerRMSTToverC=powerRMSTToverC,powerRMSTCoverT=powerRMSTCoverT,
                 powerLRToverC=powerLRToverC, powerLRCoverT=powerLRCoverT,
                 powerLRtauToverC=powerLRtauToverC, powerLRtauCoverT=powerLRtauCoverT,
                 pKME=pKME)
    if (plot == T){
        par(mfrow = c(2,2), mar = c(5,4,1,2))
        n_vec<-seq(from = n*.5, to = n*1.5,length.out =  20)
        k1_vec<-seq(from = max(0,tau-k2), to = k1*1.5,length.out =  20)
        if (k1==0) k1_vec<-seq(from = 0, to = k2,length.out =  20)
        k2_vec<-seq(from = max(0,tau-k1), to = k2*1.5,length.out =  20)
        if (k2==0) k2_vec<-seq(from = 0, to = k1,length.out =  20)
        tau_vec<-seq(from = tau*.5, to = min(k1+k2,tau*1.5),length.out =  20)

        plotsurvdef(survdefC, survdefT, xupper = k1+k2)
        legend((k1+k2)*.5, 0.25, c('Control','Treatment'),cex = .5, lty = c(1,2), bty='n')

        if (two.sided==F) {
            plot(n_vec, sapply(n_vec, function(x) powfn(survdefC, survdefT, k1, k2, tau, x, alpha)),
                 xlab = 'Sample size', ylab = 'Power', ylim = c(0,1), type ='l')
            lines(n_vec, sapply(n_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,x, alpha)), col = 'red')
            lines(n_vec, sapply(n_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,x, alpha, tau)), col = 'blue')
            legend(n_vec[10], 0.25, c('RMST','Log-rank', 'Log-rank (tau)'), col = c(1,2,4),cex = .5, lty = 1, bty='n')

            plot(k1_vec, sapply(k1_vec, function(x) powfn(survdefC, survdefT, x, k2, tau, n, alpha)),
                 xlab = 'Accrual time', ylab = 'Power', ylim = c(0,1), type ='l')
            lines(k1_vec, sapply(k1_vec, function(x) LRpow_nonPH(survdefC, survdefT, x,k2,n, alpha)), col = 'red')
            lines(k1_vec, sapply(k1_vec, function(x) LRpow_nonPH(survdefC, survdefT, x,k2,n, alpha, tau)), col = 'blue')
            legend(k1_vec[10], 0.25, c('RMST','Log-rank', 'Log-rank (tau)'), col = c(1,2,4),cex = .5, lty = 1, bty='n')

            plot(k2_vec, sapply(k2_vec, function(x) powfn(survdefC, survdefT, k1, x, tau, n, alpha)),
                 xlab = 'Followup time', ylab = 'Power', ylim = c(0,1), type ='l')
            lines(k2_vec, sapply(k2_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,x,n, alpha)), col = 'red')
            lines(k2_vec, sapply(k2_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,x,n, alpha, tau)), col = 'blue')
            legend(k2_vec[10], 0.25, c('RMST','Log-rank', 'Log-rank (tau)'), col = c(1,2,4),cex = .5, lty = 1, bty='n')

           # plot(tau_vec, sapply(tau_vec, function(x) powfn(survdefC, survdefT, k1, k2, x, n, alpha)),
        #         xlab = 'Tau', ylab = 'Power', ylim = c(0,1), type ='l')
         #   lines(tau_vec, sapply(tau_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha)), col = 'red')
          #  lines(tau_vec, sapply(tau_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha, x)), col = 'blue')
          #  legend(tau_vec[10], 0.25, c('RMST','Log-rank', 'Log-rank (tau)'), col = c(1,2,4),cex = .5, lty = 1, bty='n')
            }
        else {
            plot(n_vec, sapply(n_vec, function(x) powfn(survdefC, survdefT, k1, k2, tau, x, alpha/2)),
                 xlab = 'Sample size', ylab = 'Power', ylim = c(0,1), type ='l')
            lines(n_vec, sapply(n_vec, function(x) powfn(survdefT, survdefC, k1, k2, tau, x, alpha/2)), lty=2)
            lines(n_vec, sapply(n_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,x, alpha/2)), col = 'red')
            lines(n_vec, sapply(n_vec, function(x) LRpow_nonPH(survdefT, survdefC, k1,k2,x, alpha/2)), col = 'red', lty = 2)
            lines(n_vec, sapply(n_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,x, alpha/2, tau)), col = 'blue')
            lines(n_vec, sapply(n_vec, function(x) LRpow_nonPH(survdefT, survdefC, k1,k2,x, alpha/2, tau)), col = 'blue', lty = 2)
            legend(n_vec[10], 0.25, c('RMST','Log-rank', 'Log-rank (tau)'), col = c(1,2,4),cex = .5, lty = 1, bty='n')

            plot(k1_vec, sapply(k1_vec, function(x) powfn(survdefC, survdefT, x, k2, tau, n, alpha/2)),
                 xlab = 'Accrual time', ylab = 'Power', ylim = c(0,1), type ='l')
            lines(k1_vec, sapply(k1_vec, function(x) powfn(survdefT, survdefC, x, k2, tau, n, alpha/2)), lty=2)
            lines(k1_vec, sapply(k1_vec, function(x) LRpow_nonPH(survdefC, survdefT, x,k2,n, alpha/2)), col = 'red')
            lines(k1_vec, sapply(k1_vec, function(x) LRpow_nonPH(survdefT, survdefC, x,k2,n, alpha/2)), col = 'red', lty = 2)
            lines(k1_vec, sapply(k1_vec, function(x) LRpow_nonPH(survdefC, survdefT, x,k2,n, alpha/2, tau)), col = 'blue')
            lines(k1_vec, sapply(k1_vec, function(x) LRpow_nonPH(survdefT, survdefC, x,k2,n, alpha/2, tau)), col = 'blue', lty = 2)
            legend(k1_vec[10], 0.25, c('RMST','Log-rank', 'Log-rank (tau)'), col = c(1,2,4),cex = .5, lty = 1, bty='n')

            plot(k2_vec, sapply(k2_vec, function(x) powfn(survdefC, survdefT, k1, x, tau, n, alpha/2)),
                 xlab = 'Followup time', ylab = 'Power', ylim = c(0,1), type ='l')
            lines(k2_vec, sapply(k2_vec, function(x) powfn(survdefT, survdefC, k1, x, tau, n, alpha/2)), lty=2)
            lines(k2_vec, sapply(k2_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,x,n, alpha/2)), col = 'red')
            lines(k2_vec, sapply(k2_vec, function(x) LRpow_nonPH(survdefT, survdefC, k1,x,n, alpha/2)), col = 'red', lty = 2)
            lines(k2_vec, sapply(k2_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,x,n, alpha/2, tau)), col = 'blue')
            lines(k2_vec, sapply(k2_vec, function(x) LRpow_nonPH(survdefT, survdefC, k1,x,n, alpha/2, tau)), col = 'blue', lty = 2)
            legend(k2_vec[10], 0.25, c('RMST','Log-rank', 'Log-rank (tau)'), col = c(1,2,4),cex = .5, lty = 1, bty='n')

          #  plot(tau_vec, sapply(tau_vec, function(x) powfn(survdefC, survdefT, k1, k2, x, n, alpha/2)),
          #       xlab = 'Tau', ylab = 'Power', ylim = c(0,1), type ='l')
          #  lines(tau_vec, sapply(tau_vec, function(x) powfn(survdefT, survdefC, k1, k2, x, n, alpha/2)), lty=2)
          #  lines(tau_vec, sapply(tau_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha/2)), col = 'red')
          #  lines(tau_vec, sapply(tau_vec, function(x) LRpow_nonPH(survdefT, survdefC, k1,k2,n, alpha/2)), col = 'red', lty = 2)
          #  lines(tau_vec, sapply(tau_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha/2, x)), col = 'blue')
          #  lines(tau_vec, sapply(tau_vec, function(x) LRpow_nonPH(survdefT, survdefC, k1,k2,n, alpha/2, x)), col = 'blue', lty = 2)
          #  legend(tau_vec[10], 0.25, c('RMST','Log-rank', 'Log-rank (tau)'), col = c(1,2,4),cex = .5, lty = 1, bty='n')
        }
    }
    if (sim == T) {
        if (n%%2 !=0) stop('Total sample size n must be even to do simulations.')
        print('Simulating datasets...')
        simout<-list(emppowRMSTToverC=NA, emppowRMSTCoverT=NA,
                     emppowLRToverC=NA, emppowLRCoverT=NA,
                     emppowLRtauToverC=NA, emppowLRtauCoverT=NA,
                     emppKME=NA, meandiff=NA, SDdiff=NA, meantrunc=NA, SDtrunc=NA)
        Fc_inv<-function(u){
            findroot<-function(x) 1-survdefC$S(x)-u
            return(uniroot(findroot, lower = 0, upper = (k1+k2)*1000)$root)
        }
        Ft_inv<-function(u){
            findroot<-function(x) 1-survdefT$S(x)-u
            return(uniroot(findroot, lower = 0, upper = (k1+k2)*1000)$root)
        }
        out<-data.frame(matrix(NA, nrow = M, ncol = 10))
        colnames(out)<-c('last_censT','last_censC','tau_temp','RMST_diff',
                         'RMST_rejectToverC','RMST_rejectCoverT',
                         'LR_rejectToverC','LR_rejectCoverT',
                         'LRtau_rejectToverC','LRtau_rejectCoverT')
        for (m in 1:M){
            out$tau_temp[m]<-tau
            grp<-c(rep('c', n/2), rep('t', n/2))
            enrol<-runif(n, min=0, max = k1)
            con<-sapply(runif(n/2, min=0, max=1),Fc_inv)
            trt<-sapply(runif(n/2, min=0, max=1),Ft_inv)
            event<-c(con, trt)+enrol< k1+k2
            time<-pmin(c(con, trt),  k1+k2-enrol)
            surv_res<-survfit(Surv(time, event)~grp, timefix = F)

            c.time<-surv_res$time[1:surv_res$strata[1]]
            c.cens<-surv_res$n.censor[1:surv_res$strata[1]]
            c.ev<-surv_res$n.event[1:surv_res$strata[1]]
            c.n.risk<-surv_res$n.risk[1:surv_res$strata[1]]
            t.time<-surv_res$time[surv_res$strata[1]+1:surv_res$strata[2]]
            t.cens<-surv_res$n.censor[surv_res$strata[1]+1:surv_res$strata[2]]
            t.ev<-surv_res$n.event[surv_res$strata[1]+1:surv_res$strata[2]]
            t.n.risk<-surv_res$n.risk[surv_res$strata[1]+1:surv_res$strata[2]]

            out$last_censC[m]<-ifelse(c.cens[length(c.cens)]==0 | tau <= max(c.time), F, T)
            out$last_censT[m]<-ifelse(t.cens[length(t.cens)]==0 | tau <= max(t.time), F, T)

            if (method == 'tau_star' & (out$last_censC[m]|out$last_censT[m])){
                out$tau_temp[m]<-min(c(max(t.time), max(c.time)), na.rm =T )
            }
            else if (method == 'efron' & (out$last_censC[m]|out$last_censT[m])){
                if (out$last_censC[m]) event[grp=='c'][which.max(time[grp=='c'])]<-T
                if (out$last_censT[m]) event[grp=='t'][which.max(time[grp=='t'])]<-T
                surv_res<-survfit(Surv(time, event)~grp, timefix = F)
            }
            else if (method %in% c('risk1', 'risk2', 'risk5')){
                if (method == 'risk1') out$tau_temp[m]<-min(c(max(t.time), max(c.time)), na.rm =T )
                if (method == 'risk2') out$tau_temp[m]<-min(c(max(t.time[t.n.risk>=2]),
                                                              max(c.time[c.n.risk>=2])), na.rm =T )
                if (method == 'risk5') out$tau_temp[m]<-min(c(max(t.time[t.n.risk>=5]),
                                                              max(c.time[c.n.risk>=5])), na.rm =T )
            }

            summ<-summary(surv_res, rmean = out$tau_temp[m])$table
            out$RMST_diff[m]<-summ[2,5]-summ[1,5]
            survtest<-survdiff(Surv(time, event)~grp)#, timefix = F)
            time2<-pmin(time, tau)
            event2<-ifelse(time<=tau, event, F)
            survtest2<-survdiff(Surv(time2, event2)~grp)#, timefix = FALSE)
            if (two.sided == F) {
                out$RMST_rejectToverC[m]<-out$RMST_diff[m]/(sqrt(summ[2,6]^2 + summ[1,6]^2))>qnorm(1-alpha)
                out$RMST_rejectCoverT[m]<-NA
                out$LR_rejectToverC[m]<-((1 - pchisq(survtest$chisq, 1))<alpha &
                                       survtest$obs[1]>survtest$exp[1])
                out$LR_rejectCoverT[m]<-NA
                out$LRtau_rejectToverC[m]<-((1 - pchisq(survtest2$chisq, 1))<alpha &
                                          survtest2$obs[1]>survtest2$exp[1])
                out$LRtau_rejectCoverT<-NA
            }
            else {
                out$RMST_rejectToverC[m]<-out$RMST_diff[m]/(sqrt(summ[2,6]^2 + summ[1,6]^2))>qnorm(1-alpha/2)
                out$RMST_rejectCoverT[m]<-(out$RMST_diff[m])/(sqrt(summ[2,6]^2 + summ[1,6]^2))<qnorm(alpha/2)
                out$LR_rejectToverC[m]<-((1 - pchisq(survtest$chisq, 1))<alpha/2 &
                                             survtest$obs[1]>survtest$exp[1])
                out$LR_rejectCoverT[m]<-((1 - pchisq(survtest$chisq, 1))<alpha/2 &
                                             survtest$obs[1]<survtest$exp[1])
                out$LRtau_rejectToverC[m]<-((1 - pchisq(survtest2$chisq, 1))<alpha/2 &
                                                survtest2$obs[1]>survtest2$exp[1])
                out$LRtau_rejectCoverT[m]<-((1 - pchisq(survtest2$chisq, 1))<alpha/2 &
                                                survtest2$obs[1]<survtest2$exp[1])
            }
        }
        simout$emppowRMSTToverC<-sum(out$RMST_rejectToverC)/M
        simout$emppowRMSTCoverT<-sum(out$RMST_rejectCoverT)/M
        simout$emppowLRToverC<-sum(out$LR_rejectToverC)/M
        simout$emppowLRCoverT<-sum(out$LR_rejectCoverT)/M
        simout$emppowLRtauToverC<-sum(out$LRtau_rejectToverC)/M
        simout$emppowLRtauCoverT<-sum(out$LRtau_rejectCoverT)/M
        simout$emppKME<-sum(out$last_censT==F & out$last_censC==F)/M
        simout$meandiff<-mean(out$RMST_diff)
        simout$SDdiff<-sd(out$RMST_diff)
        simout$meantrunc<-mean(out$tau_temp, na.rm = T)
        simout$SDtrunc<-sd(out$tau_temp, na.rm = T)
        to_ret[['simout']]<-simout
    }
    return(to_ret)
}
