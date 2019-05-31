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
#' @param plot if T, plots of power versus each design parameter will be produced. Default is F.
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
#' @param alpha one-sided type I error level. Default is 0.025.
#'
#' @return a list with components
#' \item{n}{the user specified n, or if n was left blank, the n needed to achieve the user-specified power.}
#' \item{powerRMST}{the user specified power, or if power was left blank, the asymptotic power of the RMST test.}
#' \item{powerLR}{the asymptotic power of the log-rank test.}
#' \item{pKME}{the probability that you will be able to estimate RMST difference
#'  at time tau using the standard Kaplan-Meier estimator. If the last observation
#'  in either group is censored, and the censoring time is less than tau, the
#'  Kaplan-Meier estimate is not defined through time tau, and the RMST difference
#'  cannot be estimated using the standard area under the Kaplan-Meier curve. A
#'  modified estimator must be used.}
#' \item{plotvals}{a list used for generating plots, returned if `plot = T`,
#' with components:}
#' \itemize{
#' \item{`n_vec`}{ a vector of sample sizes near the selected sample size}
#' \item{`n_pow`}{ power of the RMST based test for the values in n_vec, holding k1, k2 and tau constant  }
#' \item{`k1_vec`}{ a vector of accrual periods near the specified k1}
#' \item{`k1_pow`}{ power of the RMST based test for the values in k1_vec, holding n, tau and total trial length constant}
#' \item{`k2_vec`}{ a vector of follow-up periods near the specified k2}
#' \item{`k2_pow`}{ power of the RMST based test for the values in k2_vec, holding n, k1 and tau constant  }
#' \item{`tau_vec`}{ a vector of restriction times near the specified tau}
#' \item{`tau_pow`}{ power of the RMST based test for the values in tau_vec, holding n, k1 and k2 constant}
#'}
#' \item{simout}{a list returned if `sim = T`, with components:}
#'
#' \itemize{
#' \item{`emppowRMST`}{ empirical power of the RMST test.}
#' \item{`emppowLR`}{ empirical power of the log-rank test.}
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
                  plot = F, sim = F, M = 1000, method = 'tau_star', alpha = 0.025){
    if (is.na(power)+is.na(n)!=1) stop('One of n, power must be missing.')
    if (is.na(n)){
        RMST_truediff<-integrate(function(x) survdefT$S(x)-survdefC$S(x),
                                 lower = 0, upper = tau)$value
        RMST_trueSE <- RMST_truediff/(qnorm(1-alpha)-qnorm(1-power))
        find_n<-function(N) sqrt(evar(survdefT, k1, k2, tau, k1+k2, N)+
                                     evar(survdefC, k1, k2, tau, k1+k2, N))-RMST_trueSE
        N<- uniroot(find_n, lower = 1, upper = 10000)$root
        n <- 2*ceiling(N)
    }
    powerRMST <- powfn(survdefC, survdefT, k1, k2, tau, n, alpha)
    powerLR <- LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha)
    pKME <- RMSTeval(survdefC, survdefT, k1, k2, tau, n)
    to_ret<-list(n=n, powerRMST=powerRMST, powerLR=powerLR, pKME=pKME)
    if (plot == T){
        par(mfrow = c(2,2), mar = c(5,4,1,2))
        plot_val<-list(n_vec=NA, n_pow=NA, k1_vec=NA, k1_pow=NA,
                       k2_vec=NA, k2_pow=NA, tau_vec=NA, tau_pow = NA)
        print('Calculating power for different values of n...')
        plot_val$n_vec<-seq(from = n*.5, to = n*1.5,length.out =  20)
        plot_val$n_pow<-sapply(plot_val$n_vec, function(x) powfn(survdefC, survdefT, k1, k2, tau, x, alpha))
        print('Calculating power for different values of k1...')
        plot_val$k1_vec<-seq(from = max(0,tau-k2), to = k1*1.5,length.out =  20)
        if (k1==0) plot_val$k1_vec<-seq(from = 0, to = k2,length.out =  20)
        plot_val$k1_pow<-sapply(plot_val$k1_vec, function(x) powfn(survdefC, survdefT, x, k2, tau, n, alpha))
        print('Calculating power for different values of k2...')
        plot_val$k2_vec<-seq(from = max(0,tau-k1), to = k2*1.5,length.out =  20)
        if (k2==0) plot_val$k2_vec<-seq(from = 0, to = k1,length.out =  20)
        plot_val$k2_pow<-sapply(plot_val$k2_vec, function(x) powfn(survdefC, survdefT, k1, x, tau, n, alpha))
        print('Calculating power for different values of tau...')
        plot_val$tau_vec<-seq(from = tau*.5, to = min(k1+k2,tau*1.5),length.out =  20)
        plot_val$tau_pow<-sapply(plot_val$tau_vec, function(x) powfn(survdefC, survdefT, k1, k2, x, n, alpha))

        plot(plot_val$n_vec, plot_val$n_pow, xlab = 'Sample size', ylab = 'Power', ylim = c(0,1), type ='l')
        lines(plot_val$n_vec, sapply(plot_val$n_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,k2,x, alpha)), col = 'red')
        legend(plot_val$n_vec[12], 0.25, c('RMST','Log-rank'), col = 1:2,cex = .5, lty = 1, bty='n')
        plot(plot_val$k1_vec, plot_val$k1_pow, xlab = 'Accrual time', ylab = 'Power', ylim = c(0,1), type ='l')
        lines(plot_val$k1_vec, sapply(plot_val$k1_vec, function(x) LRpow_nonPH(survdefC, survdefT, x,k1+k2-x,n, alpha)), col = 'red')
        legend(plot_val$k1_vec[12], 0.25, c('RMST','Log-rank'), col = 1:2,cex = .5, lty = 1, bty='n')
        plot(plot_val$k2_vec, plot_val$k2_pow, xlab = 'Followup time', ylab = 'Power', ylim = c(0,1), type ='l')
        lines(plot_val$k2_vec, sapply(plot_val$k2_vec, function(x) LRpow_nonPH(survdefC, survdefT, k1,x,n, alpha)), col = 'red')
        legend(plot_val$k2_vec[12], 0.25, c('RMST','Log-rank'), col = 1:2,cex = .5, lty = 1, bty='n')
        plot(plot_val$tau_vec, plot_val$tau_pow, xlab = 'Tau', ylab = 'Power', ylim = c(0,1), type ='l')
        abline(h=LRpow_nonPH(survdefC, survdefT, k1,k2,n, alpha), col = 'red')
        legend(plot_val$tau_vec[12], 0.25, c('RMST','Log-rank'), col = 1:2,cex = .5, lty = 1, bty='n')
        to_ret[['plot_val']]<-plot_val
    }
    if (sim == T) {
        print('Simulating datasets...')
        simout<-list(emppowRMST=NA, emppowLR=NA, emppKME=NA, meandiff=NA, SDdiff=NA, meantrunc=NA, SDtrunc=NA)
        Fc_inv<-function(u){
            findroot<-function(x) 1-survdefC$S(x)-u
            return(uniroot(findroot, lower = 0, upper = (k1+k2)*1000)$root)
        }
        Ft_inv<-function(u){
            findroot<-function(x) 1-survdefT$S(x)-u
            return(uniroot(findroot, lower = 0, upper = (k1+k2)*1000)$root)
        }
        out<-data.frame(matrix(NA, nrow = M, ncol = 6))
        colnames(out)<-c('last_censT','last_censC','tau_temp','RMST_diff','RMST_reject', 'LR_reject')
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
            out$RMST_reject[m]<-out$RMST_diff[m]/(sqrt(summ[2,6]^2 + summ[1,6]^2))>qnorm(1-alpha)
            survtest<-survdiff(Surv(time, event)~grp) #, timefix = FALSE)
            out$LR_reject[m]<-((1 - pchisq(survtest$chisq, 1))<alpha &
                                   survtest$obs[1]>survtest$exp[1])
        }
        simout$emppowRMST<-sum(out$RMST_reject)/M
        simout$emppowLR<-sum(out$LR_reject)/M
        simout$emppKME<-sum(out$last_censT==F & out$last_censC==F)/M
        simout$meandiff<-mean(out$RMST_diff)
        simout$SDdiff<-sd(out$RMST_diff)
        simout$meantrunc<-mean(out$tau_temp, na.rm = T)
        simout$SDtrunc<-sd(out$tau_temp, na.rm = T)
        to_ret[['simout']]<-simout
    }
    return(to_ret)
}
