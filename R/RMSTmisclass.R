#' Estimate Treatment and Biomarker Effects on RMST in a Biomarker Stratified Clinical Trial with Misclassification
#'
#' @param dat the dataset that will be used to calculate estimates. The dataset needs to include columns called
#' `dtime` and `dstatus` where `dtime` is time to event/censoring and `dstatus` is 1 if the participant
#' had an event, 0 if they were censored. It also needs to include columns called `trt` (1 for treatment, 0 for control)
#' and `test` (0 for those who tested negative for the biomarker, 1 for those who tested positive).
#' @param sens the sensitivity of the biomarker test.
#' @param spec the specificity of the biomarker test.
#' @param prev the prevalence of the biomarker.
#' @param tau restriction time for the RMST.
#'
#' @return a dataframe with estimated effects on RMST up to time tau, including
#' treatment effects in the marker positive (M+) and negative (M-) groups,
#' marker effects in the treatment (T+) and control (T-) groups,
#' and the interaction effect (treatment effect in M+ group minus treatment effect in M- group).
#' The dataframe also contains standard errors and p-values testing the null hypothesis
#' that the difference in RMST = 0.
#' @export
#'
#' @examples
#' set.seed(23)
#' dat<- data.frame(trt = rep(0:1, each = 50), marker = rbinom(100, 1, prob = .3), cens = runif(100, 5, 15))
#' dat$ttime <- rexp(100, rate = exp(.1*dat$trt + .1*dat$marker - 0.7*dat$trt*dat$marker))
#' dat$dtime <- pmin(dat$ttime, dat$cens)
#' dat$dstatus <- (dat$ttime < dat$cens)
#' dat$test <- rbinom(100, 1, ifelse(dat$marker == 1, 0.8, 0.2))
#' RMSTmisclass(dat, sens=.8, spec=.8, prev=.3, tau = 10)

RMSTmisclass <- function(dat, sens, spec, prev, tau){
    ppv<-prev*sens/(prev*sens+(1-prev)*(1-spec))
    npv<-(1-prev)*spec/(prev*(1-sens)+(1-prev)*spec)

    km <- survfit(Surv(dtime, dstatus)~trt, data= dat[dat$test == 1, ])
    tab <- summary(km, rmean = tau)$table
    naive_m1_t0<- tab[1,5] # estimated RMST in control group
    naive_m1_t1<- tab[2,5] # estimated RMST in trt group
    naive_m1_t0_se<- tab[1,6]
    naive_m1_t1_se<- tab[2,6]

    km <- survfit(Surv(dtime, dstatus)~trt, data= dat[dat$test == 0, ])
    tab <- summary(km, rmean = tau)$table
    naive_m0_t0<- tab[1,5] # estimated RMST in control group
    naive_m0_t1<- tab[2,5] # estimated RMST in trt group
    naive_m0_t0_se<- tab[1,6]
    naive_m0_t1_se<- tab[2,6]

    adj_m1_t0 <- (npv*(naive_m1_t0)-(1-ppv)*naive_m0_t0)/(npv*ppv-(1-npv)*(1-ppv))
    adj_m1_t1 <- (npv*(naive_m1_t1)-(1-ppv)*naive_m0_t1)/(npv*ppv-(1-npv)*(1-ppv))
    adj_m0_t0 <- (ppv*(naive_m0_t0)-(1-npv)*naive_m1_t0)/(npv*ppv-(1-npv)*(1-ppv))
    adj_m0_t1 <- (ppv*(naive_m0_t1)-(1-npv)*naive_m1_t1)/(npv*ppv-(1-npv)*(1-ppv))

    adj_m1_t0_se <- sqrt(npv^2*(naive_m1_t0_se)^2+(1-ppv)^2*naive_m0_t0_se^2)/(npv*ppv-(1-npv)*(1-ppv))
    adj_m1_t1_se <- sqrt(npv^2*(naive_m1_t1_se)^2+(1-ppv)^2*naive_m0_t1_se^2)/(npv*ppv-(1-npv)*(1-ppv))
    adj_m0_t0_se <- sqrt(ppv^2*(naive_m0_t0_se)^2+(1-npv)^2*naive_m1_t0_se^2)/(npv*ppv-(1-npv)*(1-ppv))
    adj_m0_t1_se <- sqrt(ppv^2*(naive_m0_t1_se)^2+(1-npv)^2*naive_m1_t1_se^2)/(npv*ppv-(1-npv)*(1-ppv))

    res <- data.frame(est = c(adj_m1_t1-adj_m1_t0, adj_m0_t1-adj_m0_t0, adj_m1_t1-adj_m0_t1, adj_m1_t0-adj_m0_t0, (adj_m1_t1-adj_m1_t0)-(adj_m0_t1-adj_m0_t0)),
                      se = c(sqrt(adj_m1_t1_se^2+adj_m1_t0_se^2), sqrt(adj_m0_t1_se^2+adj_m0_t0_se^2),
                             (sqrt(naive_m1_t1_se^2+naive_m0_t1_se^2)/(npv*ppv-(1-npv)*(1-ppv))),(sqrt(naive_m1_t0_se^2+naive_m0_t0_se^2)/(npv*ppv-(1-npv)*(1-ppv))),
                             (sqrt(naive_m1_t1_se^2+naive_m1_t0_se^2+naive_m0_t1_se^2+naive_m0_t0_se^2)/(npv*ppv-(1-npv)*(1-ppv)))))
    res$pval <-  2*pnorm(abs(res$est)/res$se,lower.tail = F)

    colnames(res) <- c('Estimate', 'Standard error','P-value')
    rownames(res) <- c('Treatment in M+','Treatment in M-','Marker in T+','Marker in T-','Interaction')

    return(res)
}
