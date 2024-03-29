% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMSTpow.R
\name{RMSTpow}
\alias{RMSTpow}
\title{Sample Size and Power for the Test of the Difference in Restricted Mean Survival Time}
\usage{
RMSTpow(
  survdefC,
  survdefT,
  k1,
  k2,
  tau,
  n = NA,
  power = NA,
  plot = F,
  sim = F,
  M = 1000,
  method = "tau_star",
  alpha = NA,
  two.sided = F
)
}
\arguments{
\item{survdefC}{the survival distribution of the control group, as a list in the form output by \code{survdef}.}

\item{survdefT}{the survival distribution of the treatment group, as a list in the form output by \code{survdef}.}

\item{k1}{length of the accrual period. We assume subjects will accrue
uniformly over the interval \verb{(0, k1)} and then be followed until trial time \code{k1+k2}.}

\item{k2}{length of the follow-up period.}

\item{tau}{restriction time for the RMST.}

\item{n}{total sample size for both groups. 1:1 randomization is assumed.
Either \code{n} or \code{power} can be specified, and the other value will be calculated.}

\item{power}{the desired power.}

\item{plot}{if T, plots of the assumed survival distributions and power
as a function of sample size, accrual time ka1and follow-up time k2 will be produced.
Default is F.
The power of the RMST test, the log-rank test using all available followup and
the log-rank test using only followup to time tau after randomization will be
displayed. If two-sided=T, the power of the test for superiorty (treatment over
control) and inferiority (control over treatment), are represented with solid
and dashed lines, respectively.}

\item{sim}{if T, simulations will be conducted and empirical power and other
summary statistics will be provided. Default is F. The hypothesis tests are
carried out based on the normal approximation with the variance estimated
according to the Greenwood plug-in/infinitesimal jackknife method.
Specifying situations where survival doesn't go to zero in a reasonable
amount of time (trial length times 1000) will lead to problems if the
\code{sim=T} option is used.}

\item{M}{number of simualations. Default is 1000.}

\item{method}{modification to be used in simulations if the Kaplan-Meier
estimate is not defined at time \code{tau} in either group. Default is 'tau_star',
which changes the restriction time to the last censoring time,
if the last observation is censored at a time earlier than \code{tau}.
Other possible values are 'gill', 'efron', 'tau_star', 'risk1', 'risk2',
and 'risk5'. The riskX' options indicate estimating RMST difference at the
latest time at which at least X people are at risk in each group,
irrespective of the value of \code{tau}.}

\item{alpha}{type I error level. Default is 0.025 if 'two.sided'=F and 0.05
if 'two.sided'=T.}

\item{two.sided}{whether a two-sided test is desired. Default is F, meaning that
all reported power values correspond to a test of the superiority of treatment
over control. If set to T, the power for a test of superiority (treatment over control)
and inferiority (control over treatment) will be reported separately in the results;
the power of a two-sided test is the sum of two.}
}
\value{
a list with components
\item{n}{the user-specified n, or if n was left blank, the n needed to achieve the user-specified power.}
\item{powerRMST}{the user-specified power, or if power was left blank, the asymptotic power of the RMST test.
If \code{one-sided=T}, \code{powerRMST} is equivalent to \code{powerRMSToverC}.
If \code{one-sided=F}, \code{powerRMST} is equivalent to the sum of the power of a one-sided test in each direction, i.e.
\code{powerRMSToverC + powerRMSCoverT}.}
\item{powerRMSToverC}{the asymptotic power for a test of superiority of treatment over control.}
\item{powerRMSCoverT}{the asymptotic power for a test of superiority of control over treatment.
If a one-sided test is specified, this is set to NA.}
\item{powerLRToverC}{the asymptotic power of the log-rank test of superiority
of treatment over control.}
\item{powerLRCoverT}{the asymptotic power of the log-rank test of superiority
of control over treatment. If a one-sided test is specified, this is set to NA.}
\item{powerLRtauToverC}{the asymptotic power of the log-rank test of superiority
of treatment over control, using only data up to time tau after randomization.}
\item{powerLRtauCoverT}{the asymptotic power of the log-rank test of superiority
of control over treatment, using only data up to time tau after randomization.
If a one-sided test is specified, this is set to NA.}
\item{pKME}{the probability that you will be able to estimate RMST difference
at time tau using the standard Kaplan-Meier estimator. If the last observation
in either group is censored, and the censoring time is less than tau, the
Kaplan-Meier estimate is not defined through time tau, and the RMST difference
cannot be estimated using the standard area under the Kaplan-Meier curve. A
modified estimator must be used.}
\item{simout}{a list returned if \code{sim = T}, with components:}

\itemize{
\item{\code{emppowRMSTToverC}}{ empirical power of the RMST test for the superiority of
treatment over control.}
\item{\code{emppowRMSTCoverT}}{ empirical power of the RMST test for the superiority of
control over treatment. If a one-sided test is specified, this is set to NA.}
\item{\code{emppowLRToverC}}{ empirical power of the log-rank test for the superiority of
treatment over control.}
\item{\code{emppowLRCoverT}}{ empirical power of the log-rank testthe superiority of
control over treatment. If a one-sided test is specified, this is set to NA.}
\item{\code{emppowLRtauToverC}}{ empirical power of the log-rank test for the superiority of
treatment over control, using only data up to time tau after randomization.}
\item{\code{emppowLRtauCoverT}}{ empirical power of the log-rank testthe superiority of
control over treatment, using only data up to time tau after randomization.
If a one-sided test is specified, this is set to NA.}
\item{\code{emppKME}}{ proportion of simulations where the standard KM estimator was used.}
\item{\code{meandiff}}{ mean estimated difference in RMST across the simulated datasets.}
\item{\code{SDdiff}}{ standard deviation of the estimated difference in RMST across the simulated datasets.}
\item{\code{meantrunc}}{ mean truncation time used in the simulated datasets (may be smaller than tau if method = 'tau_star' or 'riskX' options are used).}
\item{\code{SDtrunc}}{ standard deviation of the truncation time used in the simulated datasets.}
}
}
\description{
Determine the asymptotic power of the test of RMST under a given trial design, or
calculate the samples size needed to achieve a desired power.
}
\examples{
con<-survdef(times = 3, surv = 0.5)
trt<-survdef(haz = 0.67*con$h(1))
RMSTpow(con, trt, k1 = 0, k2 = 3, tau = 3, power = 0.8)
RMSTpow(con, trt, k1 = 0, k2 = 3, tau = 3, n = 552)
}
