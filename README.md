
<!-- README.md is generated from README.Rmd. Please edit that file -->
RMSTdesign
==========

<!-- badges: start -->
<!-- badges: end -->
The goal of `RMSTdesign` is to make it easy to design clinical trials with the difference in restricted mean survival time as the primary endpoint. This includes calcuating asymptotic power or sample size under user-specified design parameters and performing simulations to determine the operating characteristics of the test of restricted mean survival time.

Installation
------------

You can install `RMSTdesign` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("anneae/RMSTdesign")
```

Documentation
-------------

Documentation and vignettes are available [here](https://anneae.github.io/RMSTdesign/).

Example
-------

In this basic example, we specify that the control group is expected to have a median survival of 3 years, the treatment group is expected to have a median survival of 4 years, and we plan to accrue patients over one year and follow them until the last patient has 3 years of followup. Also, we select *τ* = 3 years as a clinically meaningful restriction time. The code below determines the required sample size for a trial with 80% power.

``` r
library(RMSTdesign)
con<-survdef(times = 3, surv = 0.5) 
trt<-survdef(times = 4, surv = 0.5)

RMSTpow(con, trt, k1 = 1, k2 = 3, tau = 3, power = 0.8)
#> $n
#> [1] 1024
#> 
#> $powerRMST
#> [1] 0.8004313
#> 
#> $powerRMSTToverC
#> [1] 0.8004313
#> 
#> $powerRMSTCoverT
#> [1] NA
#> 
#> $powerLRToverC
#> [1] 0.9033999
#> 
#> $powerLRCoverT
#> [1] NA
#> 
#> $powerLRtauToverC
#> [1] 0.8707378
#> 
#> $powerLRtauCoverT
#> [1] NA
#> 
#> $pKME
#> [1] 1
```
