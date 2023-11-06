### Define a survdef object for when a group is a mix of two groups

#### Allows each argument to be a survdef object, so true distributions in each
#### group can be piecewise constant hazard or Weibull!
#### If the curves are both exponential, add "exponential = T" which makes the
#### code run faster.

library(numDeriv)

survdefMix <- function(survdef1, survdef2, w1, w2, exponential = F){
    if (w1+w2 !=1) stop('w1 and w2 should sum to 1.')
    if (exponential == T){
        h1<-survdef1$h(1)
        h2<-survdef2$h(1)
        S <- function(x) sapply(x, function(y) w1*exp(-h1*y) + w2*exp(-h2*y))
        h <- function(x) sapply(x, function(y) (h1*w1*exp(-h1*y) + h2*w2*exp(-h2*y))/S(y))
    }
    else{
        S <- function(x) sapply(x, function(y) w1*survdef1$S(y) + w2*survdef2$S(y))
        h <- function(x) sapply(x, function(y) -(w1*grad(survdef1$S, y) + w2*grad(survdef2$S, y))/S(y))
    }
    return(list(S = S, h = h))
}
