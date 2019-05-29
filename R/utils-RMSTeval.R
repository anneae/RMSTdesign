RMSTeval<-function(survdefC, survdefT, k1, k2, tau, n){
    if (tau<=k2) return(1)
    P_x_followup<-function(x){ # density for exactly x amount of followup
        if (x<k2) 0
        else if (x<=(k1+k2)) 1/(k1)
        else 0
    }
    CDF_x_followup<-function(x){ # CDF for x amount of followup
        if (x<k2) 0
        else if (x<=(k1+k2)) (x-k2)/(k1)
        else 1
    }
    # control group
    toint<-function(u) sapply(u, function(t) (1-(1-CDF_x_followup(t))*survdefC$S(t))^(n/2-1)*P_x_followup(t)*survdefC$S(t))
    peval_c<-1-(n/2)*integrate(toint, lower = k2, upper = min(k1+k2, tau))$value
    # trt group
    toint<-function(u) sapply(u, function(t) (1-(1-CDF_x_followup(t))*survdefT$S(t))^(n/2-1)*P_x_followup(t)*survdefT$S(t))
    peval_t<-1-(n/2)*integrate(toint, lower = k2, upper = min(k1+k2, tau))$value
    return(peval_c*peval_t)
}
