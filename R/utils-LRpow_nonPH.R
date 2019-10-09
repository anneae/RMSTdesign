LRpow_nonPH<-function(survdefC, survdefT, k1, k2, n, alpha, tau = NA){
    myk1<-round(k1, 3)
    myk2<-round(k2, 3)
    tot<-myk1+myk2
##    L<-floor((tot)*5000) # L is integer due to round above
    L<-10000
    tis<-seq(from =0, to = tot, length.out = L)
    # Assess the hazard ratio and scaled hazard at the middle of each interval.
#    xi<-survdefT$h(seq(from =0, to = tot, length.out = L)+1/10000)/
#        survdefC$h(seq(from =0, to = tot, length.out = L)+1/10000)
    xi<-survdefT$h(tis+tot/(2*(L-1)))/survdefC$h(tis+tot/(2*(L-1)))
    h_scl<-survdefC$h(tis+tot/(2*(L-1)))*tot/(L-1)
#    h_scl<-survdefC$h(seq(from =0, to = tot, length.out = L)+1/10000)/5000
    nt<-nc<-rep(NA, L)
    nt[1]<-nc[1]<-n/2
    for (i in 2:L){
#        if ((i-1)<myk2*5000){
        if (tis[i-1]<myk2){
            nc[i]<-nc[i-1]*(1-h_scl[i-1])
            nt[i]<-nt[i-1]*(1-h_scl[i-1]*xi[i-1])
        }
        else {
            nc[i]<-nc[i-1]*(1-h_scl[i-1]-(1/(L-(i-1))))
            nt[i]<-nt[i-1]*(1-h_scl[i-1]*xi[i-1]-(1/(L-(i-1))))
        }
    }
    nc[L]<-nt[L]<-0
    p<-nt/nc
    p[L]<-0
    D<-nc*h_scl + nt*h_scl*xi
    if (is.na(tau)) phi<-sum(D*(xi*p/(1+xi*p)-p/(1+p)))/sqrt(sum(D*p/(1+p)^2))
    else {
        if (tau<=0) stop('Tau must be greater than zero.')
#        mytau <- floor(round(tau, 3)*5000)
#        mytau<-min(mytau, L)
#        phi<-sum((D*(xi*p/(1+xi*p)-p/(1+p)))[1:mytau])/sqrt(sum((D*p/(1+p)^2)[1:mytau]))
        phi<-sum((D*(xi*p/(1+xi*p)-p/(1+p)))[tis<=tau])/sqrt(sum((D*p/(1+p)^2)[tis<=tau]))
    }
#    print(pnorm(-qnorm(1-alpha)-phi)) # power for superiority of trt
#    print(pnorm(qnorm(1-alpha)-phi, lower.tail = F)) # power for inferiority
    return(pnorm(-qnorm(1-alpha)-phi))
}
