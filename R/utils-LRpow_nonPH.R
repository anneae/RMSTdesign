LRpow_nonPH<-function(survdefC, survdefT, k1, k2, n, alpha){
    myk1<-round(k1, 3)
    myk2<-round(k2, 3)
    L<-floor((myk1+myk2)*1000)
    xi<-survdefT$h(seq(from =1, to = myk1+myk2, length.out = L))/
        survdefC$h(seq(from =1, to = myk1+myk2, length.out = L))
    h_scl<-survdefC$h(seq(from =1, to = myk1+myk2, length.out = L))/1000
    nt<-nc<-rep(NA, L)
    nt[1]<-nc[1]<-n/2
    for (i in 2:L){
        if ((i-1)<myk2*1000){
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
    phi<-sum(D*(xi*p/(1+xi*p)-p/(1+p)))/sqrt(sum(D*p/(1+p)^2))
    return(pnorm(-qnorm(1-alpha)-phi))
}
