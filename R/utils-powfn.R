powfn<-function(survdefC, survdefT, k1, k2, tau, n, alpha){
    RMST_truediff<-integrate(function(x) survdefT$S(x)-survdefC$S(x),
                             lower = 0, upper = tau)$value
    RMST_trueSE<-sqrt(evar(survdefT, k1, k2, tau, k1+k2, n/2)+
                          evar(survdefC, k1, k2, tau, k1+k2, n/2))
    return(1-pnorm(qnorm(1-alpha)-RMST_truediff/RMST_trueSE))
}
