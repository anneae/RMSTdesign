evar<-function(survdefo,  k1, k2, tau, ttime, n){
    P_enroll<-function(t,x) {
        if (t<=x) 0
        else if (x<t & t<=(k1+x))  (t-x)/k1
        else if (k1+x <t & t<= k1+k2) 1
        else if (k1+x <t & all.equal(t, k1+k2)) 1}
    toint<-function(u) sapply(u, function(x) integrate(survdefo$S, lower = x, upper = tau)$value^2*survdefo$h(x)/(n*survdefo$S(x)*P_enroll(ttime, x)))
    integrate(toint, lower = 0, upper = tau)$value
}
