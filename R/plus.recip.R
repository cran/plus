plus.recip <-
function(Q, eps=1e-15){
### compute reciprocal of reals    
Q[ (Q>=0) & (Q <eps) ] <- eps
Q[ (Q<0) & (Q> - eps) ] <- - eps
return(1/Q)
}

