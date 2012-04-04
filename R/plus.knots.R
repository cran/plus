plus.knots <-
function(t=0,m=0) {
### compute the function t() as in (2.6) without infinity and - infinity
### note that t_1=t(1)=t(0) in (2.6)
if (m==0) { # default MC penalty
  m <- 2 
  t <- c(0,3)
}  
t2 <- rep(0,2*m)
t2[1:m] <- - t[m:1]
t2[(m+1):(2*m)] <- t[1:m]
return(c(t2,0))
}

