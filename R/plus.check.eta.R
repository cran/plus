plus.check.eta <-
function(old.eta, new.eta, etatil, xi, z, v, eps=1e-15, e) { 
### test a possible new eta and xi combination
## input
# old.eta, 1 X p: old parallelepiped indicator; new.eta: the new one to be tried 
# etatil, 2 X (|C|+1): label and new indicator for max crossing, with a dummy variable p+1
# xi: the sign of new.tau - old.tau
# z: t(x)y/n; v: 2nd derivative of penalty; eps: numerical zero; e: environment including x, y
## out put:  return(xi, s, singular.Q, invalid.cross)
# singular.Q: error for solve(Q, z[a.set])
# invalid.cross: either singular.Q or fail to cross to new.eta
 p <- dim(e$x)[2]
 a.set <- rep(TRUE,p)
 a.set[new.eta==0] <- FALSE
 ##  update e$Sigma, e$var.list and e$a.set.var
 ##  e$Sigma[,e$a.set.var] will equal to t(x)%*%x[,a.set]/n when done 
 if (e$use.Gram) e$a.set.var <- a.set # e$Sigma is complete
 else if ( sum(e$var.list) == 0 ) { # the initial case of e$Sigma = 0 
   e$var.list <- order(!a.set)[1:sum(a.set)]
   e$a.set.var <- 1:sum(a.set)
   e$Sigma[, e$a.set.var] <- e$tx %*% e$x[,a.set]
 } 
 else { # incrementally add new columns to e$Sigma if necessary 
   new.var <- a.set
   new.var[e$var.list] <- FALSE # new.var are those in a.set but nor in e$var.list 
   if (sum(new.var)>0) { # need to add 
     e$var.list <- c(e$var.list, order(!new.var)[1:sum(new.var)])
     if (length(e$var.list) > (dim(e$Sigma)[2])) { # need to make e$Sigma larger 
       tmp.Sigma <- e$Sigma
       e$Sigma <- matrix(0,p,dim(tmp.Sigma)[2]+100)
       e$Sigma[,1:(dim(tmp.Sigma)[2])] <- tmp.Sigma
     }
     e$Sigma[,(length(e$var.list) - sum(new.var)+1):length(e$var.list)] <- e$tx %*% e$x[,new.var]
   }
   e$a.set.var <- (1:length(e$var.list))[a.set[e$var.list]] # the equivalent labels of a.set
   e$a.set.var <- e$a.set.var[order(e$var.list[ e$a.set.var])]  
 }
 ## calculate for return(xi, singular.Q, s, invalid.cross), xi is given 
 if (sum(a.set)==1) Q <- e$Sigma[a.set,e$a.set.var] - v[ abs(new.eta[a.set]) ]
 else Q <- e$Sigma[a.set,e$a.set.var] - diag( v[ abs(new.eta[a.set]) ] )
 singular.Q <- FALSE
 # LINPACK=T generated an incorrect calculation, due to the use of single precision
 # s1 <- try(solve(Q, z[a.set], LINPACK = (length(Q)>100)), silent=T)
 # s1 <- try(solve(Q, z[a.set], symmetry=T), silent=T)  # option has no effect 
 s1 <- try(solve(Q, z[a.set]), silent=TRUE)
 if (inherits(s1,"try-error")==TRUE) {
   s1 <- rep(0, sum(a.set))
   singular.Q <- TRUE
 }
 s <- rep(0,p)
 if (singular.Q==FALSE) {
   s[ a.set ] <- s1
   ## check the new s
   etatil2 <- c(old.eta,0) # etatil in full length
   etatil2[ etatil[1,] ] <- etatil[2,]
   etatil2 <- etatil2[1:p] # remove dummy
   #  e$g.prime <- z - as.matrix((e$Sigma[,e$a.set.var]))%*%s1  # Eq. (2.13)
   if (length(s1)==1)  e$g.prime <- z - (e$Sigma[,e$a.set.var])*s1  # Eq. (2.13)
   else  e$g.prime <- z - (e$Sigma[,e$a.set.var])%*%s1 
   flag <- rep(TRUE,p)  # T means violation 
   flag[old.eta==etatil2] <- FALSE # noncritical j are fine
   # check according to (2.15) with eps = numerical zero
   flag[ (old.eta != new.eta) & (new.eta != 0) & (xi*(new.eta-old.eta)*s > -eps) ] <- FALSE
   flag[ (old.eta == new.eta) & (new.eta != 0) & (etatil2 != old.eta) & (xi*(etatil2-old.eta)*s < eps) ] <- FALSE 
   flag[ (old.eta != 0) & (new.eta == 0) & (xi*old.eta*e$g.prime < eps)] <- FALSE
   flag[ (old.eta == 0) & (new.eta == 0) & (etatil2 != 0) & (xi*etatil2*e$g.prime < eps)] <- FALSE
   invalid.cross <- sum(flag) > 0 # T here iff the new eta and xi pair is invalid and singular.Q==FALSE
 }
 else invalid.cross <- TRUE # invalid.cross is T if singular.Q == T
 return(list(xi=xi, singular.Q=singular.Q, s=s, invalid.cross=invalid.cross))
}

