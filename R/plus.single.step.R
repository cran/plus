plus.single.step <-
function(k, etatil,  z, b, tau, m, v, t, eps=1e-15, e) { 
## input 
# k:  number of parallelepiped indicators already found, including 0, k-1 steps completed 
# etatil, 2 X (|C|+1): label and new indicator for max crossing with dummy variable p+1
# NOTE: |C| > 0 required here
# z: t(x)y/n; b, tau: old turning point b at old tau; m, v, t: penalty; eps: numerical zero, 
## out put:   return(new.eta, etatil, s, tau,b, singular.Q, forced.stop, full.path) 
#   including info at the end of segment new.eta, e.g. etatil for the eta beyond new.eta 
 n <- dim(e$x)[1] # no dummy
 p <- dim(e$x)[2]
 # this is in step k, eta^{(0)},...,eta^{(k-1)} have already been found 
 sides <- length(etatil)/2-1 # number of boundaries for possible crossing
 singular.Q <- FALSE # this and next values are needed in case plus.check.eta not called 
 forced.stop <- TRUE 
 if (sides == 1) { # one-at-a-time scenario 
   new.eta <- plus.new.eta(e$eta[k,],etatil,rep(TRUE,2),p) # T for crossing the boundary and dummy 
   ## always a new eta in theory, only check when k is a multiplier of 50 to prevent numerical loop    
   if (k != round(k/50)*50) check.new.eta <- TRUE
   else check.new.eta <- plus.not.in.loop(k, new.eta, e=e) 
   if (check.new.eta==TRUE) { # do not run plus.check.eta if eta is not new, try xi = 1 first 
     tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, 1, z, v, eps=eps,e=e) 
     ## tmp holds return(xi, singular.Q, s, invalid.cross)      
     singular.Q <- tmp$singular.Q  # save if invalid cross is due to singular Q, an almost nonevent
     if ( (m>1) & ((tmp$invalid.cross) == TRUE) ) { # try xi = -1
       tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, -1, z, v, eps=eps,e=e)
       if (singular.Q == FALSE) singular.Q <- tmp$singular.Q
     }
     forced.stop <- tmp$invalid.cross # due to singular Q or not 
   } # if (check.new.eta==T)
 }  # end of the one-at-a-time case 
 else { # sides > 1, first check ties
   etatil <- etatil[, c((! e$ties[ etatil[1,1:sides] ]),TRUE) ] # remove known ties 
   sides <- length(etatil)/2-1
   if (sides > 1) { # remove new ties 
     if (e$use.Gram) Sgm.B <- e$Sigma[etatil[1,1:sides], etatil[1,1:sides]]
     else Sgm.B <- e$tx[etatil[1,1:sides],] %*% e$x[,etatil[1,1:sides]]        
     for (i in  1:(sides-1)) for (j in (i+1):sides) 
       if (Sgm.B[i,j] > 1- eps) {
         e$ties[j] <- TRUE
        }
     etatil <- etatil[, c(! e$ties[ etatil[1,1:sides] ],TRUE) ] # remove known ties 
     sides <- length(etatil)/2-1
   }
   if (sides == 1) {# go back to the one-at-time case
     new.eta <- plus.new.eta(e$eta[k,],etatil,rep(TRUE,2),p) # T for crossing the boundary and dummy 
     ## always a new eta in theory, only check when k is a multiplier of 50 to prevent numerical loop    
     if (k != round(k/50)*50) check.new.eta <- TRUE
     else check.new.eta <- plus.not.in.loop(k, new.eta, e=e) 
     if (check.new.eta==TRUE) { # do not run plus.check.eta if eta is not new, try xi = 1 first 
       tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, 1, z, v, eps=eps,e=e) 
       ## tmp holds return(xi, singular.Q, s, invalid.cross)      
       singular.Q <- tmp$singular.Q  # save if invalid cross is due to singular Q, an almost nonevent
       if ( (m>1) & ((tmp$invalid.cross) == TRUE) ) { # try xi = -1
         tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, -1, z, v, eps=eps,e=e)
         if (singular.Q == FALSE) singular.Q <- tmp$singular.Q
       }
       forced.stop <- tmp$invalid.cross # due to singular Q or not 
     } # if (check.new.eta==T)
   }  # end of the one-at-a-time case 
   else if (sides > 1) { # general case, (more than) two-at-a-time 
     M <- matrix(FALSE, 2^sides, sides + 1) # generate all possible candidates for new eta
     M[, sides+1] <- TRUE # last T for dummy
     M[ (2^(sides-1)+1):(2^sides), 1 ] <- TRUE # T for cross the indicated side
     for (i in (2:sides) ) M[,i] <- M[(2^(sides-i)+1):(3*2^(sides-i)),i-1]
     M <- M[order(apply(M,1,sum)),] # closer eta (smaller number of crossings) has higher priority  
     # initial while loop look for new eta 
     continue.while <- TRUE 
     j <- 1
     while (continue.while == TRUE) {
       new.eta <- plus.new.eta(e$eta[k,],etatil,M[j+1,],p) # skip the first M[1,]=F...FT
       new.is.new <- plus.not.in.loop(k, new.eta, e=e) 
       if (new.is.new == TRUE) { # try xi =1
         tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, 1,  z, v, eps=eps, e=e)
         ## tmp holds return(xi, singular.Q, s, invalid.cross)      
         ## invalid.cross == T iff  either (singular.Q == T) or (xi and s fail to cross to new.eta)
         if (singular.Q == FALSE) singular.Q <- tmp$singular.Q
         if ( (m>1) & ((tmp$invalid.cross) == TRUE)) { # try xi = -1
           tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, -1, z, v, eps=eps, e=e)
           if (singular.Q == FALSE) singular.Q <- tmp$singular.Q  # singular.Q should be T if happens once 
         }
         forced.stop <- tmp$invalid.cross # forced.stop is true if all "not new" or invalid
         continue.while <- tmp$invalid.cross # exit while loop for valid cross to new eta 
       } # end if (new.is.new == T)
       j <- j+1 # continue.while does not change when new.is.new == F
       if (j==2^sides) continue.while <- FALSE
     } # end while, if forced.stop == F, then tmp holds info for a valid new eta 
   } # end of (more than) two-at-a-time
   else if (sides < 1) new.eta <- rep(m+2,p) 
 }
 full.path <- FALSE
 if (forced.stop == TRUE) { # terminate the program, need value for the 9 out of 10 var
   ## dummy output for return(new.eta, etatil, s, tau,b, singular.Q, forced.stop, full.path)
   ##  values have already been assigned to new.eta, singular.Q, forced.stop, full.path
   etatil <- c(p+1,m+2)  # etatil <- as.matrix(c(p+1,m+2)) 
   s <- rep(0,p)
   tau <- 0
   b <- rep(0,p)
 }
if (forced.stop==FALSE) { # valid new eta has been found as new.eta
  ## tmp holds return(xi, singular.Q, s, invalid.cross), singular.Q = invalid.cross = F      
  xi <- tmp$xi
  s <- tmp$s
  ## need values for return(etatil, tau,b, full.path); new.eta, singular.Q, forced.stop done 
  a.set <- (new.eta!=0)
  if ( k - 50*round(k/50) == 1) # directly calculate the gradient to prevent accumulation of error 
  # e$grad <- tau*z - as.matrix(e$Sigma[,e$a.set.var])%*%b[a.set]    
  if (sum(a.set)<2) 
     e$grad <- tau*z - (e$Sigma[,e$a.set.var])*b[a.set]          # gradient
  else
     e$grad <- tau*z - (e$Sigma[,e$a.set.var])%*%b[a.set]     
  ## compute DeltaJ and new etatil   
  DeltaJ <- rep(-1,p) # -1 represents infinity     
  new.etatil <- rep(m+2,p) # m+2 is the dummy value for etatil calculation 
  ind <- (xi*s>=eps) & (new.eta != 0) & (new.eta < m)
  if (sum(ind)>0) { # move towards higher eta value
    DeltaJ[ind] <- xi*(e$t.fun[new.eta[ind]+1+m]-b[ind])*plus.recip(s[ind],eps=eps)
    new.etatil[ind] <- new.eta[ind] + 1
  }
  ind <- (xi*s<= -eps) & (new.eta != 0) & (new.eta > -m)
  if (sum(ind)>0) { # move towards lower eta value
    DeltaJ[ind] <- xi*(e$t.fun[new.eta[ind]+m]-b[ind])*plus.recip(s[ind],eps=eps)
    new.etatil[ind] <- new.eta[ind] - 1
  }
  ind <- (xi*e$g.prime>= eps) & (new.eta == 0)
  if (sum(ind)>0) { # new variable to with positive b
    DeltaJ[ind] <- xi*(1-e$grad[ind])*plus.recip(e$g.prime[ind],eps=eps)
    new.etatil[ind] <- 1
  }   
  ind <- (xi*e$g.prime<= -eps) & (new.eta == 0)
  if (sum(ind)>0) { # new variable to with negative b 
    DeltaJ[ind] <- xi*(-1-e$grad[ind])*plus.recip(e$g.prime[ind],eps=eps)
    new.etatil[ind] <- -1
  }
  # compute Delta, new.etatil, full.path
  old.etatil <- rep(m+2,p+1) # again m+2 is the dummy value
  old.etatil[etatil[1,]] <- etatil[2,]
  old.etatil <- old.etatil[1:p]  # old etatil in long format 
  if ( all(DeltaJ == (-1)) ) { # all infinity
    Delta <- -1                              # lse has attained for m > 1 and will be the next for lasso  
    full.path <- TRUE
    index.new.etatil <- rep(FALSE,p) # no new sides, all done 
  }
  else { # some stopping time is finite
    Delta <- min(DeltaJ[DeltaJ != (-1)])
    # index.new.etatil <- (DeltaJ < (Delta + eps)) & (DeltaJ !=(-1)) # does not work well 
    index.new.etatil <- DeltaJ==Delta # sides hit, numerically delicate 
    ind <- (abs(s) < eps) & (new.eta != 0) & (old.etatil !=m+2)  
    if (sum(ind)>0) { # sides not moving 
      new.etatil[ind] <- old.etatil[ind]
      index.new.etatil[ind] <- TRUE
    }
    ind <- (abs(e$g.prime) < eps) & (new.eta == 0) & (old.etatil !=m+2) 
    if (sum(ind)>0) { # sides not moving  
      new.etatil[ind] <- old.etatil[ind]
      index.new.etatil[ind] <- TRUE
    }
  } # end of else 
  # calculate the new etatil
  index.new.etatil <- c(index.new.etatil,TRUE) # add dummy
  e$etatil2[2,1:p] <- new.etatil  # other entries of e$etatil2 never change value 
  etatil <- e$etatil2[,index.new.etatil] 
  # etatil is m+2 at dummy column p+1 when full.path == T 
  # still need to compute return(tau,b) 
  if (Delta != -1) {
    tau <- tau + xi*Delta
    b <- b + xi*Delta*s  # update b, b[k+1,] <- b[k,] + ...
    e$grad <- e$grad + xi*Delta*e$g.prime 
    if (tau * eps > 1/5) {  # end at perfect fit with zero gradient numerically
      full.path <- TRUE
      index.new.etatil <- rep(FALSE,p) # 
    } 
    if (tau < e$tau1 - eps) { # shall not return to the beginning
      forced.stop <- TRUE
    }
  }
  else {
    tau <- xi
    b <- s
  }
} # end if forced.stop == F 
return(list(new.eta=new.eta, etatil=etatil, s=s, tau=tau, b=b, singular.Q=singular.Q, 
                   forced.stop = forced.stop, full.path = full.path)) 
}

