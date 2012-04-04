plus <-
function(x,y, method = c("lasso", "mc+", "scad", "general"), m=2, gamma,v,t, 
  monitor=FALSE, normalize = TRUE, intercept = TRUE, 
  Gram, use.Gram = FALSE, eps=1e-15, max.steps=500, lam)
{ 
 ## fill missing values
 if (missing(gamma)) gamma <- 0
 if (missing(v)) v <- 0
 if (missing(t)) t <- 0
 if (missing(Gram)) Gram <- 0
 if (missing(lam)) lam <- -1
 lam.min <- min(lam)
 ## check data
 data.dim.error <- FALSE
 if (! is.matrix(x)) data.dim.error <- TRUE
 else {
   n <- dim(x)[1]
   if ( (n<=1) | (n!=length(y)) ) data.dim.error <- TRUE
   p <- dim(x)[2]
   if (p<= 1) data.dim.error <- TRUE
 }
 if (data.dim.error == TRUE){ 
   cat("Data dimention error.", "\n")
   return(list(data.dim.error = TRUE))
 }
 ## standardize x 
 orig.x <- x
 if (intercept)  x <- t(t(x)-apply(x,2,mean)) # center to mean zero
 if (normalize) {
   normalize.factor <- sqrt(apply(x^2,2,mean))
   x <- t(t(x)/normalize.factor)
 }
 ## compute Gram matrix if necessary
 ## if ((!use.Gram) & (p<=200)) use.Gram <- TRUE
 if ((use.Gram) & (p > 1000)) use.Gram <- FALSE
 if (use.Gram)
   if (! is.matrix(Gram)) Gram <- t(x)%*%x
   else if ((dim(Gram)[1]!=p) || (dim(Gram)[2]!=p)) Gram <- t(x)%*%x
 ## determine method and compute penalty        
 if (length(method)==1)
   if (method == "lasso") m <- 1
   else if (method == "mc+") m <- 2
   else if (method == "scad") m <- 3
 if (m==1) method <- "LASSO"
 else if (m==2) method <- "MC+"
 else if (m==3) method <- "SCAD"
 else method <- "PLUS" 
 if (monitor == TRUE) {
   cat("\n")
   cat("Sequence of", method, "moves:", "\n","\n")
 }
 if (m<4) { # provide default penalty for lasso, MCP and SCAD 
   if ((m>1) & (use.Gram) & (gamma == 0)) {
     max.corr <- max(abs(Gram[ row(Gram) != col(Gram) ]))/n
     if ((1 - max.corr) >= (1e-5) )   gamma <- 2/(1- max.corr)
   }
   if ((m==3) & (gamma <=1)) gamma <- 0  
   tmp <- plus.penalty(m,a=gamma)
   m <- tmp$m
   gamma <- tmp$gamma
   v <- tmp$v
   t <- tmp$t
 }
 ## declare variables 
 names <- dimnames(x)
 b <- matrix(0, max.steps+1,p)
 etatil <- matrix(0,2,p+1)
 tau <- rep(0, max.steps+1)
 my.env <-environment()
 my.env$x <-x
 my.env$y <- y
 my.env$tx <- t(x)/n
 my.env$eta <- b
 my.env$g.prime <- rep(0,p)
 my.env$grad <- rep(0,p)
 my.env$etatil2 <- etatil
 my.env$etatil2[1,] <- 1:(p+1)
 if (use.Gram) my.env$Sigma <- Gram/n
 else { 
   my.env$Sigma <- matrix(0,p,100)  # e$a.set.var will match a.set
   my.env$var.list <- 0
 }
 my.env$use.Gram <- use.Gram
 my.env$ties <- rep(FALSE,p)
 my.env$t.fun <- plus.knots(t,m)
 z <- my.env$tx %*%y 
 my.env$tau1 <- 1/max(abs(z))
 ## initialization: for the initial step k=0, tau^{(0)} is needed, b[1,] <- 0 done already
 tau[1] <- my.env$tau1   
 etatil[1,] <- 1:(p+1)   # the first etatil
 etatil[2,1:p] <- sign(z)
 etatil <- etatil[ , c(abs(z) >= max(abs(z)) - eps,TRUE) ]
 k <-1 
 exit.while <- FALSE
 last.valid.s <- rep(0,p)
 ## the main loop
 while (exit.while==FALSE) { # segment k: model eta[k+1,] begins with b[k,] and ends with b[k+1,]
   tmp <- plus.single.step(k, etatil, z, b[k,], tau[k], m, v, t, eps=eps,e=my.env)
   ## tmp holds return(new.eta, etatil, s, tau, b, singular.Q, forced.stop, full.path) 
   etatil <- tmp$etatil     # etatil = 0 if forced.stop, length(etatil)=2 if Ise.termination
   exit.while <- tmp$forced.stop # cannot find a valid new.eta
   full.path <- tmp$full.path 
   if (exit.while == FALSE) { # save data from tmp if not forced to stop 
     eta[k+1,] <- tmp$new.eta
     if (monitor == TRUE) { # print info to monitor the iterations
       var.change <- sign(abs(sign(eta[k+1,]))-abs(sign(eta[k,])))
       if (sum(abs(var.change))==1) {
         if(sum(var.change)==1)
           action <-"added "
         else
           action <-"removed " 
           var_name <- (order(-abs(var.change))[1])
         if (!is.null(names))
      var_name <- names[[2]][var_name]      
         cat("Step",k,": ",action,var_name,"\n")
       }  
     }
     last.valid.s <- tmp$s  # save the s for the last valid eta 
     tau[k+1]<- tmp$tau
     b[k+1,] <- tmp$b  # b[k+1] is at the end of segment eta[k+1,] 
     b[k+1,eta[k+1,]==0] <- 0 # more accurate b=0, b[k,eta[k+1,]==0] <- 0 is also ok
     k <- k+1
     if ((k > max.steps)|(full.path)|(tau[k]*lam.min > 1)) { # stop with a valid eta
       exit.while <- TRUE
       if ((full.path) & (m>1)) k <- k-1
     } 
   } # if (exit.while == FALSE) 
 } # end the main loop
 ## tmp holds return(new.eta, etatil, s, tau, b, singular.Q, forced.stop, full.path) 
 singular.Q <- tmp$singular.Q
 forced.stop <- tmp$forced.stop
 lam.path <- 1/tau[1:k]
 beta.path <- lam.path * b[1:k,]
 if (normalize) {
    beta.path <- t(t(beta.path)/normalize.factor)
 }
 lam.path[k] <- max(abs(my.env$tx %*%(y-x%*%beta.path[k,])))
 if (full.path) lam.path[k] <- 0
 eta <- eta[1:k,]
 total.steps <- k-1
 ## calculate the estimator beta at the given lam or the set of values of lam.path
 x <- orig.x
 if (sum(lam) == -1) lam <- sort(lam.path, decreasing=TRUE)
 tmp <- plus.hit.points(lam,lam.path) 
 total.hits <- tmp$total.hits
 if (total.hits > 0) {
   beta <- tmp$w1*beta.path[tmp$k1,] + (1-tmp$w1)*beta.path[tmp$k2,]
   lam <- lam[1:total.hits]
   dim <- apply(abs(beta)>0, 1, sum)
   if (total.hits == 1)
     r.square <- 1- sum((y - x%*%beta)^2)/sum(y^2)
   else
     r.square <- 1- apply((x%*%t(beta) - y)^2, 2, sum)/sum(y^2)   
     obj<-list(x=x,y=y,eta=eta, beta.path=beta.path, lam.path=lam.path, 
    beta=beta, lam=lam, dim=dim, r.square = r.square, total.hits=total.hits, 
  method = method, gamma = gamma, total.steps=total.steps, max.steps=max.steps, 
  full.path=full.path, forced.stop=forced.stop, singular.Q=singular.Q)
 }
 else 
   obj<-list(x=x,y=y,eta=eta, beta.path=beta.path, lam.path=lam.path, total.hits=total.hits, 
  method = method, gamma = gamma, total.steps=total.steps, max.steps=max.steps, 
  full.path=full.path, forced.stop=forced.stop, singular.Q=singular.Q)  
 class(obj) <- "plus"
 return(obj)  
}

