predict.plus <-
function(object,lam,newx, ...) { 
  if (missing(newx)) {
    flag <- FALSE 
    cat("Warning message: ")
    cat("no newx argument; newy not produced", "\n") 
  } 
  else flag <- TRUE
  if (missing(lam)) {
    lam <- sort(object$lam.path, decreasing=TRUE)
    cat("Warning message: ")
    cat("no lam argument; object$lam.path used", "\n") 
  }
  else if (max(lam) < min(object$lam.path)) {
    lam <- sort(object$lam.path, decreasing=TRUE)
    cat("Warning message: ")
    cat("lam not reached by the plus path; object$lam.path used", "\n") 
  }
  cat("\n")
  tmp <- plus.hit.points(lam,object$lam.path)
  total.hits <- tmp$total.hits
  beta <- tmp$w1* object$beta.path[tmp$k1,] + (1-tmp$w1)* object$beta.path[tmp$k2,]
  lam <- lam[1:total.hits]
  dim <- apply(abs(beta)>0, 1, sum)
  if (total.hits == 1) {
    r.square <- 1- sum((object$y - object$x%*%beta)^2)/sum(object$y^2)
    if (flag) {
      if (is.matrix(newx)) newy <- newx%*% beta
      else newy <- sum(beta * newx)
    }
  }
  else {
    r.square <- 1- apply((object$x%*%t(beta) - object$y)^2, 2, sum)/sum(object$y^2)
    if (flag) {
      if (is.matrix(newx)) newy <- beta %*% t(newx)
      else newy <- beta %*% newx
    }
  }
  if (flag) 
    return(list(lambda=lam, coefficients=beta, dimension=dim, r.square = r.square, step = tmp$k2-1,  
  method = object$method, newy = newy))
  else
    return(list(lambda=lam, coefficients=beta, dimension=dim, r.square = r.square, step = tmp$k2-1,  
  method = object$method))
}

