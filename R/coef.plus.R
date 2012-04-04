coef.plus <-
function(object,lam, ...){ 
  if (missing(lam)) {
    lam <- sort(object$lam.path, decreasing=TRUE)
    cat("Warning message: ")
    cat("no lam argument; object$lam.path used", "\n") 
  }
  else if (max(lam) < min(object$lam.path)) {
    lam <- sort(object$lam.path, decreasing=TRUE)
    cat("Warning message: ")
    cat("lam not reached by the plus path; x$lam.path used", "\n") 
  }
  cat("\n")
  tmp <- plus.hit.points(lam,object$lam.path)
  beta <- tmp$w1* object$beta.path[tmp$k1,] + (1-tmp$w1)* object$beta.path[tmp$k2,]
  return(beta)
}

