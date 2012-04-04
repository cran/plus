plot.plus <-
function(x, xvar=c("lam","step"), yvar=c("coef","newy","lam","dim","R-sq"), 
        newx, step.interval, lam.interval, predictors, ...) {
  # determin xvar in the plot
  if (length(unique(c(xvar,"lam")))==length(unique(c(xvar)))) xvar <- "lambda"
  else if (length(unique(c(xvar,"step")))==length(unique(c(xvar)))) xvar <- "step"
  else xvar <- "lambda"
  # determin yvar in the plot
  if (length(unique(c(yvar,"coef")))==length(unique(c(yvar)))) yvar <- "coefficients"
  else if (length(unique(c(yvar,"newy")))==length(unique(c(yvar)))) yvar <- "newy"
  else if (length(unique(c(yvar,"lam")))==length(unique(c(yvar)))) yvar <- "lambda"
  else if (length(unique(c(yvar,"dim")))==length(unique(c(yvar)))) yvar <- "dimension"
  else if (length(unique(c(yvar,"R-sq")))==length(unique(c(yvar)))) yvar <- "R-square"
  else yvar <- "coefficients"
  if (yvar != "newy") newx <- matrix(0,2,dim(x$x)[2])
  if ((yvar == "newy") & missing(newx)) {
    cat("Warning message: ")
    cat("no newx argument", "\n") 
    return()
  }
  if (missing(predictors)) predictors <- 1:(dim(x$x)[2])
  ## plot 
  if (xvar == "step") {
    xtmp <-  0:(x$total.steps)
    plot.set <- rep(TRUE, length(xtmp))
    if (! missing(step.interval)) {
      plot.set[xtmp > max(step.interval)] <- FALSE
      plot.set[xtmp < min(step.interval)-1] <- FALSE
      xtmp <- xtmp[plot.set]
    }
    if (yvar == "coefficients") ytmp <- x$beta.path[plot.set, ]
    if ((yvar == "newy") & is.matrix(newx)) ytmp <- x$beta.path[plot.set, ]%*%t(newx)
    if ((yvar == "newy") & (! is.matrix(newx))) ytmp <- x$beta.path[plot.set, ]%*%newx
    if (yvar == "lambda") ytmp <- x$lam.path[plot.set]
    if (yvar == "dimension") ytmp <- apply(abs(sign(x$eta[plot.set,])),1, sum)
    if (yvar == "R-square") 
      ytmp <- 1- apply((x$x%*%t(x$beta.path[plot.set, ]) - x$y)^2, 2, sum)/sum(x$y^2)
  }    
  if (xvar == "lambda") {
    xtmp <- sort(x$lam.path,decreasing=TRUE)
    tmp <- predict(x,xtmp,newx)
    plot.set <- rep(TRUE, length(xtmp))
    if (! missing(lam.interval)) {
      plot.set[xtmp > max(lam.interval)] <- FALSE
      plot.set[xtmp < min(lam.interval)] <- FALSE
      xtmp <- - xtmp[plot.set]
    }
    else xtmp <- -xtmp
    if (yvar == "coefficients") ytmp <- tmp$coefficients[plot.set, ]      
    if ((yvar == "newy") & is.matrix(newx)) ytmp <- tmp$newy[plot.set,]
    if ((yvar == "newy") & (! is.matrix(newx))) ytmp <- tmp$newy[plot.set]
    if (yvar == "lambda") ytmp <- tmp$lambda[plot.set]    
    if (yvar == "dimension") ytmp <- tmp$dimension[plot.set]
    if (yvar == "R-square") ytmp <- tmp$r.square[plot.set]
  }
  xmax <- max(xtmp)
  xmin <- min(xtmp)
  ymax <- max(ytmp)
  ymin <- min(ytmp)
  if (xvar == "lambda") {
    plot(c(xmin,xmax),c(ymin,ymax), xaxt = "n", xlab = xvar, ylab=yvar, type="l",lty=0, main=x$method)
    axis(1, at = axTicks(1), labels = abs(axTicks(1)) )
  } 
  else plot(c(xmin,xmax),c(ymin,ymax), xlab = xvar, ylab=yvar, type="l",lty=0, main=x$method)
  if ( ((yvar=="coefficients")|(yvar=="newy")) & is.matrix(ytmp)) {
     index <- 1: (dim(ytmp)[2])
     if (yvar=="coefficients") {
       f.set <- rep(FALSE, max(length(index),max(predictors)))
       f.set[predictors] <- TRUE
       f.set <- f.set[1:length(index)]
     }
     else f.set <- rep(TRUE, length(index))
     axis(4, at = ytmp[dim(ytmp)[1],f.set], labels = index[f.set])
     for (j in 1: (dim(ytmp)[2])) 
       if (f.set[j]) lines(xtmp, ytmp[,j]) 
  }
  else lines(xtmp, ytmp)
}

