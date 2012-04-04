print.plus <-
function(x, print.moves = 20, ...) { 
  cat("\n")
  cat("Sequence of", x$method, "moves:", "\n","\n")
  print.moves <- round(print.moves)
  if (print.moves < 1) print.moves <- 20
  names <- dimnames(x$x)
  m <- 1
  k <- 1
  exit.while <- FALSE
  while (exit.while == FALSE) {
    var.change <- sign(abs(sign(x$eta[k+1,]))-abs(sign(x$eta[k,])))
    if (sum(abs(var.change))==1) {
      if(sum(var.change)==1)
        action <-"added "
      else
        action <-"removed "   
      var_name <- (order(-abs(var.change))[1])
      if (!is.null(names))
        var_name <- names[[2]][var_name]
      cat("Step",k,": ",action,var_name,"\n")
      m <- m+1
      if (m>print.moves) exit.while <- TRUE
    }
    k <- k+1
    if (k > x$total.steps) exit.while <- TRUE
  } # end while
}

