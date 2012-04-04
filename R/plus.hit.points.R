plus.hit.points <-
function(lam, lam.path) {
## find the first point the plus path hits the specified penalty levels lam 
## output:  k1, k2, w for beta <- w1*beta.path[k1,] + (1-w1)*beta.path[k2,]
## with k1 = k2 at the beginning and end, k2 = k1+1 otherwise
	if (is.null(lam) || !is.numeric(lam) || is.nan(sum(lam)))
	{        
		lam <- sort(lam.path, decreasing = TRUE)
		cat("Warning message: ")
		cat("invalid lam; lam.path used\n")
	}
	if (max(lam)< min(lam.path))
	{        
		lam <- sort(lam.path, decreasing = TRUE)
		cat("Warning message: ")
		cat("lam not reached by the plus path; lam.path used\n")
	}
	lam <- mapply(function(x) min(x,max(lam.path)),lam)
	if (min(lam) < min(lam.path))
	{
		lam <- lam[lam>=min(lam.path)]
		cat("Warning message: ")
		cat("some lam not reached by the plus path and dropped\n")
	}
	lam.lesser.path <- mapply(function(x) x<=lam.path,lam)
	lam.greater.path <- mapply(function(x) x>=lam.path,lam)
	hit <- (lam.lesser.path[1:(length(lam.path)-1),] & 
			lam.greater.path[2:length(lam.path),])
	k1 <- apply(as.matrix(hit),2,function(bool.vec) which(bool.vec)[1])
	k2 <- k1+1
	w1 <- (lam - lam.path[k2])/(lam.path[k1]-lam.path[k2])
	return(list("k1"=k1,"k2"=k2,"w1"=w1,"total.hits"=length(lam)))
}

