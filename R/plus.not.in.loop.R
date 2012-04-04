plus.not.in.loop <-
function( k, new.eta,e) {
## input: eta, (k+1)Xp; new.eta, 1Xp
## output is T iff new.eta is not any of eta[i,], i\le k
  new.is.new <- TRUE
  for (i in 1:k) 
    if( all(e$eta[i,]==new.eta)) return(FALSE)
  return(new.is.new)
}

