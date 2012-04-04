plus.new.eta <-
function(eta, etatil, cross,p) {
  eta <- c(eta,0) # put the dummy in the zero interval
  eta[ etatil[1,cross] ] <- etatil[2, cross]
  eta <- eta[1:p] # remove dummy
  return(eta)
}

