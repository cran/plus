plus.penalty <-
function(m,a=0) {
### lasso for m=1, MCP for m=2, SCAD for m=3; default gamma = a = 3.7
### vectors v and t with dummy are used for lasso to maintain vector attribute 
m <- round(m)
if (m < 1) m <- 2
if (m > 3) m <- 2
if (m == 2) {
  if (a==0) a <- 3.7
  v <- c(1/a,0)
  t <- c(0,a)
}
if (m == 1) {
  a <- 0
  v <- c(0,0) 
  t <- c(0,0)
}
if (m==3) {
  if (a==0) a <- 3.7
  v <- c(0,1/(a-1),0) 
  t <- c(0,1,a)
}
return(list(m=m, gamma=a, v=v,t=t))
}

