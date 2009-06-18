\name{plus}
\alias{plus}
\title{
Fits linear regression with a quadratic spline penalty, including the Lasso, MC+ and SCAD. 
}
\description{
The algorithm generates a piecewise linear path of coefficients and penalty levels
as critical points of a penalized loss in linear regression, starting with zero
coefficients for infinity penalty and ending with a least squares fit for zero penalty.
It is an extension of the LARS algorithm from the absolute value penalty to quadratic
spline penalties. 
}
\usage{
plus(x,y, method = c("lasso", "mc+", "scad", "general"), m=2, gamma,v,t,
   monitor=FALSE, normalize = TRUE, intercept = TRUE,
   Gram, use.Gram = FALSE, eps=1e-15, max.steps=500, lam)
}
\arguments{
\item{x}{
predictors, an n by p matrix with n > 1 and p > 1. 
}
\item{y}{
response, an n-vector with n > 1.
}
\item{method}{
c("lasso", "mc+", "scad", "general"); the LASSO penalty is specified by m = 1, MC+ is
specified by m = 2 and gamma > 0, SCAD by m = 3 and gamma > 1.  A general
quadratic penalty is specified by m-vectors v and t.
}
\item{m}{
number of knots with a quadratic spline penalty: m = 1 for Lasso, m = 2 for MC+, m = 3 for SCAD. Default is m = 2. 
}
\item{gamma}{
the largest knot of a quadratic spline penalty, say rho(.); gamma = 0 for lasso. 
}
\item{v}{
m-vector giving the negative second derivative rho(.) of the penalty between
two knots or beyond gamma. 
}
\item{t}{
m-vector giving the discontinuities of the derivatives of the penalty function rho(.) as
knots, including 0 as a knot. 
}
\item{monitor}{
If TRUE, plus prints out its progress when variables move in and out of the active set.
Default is FALSE. 
}
\item{normalize}{
If TRUE, each variable is standardized to have unit mean squares,
otherwise it is left alone. Default is TRUE. 
}
\item{intercept}{
If TRUE, an intercept is included in the model (and not penalized),
otherwise no intercept is included. Default is TRUE. 
}
\item{Gram}{
The X'X matrix; useful for repeated runs (e.g. bootstrap) where a large X'X stays the same. 
}
\item{use.Gram}{
When p is very large, you may not want PLUS to precompute the entire Gram matrix.
Default is FALSE. 
}
\item{eps}{
An effective zero. 
}
\item{max.steps}{
Limit the number of steps taken. Default is 500. There can be many more steps than
n or p since variables can be removed and added as the algorithm proceeds. Users
should check if the desired penalty level is reached if PLUS ends in the maximum step. 
}
\item{lam}{
A decreasing sequence of nonnegative numbers as penalty levels for which penalized
estimates of coefficients are generated. Default is the vector of ordered penalty levels
at the turning points of the computed path. If lam is set, the computation stops when the
path first hits the minimum of lam. The scale of lam is determined by the penalized loss
sum((y - x %*%beta)2)/(2n) + sum(lam2 * rho(abs(beta)/lam)) for p-vectors beta
}

}
\value{
A "plus" object is returned, for which print, predict, coef and plot methods exist. In addition 
to arguments x, y, max.steps, and the used values of method, gamma and lam, the object 
contains the following items: 

Some significant components of the object are:
\item{v}{matrix with rows as p-vectors indicating the parallelepipeds in which the computed path lives}
\item{beta.path}{Tmatrix with rows as p-vectors of regression coefficients at the turning points of the solution path}
\item{lam.path}{penalty levels at the turning points of the computed path. When the penalty function is concave, 
	lam.path may not be a decreasing sequence but always takes nonnegative values.}
\item{beta}{matrix with rows as p-vector of coefficients when the solution path first hits lam}
\item{lam}{the specified penalty levels hit by lam.path. This may not be the same as argument lam if 
	the minimum of the argument is not reached by the computed solution path. }
\item{dim}{the number of nonzero beta}
\item{r.square}{R-square values for beta}
\item{total.hits}{length of output lam}
\item{total.steps}{total number of steps executed, the same as the total number of segments in the 
	computed solution path. With zero as the first coefficient vector, beta.path contains 
	one more vector than total.steps. }
\item{full.path}{TRUE if zero penalty is reached.}
\item{forced.stop}{TRUE if PLUS is forced to stop due to reasons other than reaching max.steps or the 
	minimum of argument lam.} 
\item{singular.Q}{TRUE if PLUS is forced to stop when a matrix is not invertible.}  	
}
\details{
PLUS is described in detail in Zhang (2007). It computes a complete path of crititcal points of a 
penalised squared loss emcompassing from zero for infinite penalty to a lease squares fit for zero 
penalty, including possible multiple local minima for each penalty level. 
}
\references{
Zhang, C.-H. (2007). Penalized linear unbiased selection. Technical Report No. 2007-003. 
Department of Statistics, Rutgers University.
}
\author{Cun-Hui Zhang and Ofer Melnik}
\seealso{
print, plot, and predict methods
}
\examples{
data(sp500)
attach(sp500)
x <- sp500.percent[,3: (dim(sp500.percent)[2])] 
y <- sp500.percent[,1]

par(mfrow=c(2,3))
object <- plus(x,y,method="lasso")
plot(object)
plot(object, yvar="dim")
plot(object, yvar="R-sq")
object <- plus(x,y,method="mc+")
plot(object)
plot(object, yvar="dim")
plot(object, yvar="R-sq")
detach(sp500)
}
\keyword{regression}


