\name{plot.plus}
\alias{plot.plus}
\title{Plot method for plus objects}
\description{
Produce a plot from a plus solution path. 
}
\usage{
\method{plot}{plus}(x, xvar=c("lam","step"), yvar=c("coef","newy","lam","dim","R-sq"), 
	newx, step.interval, lam.interval, predictors, ...)
}
\arguments{
  \item{x}{a plus object}
  \item{xvar}{penalty level or plus step as the variable for the horizontal axis in the plot. Default is "lam".}
  \item{yvar}{paths of coefficients, predictions, penalty level, the number of nonzero coefficients or R-square 
			as the variable for the vertical axis in the plot. Default is "coef"}
  \item{newx}{x values at which the fit is required. If newx is not set and yvar is "newy", no plot is produced.}
  \item{step.interval}{lower and upper bounds of the x-axis in the plot when xvar is "step". Defult covers all 
			steps in the computed path. }      
  \item{lam.interval}{lower and upper bounds of the x-axis in the plot when xvar is "lam".  Default covers 
			all penalty levels in the computed path.}
  \item{predictors}{a subset of predictors for which coefficients are plotted. Default is the entire set of predictors.}
  \item{\dots}{Additonal arguments for generic methods}  
}
\details{
The fitted coefficients and penalty levels are linear between two consecutive turning points in the plus path 
so that exact values of "coef", "newy" and "lam" are ploted when xvar is set as "step". For concave penalties, 
the solution path is not necessarily monotone in penalty level. Since the extracted coefficients for a particular 
given penalty level is defined as the first point at which the solution path hits the given penalty level, the 
"coef" and "newy" plotted as approximations as the linear interpolation of their actual values at specifiec 
lam when xvar is set as "lam".  
}
\value{NULL}
\references{
Zhang, C.-H. (2007). Penalized linear unbiased selection. Technical Report No. 2007-003. 
Department of Statistics, Rutgers University.}

\author{Cun-Hui Zhang and Ofer Melnik}
\examples{
data(sp500)
attach(sp500)

x <- sp500.percent[,3: (dim(sp500.percent)[2])] 
y <- sp500.percent[,1]

object <- plus(x,y,method="mc+")
plot(object)
detach(sp500)
}
\keyword{hplot}
\keyword{methods}
