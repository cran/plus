\name{print.plus}
\alias{print.plus}
\title{Print plus() steps when predictors are added or removed.}
\description{
Print plus steps when predictors are added or removed, beginning with a description of the method 
such as LASSO, MC+ or SCAD. 
}
\usage{
\method{print}{plus}(x, print.moves = 20, ...)
}
\arguments{
  \item{x}{a plus object}
  \item{print.moves}{the number of lines printed. Default is 20.}
  \item{\dots}{Additonal arguments for generic methods}
}
\details{
When a predictor is added or removed in the solution path, print one line providing the 
step number, the action, and the predicor involved.
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
print(object, print.moves=30)
detach(sp500)
}
\keyword{print}
\keyword{methods}
