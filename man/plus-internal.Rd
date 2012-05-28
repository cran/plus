\name{plus-internal}
\title{Internal plus functions}
\alias{plus.single.step}
\alias{plus.check.eta}
\alias{plus.hit.points}
\alias{plus.new.eta}
\alias{plus.not.in.loop}
\alias{plus.penalty}
\alias{plus.recip}
\alias{plus.penalty}
\alias{plus.knots}
\alias{plus.penalty}
\description{Internal plus functions}
\usage{
plus.single.step(k, etatil,  z, b, tau, m, v, t, eps=1e-15, e)
plus.check.eta(old.eta, new.eta, etatil, xi, z, v, eps=1e-15, e)
plus.hit.points(lam, lam.path)
plus.new.eta(eta, etatil, cross,p)
plus.not.in.loop( k, new.eta,e)
plus.penalty(m,a=0)
plus.recip(Q, eps=1e-15)
plus.knots(t=0,m=0)
}
\author{Ofer Melnik and Cun-Hui Zhang, with contributions by Stephanie Zhang}
\details{
  These are not to be called by the user.
}
\keyword{internal}
