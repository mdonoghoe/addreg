\name{negbin1}
\alias{negbin1}
\title{
Family Functions for Negative Binomial 1 Models
}
\description{
Specifies the information required to fit a negative binomial 1 (NB1) model.
}
\usage{
negbin1(link, phi = stop("'phi' must be given"))
}
\arguments{
  \item{link}{
	included for compatibility with \code{\link[stats]{family}}.
  For \code{\link{addreg}} models, this will always be
  \code{"identity"}.
}
  \item{phi}{
	the value of the scale parameter of the NB1 distribution (see
  "Details"). This can be set to \code{NA} for initialisation,
  but during estimation the family should be updated with the
  current estimate, and must be strictly positive.
}
}
\details{
  The NB1 distribution can be parameterised in
  terms of a mean \eqn{\mu} and scale parameter \eqn{\phi} (the 
  \code{phi} argument of this function), such that if
  \eqn{Y \sim NB1(\mu, \phi)}, then \eqn{E(Y) = \mu} and
  \eqn{Var(Y) = (1+\phi)\mu}.
  
  These can be related to the \code{size} and \code{prob} arguments
  of the \code{\link{NegBinomial}} functions by \code{size} = \eqn{\mu / \phi}
  and \code{prob} = \eqn{1 / (1 + \phi)}.
}
\value{
An object of class \code{"family"}: see \code{\link[stats]{family}} for
full details. Note that when the estimate of \code{phi} is updated in
a model, this \code{family} object must be reloaded using the new
estimate.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\keyword{ models }
