\name{gamma1}
\alias{gamma1}
\title{
Family Functions for Gamma 1 Models
}
\description{
Specifies the information required to fit a gamma 1 model.
}
\usage{
gamma1(link, phi = stop("'phi' must be given"))
}
\arguments{
  \item{link}{
	included for compatibility with \code{\link[stats]{family}}.
  For \code{\link{addreg}} models, this will always be
  \code{"identity"}.
}
  \item{phi}{
	the value of the scale parameter of the gamma 1 distribution (see
  "Details"). This can be set to \code{NA} for initialisation,
  but during estimation the family should be updated with the
  current estimate, and must be strictly positive.
}
}
\details{
  The gamma 1 distribution can be parameterised in
  terms of a mean \eqn{\mu} and scale parameter \eqn{\phi} (the 
  \code{phi} argument of this function), such that if
  \eqn{Y \sim Gamma1(\mu, \phi)}, then \eqn{E(Y) = \mu} and
  \eqn{Var(Y) = \phi\mu}.
  
  These can be related to the \code{shape} and \code{scale} arguments
  of the \code{\link{GammaDist}} functions by \code{shape} = \eqn{\mu / \phi}
  and \code{scale} = \eqn{1 / \phi}.
}
\value{
An object of class \code{"family"}: see \code{\link[stats]{family}} for
full details. Note that when the estimate of \code{phi} is updated in
a model, this \code{family} object must be reloaded using the new
estimate.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}
\keyword{ models }