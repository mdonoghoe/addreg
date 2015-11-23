\name{addbin}
\alias{addbin}
\title{
Fitting Additive Binomial Regression Models
}
\description{
Workhorse function for \code{\link{addreg}} with
\code{binomial} family.
}
\usage{
addbin(y, x, start = NULL, control = list(), allref)
}
\arguments{
  \item{y}{
	binomial response. May be a single column of 0/1
  or two columns, giving the number of successes and
  failures.
}
  \item{x}{
	non-negative design matrix. Must have an
  intercept column.
}
  \item{start}{
	starting values for the parameters in the
  linear predictor.
}
  \item{control}{
	list of parameters for controlling the
  fitting process, passed to
  \code{\link{addreg.control}}.
}
  \item{allref}{
	a list of all parameterisations for this
  model, obtained from \code{\link{addreg.allref}}.
}
}
\details{
An additive binomial fit can be converted into an additive
Poisson fit via the multinomial--Poisson transformation 
(Baker, 1994). This function transforms the data as
described by Donoghoe and Marschner (2014) and passes it to \code{\link{addreg}}
with a Poisson family to get the maximum likelihood
estimate. The coefficients (and other values) from the
Poisson model are transformed back to relate to the
additive binomial model.

This is a workhorse function for \code{\link{addreg}} when
a \code{\link{binomial}} family is specified. It would not usually
be called directly.
}
\value{
A list of (most of) the components needed for an object of
class \code{"addreg"}; see \code{\link{addreg}} for
details.
}
\references{
Baker, S. G. (1994). The multinomial--Poisson
transformation. \emph{The Statistician} 43(4): 495--504.

Donoghoe, M. W. and I. C. Marschner (2014). Stable computational methods
for additive binomial models with application to adjusted risk differences. 
\emph{Computational Statistics and Data Analysis} 80: 184--196.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}

\seealso{
\code{\link{addreg}}
}
\keyword{ models }
\keyword{ regression }
