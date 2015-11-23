\name{vcov.addreg}
\alias{vcov.addreg}
\title{
Calculate Variance-Covariance Matrix for a Fitted addreg Model Object
}
\description{
Returns the variance-covariance matrix of the main parameters of a fitted \code{addreg} model object.
}
\usage{
\method{vcov}{addreg}(object, ...)
}
\arguments{
  \item{object}{
	an object of class \code{"addreg"}, usually from a call to \code{\link{addreg}}
	or \code{\link{addreg.smooth}}.
}
  \item{\dots}{
    additional arguments for method functions.
}
}
\details{
An equivalent method to \code{\link{vcov}}, to use with \code{\link{addreg}} models.
}
\value{
A matrix of the estimated covariances between the parameter estimates in the linear or 
non-linear predictor of the model. This should have row and column names corresponding 
to the parameter names given by the \code{\link{coef}} method.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}
\note{
If \code{object$boundary == TRUE}, the standard errors of the coefficients
are not valid, and a matrix of \code{NaN}s is returned.
}

\seealso{
\code{\link{summary.addreg}}, \code{\link{vcov.glm}}
}
\examples{
## For an example, see example(addreg)
}
\keyword{ models }
\keyword{ regression }
