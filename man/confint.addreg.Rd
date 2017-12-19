\name{confint.addreg}
\alias{confint.addreg}
\title{
Confidence Intervals for addreg Model Parameters
}
\description{
Computes confidence intervals for one or more parameters in
a fitted \code{\link{addreg}} model.
}
\usage{
\method{confint}{addreg}(object, parm, level = 0.95, ...)
}
\arguments{
  \item{object}{
	a fitted model object, resulting from a
  call to \code{\link{addreg}}.
}
  \item{parm}{
	a specification of which parameters are to be
  given confidence intervals, either a vector of numbers or
  a vector of names. If missing, all parameters are
  considered.
}
  \item{level}{
	the confidence level required.
}
  \item{\dots}{
	additional argument(s) passed to \code{\link{confint.default}}.
}
}
\details{
Calculates confidence intervals for model parameters assuming asymptotic normality, using \code{\link{vcov.addreg}(object)}. As
such, if the MLE is on the boundary of the parameter space,
(i.e. \code{object$boundary == TRUE}) the normality assumption 
is invalid and \code{NA} is returned.
}
\value{
A matrix (or vector) with columns giving lower and upper
confidence limits for each parameter. These will be
labelled as (1-level)/2 and 1-(1-level)/2 in \% (by default
2.5\% and 97.5\%).
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\seealso{
\code{\link{confint.default}}, \code{\link{vcov.addreg}}
}
\examples{
## For an example, see example(addreg)
}
\keyword{ models }
