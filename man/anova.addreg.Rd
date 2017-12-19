\name{anova.addreg}
\alias{anova.addreg}
\alias{anova.addreglist}
\title{
Analysis of Deviance for addreg Fits
}
\description{
Compute an analysis of deviance table for more than one
GLM fitted using \code{\link{addreg}}.
}
\usage{
\method{anova}{addreg}(object, ..., test = NULL)
}
\arguments{
  \item{object, \dots}{
	objects of class \code{"addreg"},
  typically the result of a call to \code{\link{addreg}},
  or a \code{\link{list}} of \code{objects} for the "\code{addreglist}"
  method.
}
  \item{test}{
	a character string, (partially) matching one
  of \code{"Chisq"}, \code{"LRT"}, \code{"Rao"}, \code{"F"}
  or \code{"Cp"}. See \code{\link{stat.anova}}.
}
}
\details{
Unlike \code{\link{anova.glm}}, specifying a single object
is not allowed.

The table has a row for the residual degrees of freedom and
deviance for each model. For all but the first model, the
change in degrees of freedom and deviance is also given.
(This only makes statistical sense if the models are
nested.) It is conventional to list the models from
smallest to largest, but this is up to the user.

Models where the MLE lies on the boundary of the parameter
space will be automatically removed from the list (with a
warning), because asymptotic results to not apply to such
models.

The table will optionally contain test statistics (and p-values) 
comparing the reduction in deviance for the row to
the residuals. Mallows' \emph{Cp} statistic is the residual
deviance plus twice the estimate of \eqn{\sigma^2} times
the residual degrees of freedom, which is closely related
to AIC. You can also choose \code{"LRT"} and \code{"Rao"}
for likelihood ratio tests and Rao's efficient score test.
The former is synonymous with \code{"Chisq"} (although both
have an asymptotic chi-square distribution).
}
\value{
An object of class \code{"anova"} inheriting from class
\code{"data.frame"}.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\section{Warning}{
  The comparison between two or more models will only be
  valid if they are fitted to the same dataset. This may be
  a problem if there are missing values and R's
  default of \code{na.action = na.omit} is used, and
  \code{anova} will detect this with an error.
}
\seealso{
\code{\link{addreg}}, \code{\link{anova.glm}},
\code{\link{anova}}
}
\examples{
## For an example, see example(addreg)
}
\keyword{ models }
\keyword{ regression }
