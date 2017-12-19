\name{interpret.addreg.smooth}
\alias{interpret.addreg.smooth}
\title{
Interpret an addreg.smooth Formula
}
\description{
This is an internal function of package \code{\link[=addreg-package]{addreg}}.
It is a service routine for \code{\link{addreg.smooth}} which interprets the
smooth parts of the model formula and returns modified formulas to be used
in the fitting functions.

Not normally called directly.
}
\usage{
interpret.addreg.smooth(formula)
}
\arguments{
  \item{formula}{
	A formula as supplied to \code{\link{addreg.smooth}}, which includes
	at least one \code{\link{B}} or \code{\link{Iso}} term.
}
}
\value{
A list with components:
\item{full.formula}{a \code{\link{formula}} object which is the same as
the \code{formula} supplied, but with additional arguments removed from
the smooth terms. E.g. \code{B(x, knot.range = 0:2)} would appear as
\code{B(x)} in this formula.}
\item{fake.formula}{a \code{\link{formula}} object which is the same as
the \code{formula} supplied, but with smooth terms replaced by their
covariates alone. E.g. \code{B(x, knot.range = 0:2)} would appear as 
\code{x} in this formula. Used to construct the model matrix.}
\item{smooth.spec}{a named list containing the results of \code{\link{eval}}uating
the smooth terms. See \code{\link{B}} and \code{\link{Iso}} for details.}
\item{smooth.ind}{a vector containing the indices of the smooth components in
the \code{formula}.}
\item{terms}{the result of running \code{\link{terms.formula}(formula, specials = c("B", "Iso"))}.}
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\seealso{
\code{\link{addreg.smooth}}
}
\examples{
# Specify a smooth model with knot.range
res <- interpret.addreg.smooth(y ~ B(x, knot.range = 0:2) + x2)
# The knot.range is removed from the full.formula...
print(res$full.formula)
# ...but is stored in the $smooth.spec component of the result:
print(res$smooth.spec$x$knot.range)
}
\keyword{ smooth }
\keyword{ models }
