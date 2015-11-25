\name{addreg.control}
\alias{addreg.control}
\title{
Auxiliary for Controlling addreg Fitting
}
\description{
Auxiliary function for \code{\link{addreg}} fitting.
Typically only used internally by \code{\link{nnpois}} and 
\code{\link{nnnegbin}}, but may 
be used to construct a \code{control} argument to these functions.
}
\usage{
addreg.control(bound.tol = 1e-06, epsilon = 1e-10, maxit = 10000, trace = 0)
}
\arguments{
  \item{bound.tol}{
	positive tolerance specifying the
  interior of the parameter space. If the fitted model is
  more than \code{bound.tol} away from the boundary of the
  parameter space then it is assumed to be in the interior.
  This can allow the computational method to terminate
  early if an interior maximum is found.  No early
  termination is attempted if \code{bound.tol = Inf}.
}
  \item{epsilon}{
	positive convergence tolerance
  \eqn{\epsilon}; the estimates are considered to have 
  converged when \eqn{\sqrt{
  \sum (\theta_{old} - \theta_{new})^2} / \sqrt {\sum
  \theta_{old}^2} < \epsilon}, where \eqn{\theta} is the
  vector of parameter estimates. See \code{\link{conv.test}}.
}
  \item{maxit}{
	integer giving the maximum number of EM
  algorithm iterations for a given parameterisation.
}
  \item{trace}{
	number indicating level of output that
  should be produced. >= 1 gives output for each
  parameterisation, >= 2 gives output at each iteration.
}
}
\details{
This is used similarly to \code{\link{glm.control}}. The
\code{control} argument of \code{\link{addreg}} is by
default passed to the \code{control} argument of
\code{\link{nnpois}} or \code{\link{nnnegbin}}.

When \code{trace} is greater than zero, calls to
\code{\link{cat}} produce the output. Hence,
\code{\link{options}(digits = *)} can be used to increase
the precision.
}
\value{
A list with components named as the arguments.
}

\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}

\seealso{
\code{\link{glm.control}}, the equivalent function for
\code{\link{glm}} fitting.

\code{\link{nnpois}} and \code{\link{nnnegbin}}, the
functions used to fit \code{\link{addreg}} models.
}
\examples{
## Variation on example(glm.control) :

counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
oo <- options(digits = 12)
addreg.D93X <- addreg(counts ~ outcome + treatment, family = poisson, 
  trace = 2, epsilon = 1e-2)
options(oo)
coef(addreg.D93X)
}
\keyword{ optimize }
\keyword{ models }