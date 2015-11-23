\name{addreg.allref}
\alias{addreg.allref}
\title{
Parameterisation for CEM Algorithm
}
\description{
A workhorse function for \code{\link{addreg}},
\code{addreg.allref} takes the formula and data for an
identity-link GLM and produces a list of all
parameterisations needed for the associated CEM algorithm.
}
\usage{
addreg.allref(object, data = environment(object), mono, family, start = NULL)
}
\arguments{
  \item{object}{
	a model \code{\link{formula}} or a
  \code{\link{terms}} object for the \code{\link{addreg}} model.
}
  \item{data}{
	a data frame created with
  \code{\link{model.frame}}. If another sort of object,
  \code{\link{model.frame}} is called first.
}
  \item{mono}{
	a vector indicating which terms should be
  restricted to have a monotonically non-decreasing
  relationship with the outcome.
}
  \item{family}{
	the result of a call to a family function.
  Its component \code{$family} must be one of
  \code{"poisson"}, \code{"negbin1"} or \code{"binomial"}.
}
  \item{start}{
	starting values for the parameters in the linear predictor.
}
}
\details{
In the CEM algorithms employed by \code{\link{addreg}}, the
parameter space is partitioned into a collection of
restricted parameter spaces (see Marschner, 2014).
\code{addreg.allref} finds the list of possible
parameterisations of each term in the model.

If a term \code{x} has a \code{TRUE} value for
\code{\link{is.factor}(x)}, \code{\link{is.character}(x)}
or \code{\link{is.logical}(x)}, it is considered to be a
categorical covariate. For Poisson and negative binomial models, this has a
parameterisation for each level of the factor, and for
binomial models, every permutation of the levels must be
considered (see Donoghoe and Marschner, 2014).

Otherwise the covariate is considered to be continuous, in
which case it has two possible parameterisations, relating
to the minimum and maximum observed values.

If a covariate is restricted to be monotonic via the
\code{mono} argument, it has only one parameterisation.

The \code{\link{addreg}} function considers all possible combinations of
the parameterisations of each covariate, and uses
\code{\link{addreg.design}} to create the appropriate
non-negative design matrix to be used in the EM algorithm.
}
\value{
A list with components: 
\item{allref}{a named list, with
one component for each term in the model. Each component is
itself a list, whose components are each of the
parameterisations for that term. 

If \code{start} was specified, the first component for each
term will correspond to the parameterisation specified
by \code{start}.} 
\item{terms}{the \code{\link{terms}} component of \code{object}.}
\item{data}{the object passed into the \code{data}
argument, or the result of calling \code{\link{model.frame}} with \code{data}.}
\item{monotonic}{a named logical vector indicating which
components of \code{terms} are restricted to be
monotonically non-decreasing.}
\item{start.new}{a reparameterised version of \code{start}, corresponding
to the first parameterisation in \code{allref}. \code{NULL} if \code{start}
was not supplied.}
}
\references{
Donoghoe, M. W. and I. C. Marschner (2014). Stable computational methods
for additive binomial models with application to adjusted risk differences. 
\emph{Computational Statistics and Data Analysis} 80: 184--196.

Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}

\seealso{
\code{\link{addreg}}
}

\keyword{ misc }