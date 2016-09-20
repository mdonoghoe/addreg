\name{contr.isotonic}
\alias{contr.isotonic}
\alias{contr.opisotonic}
\title{
Contrast Matrix for Isotonic Covariate
}
\description{
Return something similar to a contrast matrix for a categorical
covariate that we wish to be monotonically non-decreasing
in a specified order.
}
\usage{
contr.isotonic(n, perm, contrasts = TRUE, sparse = FALSE)

contr.opisotonic(n, perm, contrasts = TRUE, sparse = FALSE)
}
\arguments{
  \item{n}{
	a vector of levels for a factor, or the number
  of levels.
}
  \item{perm}{
	a permutation of the levels of \code{n} (or
  of the numbers \code{1:n}), which define the order in
  which the coefficients must be monotonically
  non-decreasing.
}
  \item{contrasts}{
	a logical indicating whether constrasts
  should be computed.
}
  \item{sparse}{
	included for compatibility reasons. Has no
  effect.
}
}
\details{
\code{contr.isotonic} is used in creating the design matrix
for categorical covariates with a specified order under a
particular parameterisation. For Poisson and negative binomial models, this
occurs if a categorical covariate is defined as monotonic;
for binomial models, each parameterisation defines a
permutation of the levels that must be monotonically
increasing.

For overparameterised binomial models, the design matrix for
categorical covariates must include isotonic-style dummy 
covariates for every possible permutation of the levels. This
is the function of \code{contr.opisotonic}.

In the order specified by \code{perm}, the coefficient
associated with each level is the sum of increments between
the preceding levels. That is, the first level is defined
as \eqn{0}, the second as \eqn{0 + d_2}, the third as \eqn{0 + d_2 + d_3}, and
so on.  In fitting the model, these increments are
constrained to be non-negative.

Note that these are not `contrasts' as defined in the
theory for linear models; rather this is used to define the
\code{contrasts} attribute of each variable so that
\code{\link{model.matrix}} produces the desired design
matrix.
}
\value{
A matrix with \code{n} rows and \code{k} columns, with
\code{k=n-1} if \code{contrasts} is \code{TRUE} and
\code{k=n} if \code{contrasts} is \code{FALSE}.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}
\seealso{
\code{\link{model.matrix}}, which uses
\code{contr.isotonic} to create the design matrix.

\code{\link{contr.treatment}}, \code{\link{contrasts}} for
their usual use in regression models.
}
\examples{
contr.isotonic(4,1:4)
contr.isotonic(4,c(1,3,2,4))

# Show how contr.isotonic applies within model.matrix
x <- factor(round(runif(20,0,2)))
mf <- model.frame(~x)
contrasts(x) <- contr.isotonic(levels(x), levels(x))
model.matrix(mf)
}
\keyword{ design }