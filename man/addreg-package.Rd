\name{addreg-package}
\alias{addreg-package}
\docType{package}
\title{
Additive Regression for Discrete Data
}
\description{
Methods for fitting identity-link GLMs and GAMs to discrete data, using EM-type
algorithms with more stable convergence properties than standard methods.
}
\details{
\tabular{ll}{
Package: \tab addreg\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2016-09-21\cr
License: \tab GPL (>= 2)\cr
}
This package provides methods to fit generalised linear models (GLMs) and generalised additive
models (GAMs) with identity link functions to discrete data using binomial, Poisson
and negative binomial models. It is planned that future versions will incorporate other types 
of discrete data models, such as multinomial regression.

The package has two primary functions: \code{\link{addreg}} and \code{\link{addreg.smooth}}, 
together with various supporting functions. It is useful in two main situations. The first is 
when a standard GLM routine, such as \code{\link{glm}}, fails to converge with such a model. 
The second is when a flexible semi-parametric component is desired in these models. One of the 
main purposes of this package is to provide parametric and semi-parametric adjustment of risk 
differences and rate differences.

Two approaches based on the EM algorithm can be used. These methods accommodate the required
parameter constraints and are more stable than iteratively reweighted least squares. In a 
combinatorial EM algorithm (Marschner, 2014), a collection of restricted parameter spaces is 
defined which covers the full parameter space, and the EM algorithm is applied within each
restricted parameter space in order to find a collection of restricted maxima of
the log-likelihood function, from which can be obtained the global maximum over the
full parameter space.

In the expanded EM approach, additional parameters are added to the model,
and an EM algorithm finds the MLE of this overparameterised model by imposing
constraints on each individual parameter. This requires a single application of
the EM algorithm.

In either case, the EM algorithm may be accelerated by using the capabilities of the
\pkg{\link[turboEM:turboem]{turboEM}} package.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}

Maintainer: Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}
\references{
Donoghoe, M. W. and I. C. Marschner (2014). Stable computational methods
for additive binomial models with application to adjusted risk differences. 
\emph{Computational Statistics and Data Analysis} 80: 184--196.

Donoghoe, M. W. and I. C. Marschner (2015). Flexible regression models for rate 
differences, risk differences and relative risks. \emph{International Journal of 
Biostatistics} 11(1): 91--108.

Marschner, I. C. (2010). Stable computation of maximum likelihood estimates
in identity link Poisson regression. \emph{Journal of Computational and 
Graphical Statistics} 19(3): 666--683.

Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.
}
\keyword{ package }
\keyword{ regression }
\seealso{
\code{\link{glm}}
}
\examples{
## For examples, see example(addreg) and example(addreg.smooth)
}
