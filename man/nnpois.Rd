\name{nnpois}
\alias{nnpois}
\title{
EM Algorithm for Identity-link Poisson GLM
}
\description{
Finds the maximum likelihood estimate of an identity-link Poisson GLM
using an EM algorithm,
where each of the coefficients is restricted to be non-negative.
}
\usage{
nnpois(y, x, standard, offset, start, control = addreg.control(),
       accelerate = c(""em", "squarem", "pem", "qn"),
       control.accelerate = list(list()))
}
\arguments{
  \item{y}{
	non-negative integer response vector.
}
  \item{x}{
	non-negative covariate matrix.
}
  \item{standard}{
	standardising vector, where each element is a positive constant that
	(multiplicatively) standardises the fitted value of the corresponding
	element of the response vector. The default is a vector of ones.
}
  \item{offset}{
	non-negative additive offset vector. The default is a vector of zeros.
}
  \item{start}{
	starting values for the parameter estimates. Each element must be
	greater than \code{control$bound.tol}.
}
  \item{control}{
	an \code{\link{addreg.control}} object, which controls the fitting process.
}
  \item{accelerate}{
    a character string that determines the acceleration
  algorithm to be used, (partially) matching one of \code{"em"} (no acceleration -- the default), 
  \code{"squarem"}, \code{"pem"} or \code{"qn"}. See \code{\link[turboEM]{turboem}}
  for further details. Note that \code{"decme"} is not permitted.
}
  \item{control.accelerate}{
	a list of control parameters for the acceleration algorithm. See \code{\link[turboEM]{turboem}}
	for details of the parameters that apply to each algorithm. If not specified, the defaults are used.
}
}
\details{
This is a workhorse function for \code{\link{addreg}}, and runs the EM algorithm to find the 
constrained non-negative MLE associated with an identity-link Poisson GLM. See Marschner (2010)
for full details.
}
\value{
A list containing the following components
\item{coefficients}{the constrained non-negative maximum likelihood estimate of the parameters.}
\item{residuals}{the residuals at the MLE, that is \code{y - fitted.values}}
\item{fitted.values}{the fitted mean values.}
\item{rank}{the number of parameters in the model (named ``\code{rank}" for compatibility ---
we assume that models have full rank)}
\item{family}{included for compatibility --- will always be \code{\link{poisson}(identity)}.}
\item{linear.predictors}{included for compatibility --- same as \code{fitted.values} (as this is
an identity-link model).}
\item{deviance}{up to a constant, minus twice the maximised log-likelihood.}
\item{aic}{a version of Akaike's \emph{An Information Criterion},
minus twice the maximised log-likelihood plus twice the 
number of parameters.}
\item{aic.c}{a small-sample corrected
version of Akaike's \emph{An Information Criterion}
(Hurvich, Simonoff and Tsai, 1998).}
\item{null.deviance}{the deviance for the null model, comparable with \code{deviance}.
The null model will include the offset and an intercept.}
\item{iter}{the number of iterations of the EM algorithm used.}
\item{weights}{included for compatibility --- a vector of ones.}
\item{prior.weights}{included for compatibility --- a vector of ones.}
\item{standard}{the \code{standard} vector passed to this function.}
\item{df.residual}{the residual degrees of freedom.}
\item{df.null}{the residual degrees of freedom for the null model.}
\item{y}{the \code{y} vector used.}
\item{converged}{logical. Did the EM algorithm converge 
(according to \code{\link{conv.test}})?}
\item{boundary}{logical. Is the MLE on the boundary of the parameter
space --- i.e. are any of the \code{coefficients < control$bound.tol}?}
\item{loglik}{the maximised log-likelihood.}
\item{nn.design}{the non-negative \code{x} matrix used.}
}
\references{
Hurvich, C. M., J. S. Simonoff and C.-L. Tsai (1998). Smoothing parameter 
selection in nonparametric regression using an improved Akaike 
information criterion.
\emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)} 60(2): 271--293.

Marschner, I. C. (2010). Stable computation of maximum likelihood estimates
in identity link Poisson regression. \emph{Journal of Computational and 
Graphical Statistics} 19(3): 666--683.

Marschner, I. C., A. C. Gillett and R. L. O'Connell (2012).
Stratified additive Poisson models: Computational methods
and applications in clinical epidemiology. \emph{Computational
Statistics and Data Analysis} 56(5): 1115--1130.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}.

This function is based on code from Marschner, Gillett and O'Connell (2012) written by Alexandra Gillett.
}
\keyword{ regression }
