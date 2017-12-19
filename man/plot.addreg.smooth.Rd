\name{plot.addreg.smooth}
\alias{plot.addreg.smooth}
\title{
Default addreg.smooth Plotting
}
\description{
Takes a fitted \code{addreg.smooth} object produced by \code{\link{addreg.smooth}} and plots
the component smooth functions that make it up, on the scale of the linear predictor,
for specified values of the other covariates.
}
\usage{
\method{plot}{addreg.smooth}(x, type = c("response", "link"), at = data.frame(), 
     knotlines = TRUE, nobs = 1000, ...)
}
\arguments{
  \item{x}{
	a fitted \code{addreg.smooth} object as produced by \code{\link{addreg.smooth}}.
}
  \item{type}{
	the type of prediction required. Note that, unlike \code{\link{predict.addreg.smooth}},
	\code{"terms"} is not a valid option. Also, because \code{addreg.smooth} only applies
	identity-link models, \code{"response"} and \code{"link"} will have the same results ---
	they are included for consistency.
}
  \item{at}{
	a data frame containing the values at which the prediction should be evaluated. The columns
	must contain the covariates in the model, and several rows may be provided (in which case, multiple
	lines are drawn on the same plot). Cannot be missing or \code{NULL}.
}
  \item{knotlines}{
	logical; if vertical lines should be drawn on the plot to indicate the locations
	of the knots for B-spline terms.
}
  \item{nobs}{
	the number of points which should be used to create the curve. These are placed evenly along
	the range of the observed covariate values from the original model.
}
  \item{\dots}{
	other graphics parameters to pass on to plotting commands (note: some will not work).
}
}
\details{
For each smooth covariate in the model of \code{x}, \code{\link{predict.addreg.smooth}}
is used to obtain predicted values for the range of that covariate, with the other
covariates remaining fixed at their values given in \code{at}. Several rows may be provided
in \code{at}, in which case, one curve is drawn for each, and they are coloured using
\code{\link{rainbow}(nrow(at))}. If the model contains a single smooth covariate and no other
covariates, \code{at} may be provided as an empty data frame, \code{data.frame()}.
}
\value{
The function simply generates plots.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\note{
If this function is too restrictive, it may be easier to use \code{\link{predict.addreg.smooth}}
to get predictions for the dataset of your choice, and do the plotting manually.
}
\seealso{
\code{\link{addreg.smooth}}, \code{\link{predict.addreg.smooth}}
}
\examples{
## For an example, see example(addreg.smooth)
}
\keyword{ models }
\keyword{ regression }
\keyword{ smooth }
