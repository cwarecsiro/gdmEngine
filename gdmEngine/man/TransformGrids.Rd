% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TransformGrids.R
\name{TransformGrids}
\alias{TransformGrids}
\title{Transform grids}
\usage{
TransformGrids(gdm.model, env.grids.stk, extrap.method = "Conservative",
  output.folder = NULL, output.name = "tran", verbose = TRUE)
}
\arguments{
\item{gdm.model}{(gdm model object) A GDM model object for which transformed grids of each predictor will be generated.}

\item{env.grids.stk}{(raster stack) A raster stack holding grids for all the predictors to be transformed.}

\item{extrap.method}{(string) The method used to derive transformed values when the raw values are outside the range used to train the model. Options are: "Clamped" (transformed values clamped at the extreme value in the training data); "End10" (transformed values extrapolated based on the slope of the spline function in the final (or initial) 10 percent range of values used to train the model), "WholeGrad" (transformed values extrapolated based on the slope of the spline function across the whole gradient of values used to train the model), "Conservative" (the lowest spline slope out of End10 and WholeGrad), (default = "Conservative").}

\item{output.folder}{(string) A folder to save the outputs to. If none specified, no file is written.}

\item{output.name}{(string) A name to use in saving the outputs. Default: 'pairs_table_dissim'.}

\item{verbose}{(boolean) Print messages to console. Default TRUE.}
}
\value{
Transformed grids are written to the specified output filepaths.
}
\description{
The R function to call in order to transform big grids. This in turn calls a Rcpp function, that does the magic
}
\examples{
TransformGrids(My.model, My.env.stk, output.folder = 'C:/Users/processed_data', output.name = 'My.Trans')

}
