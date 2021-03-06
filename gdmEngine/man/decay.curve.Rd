% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{decay.curve}
\alias{decay.curve}
\title{Three parameter logistic function}
\usage{
decay.curve(x, a, b, c)
}
\arguments{
\item{x}{(float) The x-value being assessed}

\item{a}{(float) Parameter specifying the minimum y-value}

\item{b}{(float) Parameter specifying the x-value at curve inflection point}

\item{c}{(float) Parameter specifying the slope of the curve (i.e. steepness of the decay)}
}
\value{
float (the y-value)
}
\description{
A three parameter logistic function for use in sampling [ y = a+(1-a)/(1+((x/b)^(c))) ]
}
\examples{
my.y = decay.curve(x=12, a=0.01, b=10, c=3)
}
