% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdm_build_single_model.R
\name{gdm_build_single_model}
\alias{gdm_build_single_model}
\title{GDM build single model}
\usage{
gdm_build_single_model(site.env.data, composition.data, predictor.names,
  geo = TRUE, train.proportion = 0.8, n.pairs.train = NULL,
  n.pairs.test = NULL, n.crossvalid.tests = 10, sample.method = "random",
  b.used.factor = 2, b.dpair.factor = 0.5, b.epair.factor = 1,
  sigma.spair = 1, spair.factor = NULL, domain.mask = NULL,
  pcs.projargs = NULL, bandwidth.geowt = NULL, bandwidth.skip = NULL,
  bandwidth.DistFact = NULL, geowt.RndProp = NULL, output.folder = NULL,
  output.name = "gdm_single_model_output", verbose = TRUE, ...)
}
\arguments{
\item{site.env.data}{(dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data. Sites are rows, Col 1 = 'xy' (character of x & y coordinates joined with underscore between), col2 = ignored (i.e. anything), col3 = decimalLongitude, col4 = decimalLatitude, col5+ = predictor data.}

\item{composition.data}{(dataframe) A dataframe holding the final species composition data: the records of all species in all grid cells to be used for modelling. Col1 = 'scientificName', col2 = 'decimalLongitude', col3 = decimalLongitude}

\item{predictor.names}{(character vector) A list of the names of the predictor variables to use in the model. These must match column names in 'site.env.data'.}

\item{geo}{(boolean) Whether geographic distance should be considered as a predictor when deriving the GDM. (default = TRUE)}

\item{train.proportion}{(float) The proportion of sites in 'site.env.data' to use in training the GDM, with the remaining proportion used to test the model. (default = 0.8)}

\item{n.pairs.train}{(integer) The number of site-pairs to use in training the GDM. If not specified, the default is to use 10 percent of the total number of pairs possible from the training sites. (default = NULL)}

\item{n.pairs.test}{(integer) The number of site-pairs to use in testing the GDM. If not specified, the default is to use the same ratio of site-pairs to availabel sites as was used to train the model. (default = NULL)}

\item{n.crossvalid.tests}{(integer) The number of cross-validation sets to use in deriving the GDM. (default = 10)}

\item{sample.method}{(string) The site-pair sample method to use. Options are 'random', 'geodist' (geographic distance), 'envdist' (environmental distance), 'geodens' (geographic density), 'geowt' (geographically weighted). (default = 'random')}

\item{b.used.factor}{(float) Multiplier for the decay curve parameter (x-value at curve inflection point) for the number of times each site is used in selected pairs. This factor is multiplied by the ratio of n.pairs.target:n.sites in site.env.data to obtain the b.used parameter (default = 2)}

\item{b.dpair.factor}{(float) For sample method 'geodist'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the distance between sites in a pair. This factor is multiplied by the mean pairwise distance to obtain the b.dpairparameter. (default = 0.5)}

\item{b.epair.factor}{(float) For sample method 'envdist'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the distance between sites in a pair. This factor is multiplied by the mean environmental distance to obtain the b.epair parameter (default = 1)}

\item{sigma.spair}{(float) The standard deviation of the isotropic smoothing kernel used in the 'density()' function from the spatstat package, which is used to determine the density of other sites around each site with compositional data. (default = NULL, in which case the value is set to 5 percent of the total x-axis extent of 'domain.mask')}

\item{spair.factor}{(float) For sample method 'geodens'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the density of other sites around each site. This factor is multiplied by the mean density of sites to obtain the b.spair parameter. (default = 1.0)}

\item{domain.mask}{(raster layer) A raster layer specifying the analysis domain}

\item{pcs.projargs}{(character) A character string of projection arguments; the arguments must be entered exactly as in the PROJ.4 documentation. Used to undertake spatial distance calculations, for example when 'domain.mask' is in geographic coordinate system. For Australian Albers, pcs.projargs="+init=epsg:3577". (default = NULL, in which case the CRS of 'domain.mask' is used for distance calculations).}

\item{bandwidth.geowt}{(float) The bandwidth to use in the 'geowt' (geographically weighted) sample function. (default = NULL, in which case bandwidth is 5 percent of the x-axis extent)}

\item{bandwidth.skip}{(float) The minimum distance (as a factor of the specified bandwidth) of any data from a geographic 'sample point'. Where all data are greater than 'b.skip' x bandwidth away from a sample point, that sample point will be not used. (default = NULL)}

\item{bandwidth.DistFact}{(float) The distance between sample points, as a factor to be multiplied by the bandwidth. (default = NULL)}

\item{geowt.RndProp}{(float) The proportion of sites relative to the geographically weighted sample that will be drawn at random from the whole region (default = NULL)}

\item{output.folder}{(string) A folder to save the outputs to. If none specified, no file is written.}

\item{output.name}{(string) A name to use in saving the outputs. (default = 'gdm_builder_output')}

\item{verbose}{(boolean) Print messages to console. Default TRUE.}
}
\value{
List, including GDM model object and cross validation stats.
}
\description{
Fit a GDM for a dataset, using a specified set of predictor variables, and specified site-pair sampling settings. Produces a final GDM object (average parameters across site-pair samples) and a single sample input table (from one of the site-pair samples).
}
\examples{
output = gdm_builder(My.site.env.data, My.composition.data, output.folder = 'C:/Users/processed_data', output.name = 'My.sitepair.data')

}
