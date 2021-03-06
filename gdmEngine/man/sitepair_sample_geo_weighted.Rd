% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sitepair_sample_geo_weighted.R
\name{sitepair_sample_geo_weighted}
\alias{sitepair_sample_geo_weighted}
\title{Geographically weighted site-pair sampler}
\usage{
sitepair_sample_geo_weighted(site.env.data, n.pairs.target, bandwidth = NULL,
  b.skip = 3, inter.sample.pt.b.factor = 1, prop.sites.background = 0.1,
  domain.mask, pcs.projargs = NULL, output.folder = NULL,
  output.name = "site_pairs_data_geowt", verbose = FALSE)
}
\arguments{
\item{site.env.data}{(dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.}

\item{n.pairs.target}{(integer) The number of site-pairs to select.}

\item{bandwidth}{(float) The bandwidth to use in the 'geowt' (geographically weighted) sample function. Units are the same x/y units as 'domain.mask' or pcs.projargs (if specified). (default = NULL, in which case bandwidth is 5 percent of the x-axis extent)}

\item{b.skip}{(float) The minimum distance (as a factor of the specified bandwidth) of any data from a geographic 'sample point'. Where all data are greater than 'b.skip' x bandwidth away from a sample point, that sample point will be not used. (default = 3)}

\item{inter.sample.pt.b.factor}{(float) The distance between sample points, as a factor to be multiplied by the bandwidth. (default=1, in which case the distance between sample points is equal to the bandwidth)}

\item{prop.sites.background}{(float) The proportion of sites relative to the geographically weighted sample that will be drawn at random from the whole region (default = 0.1 (i.e. 10 percent))}

\item{domain.mask}{(raster layer) A raster layer specifying the analysis domain}

\item{pcs.projargs}{(character) A character string of projection arguments; the arguments must be entered exactly as in the PROJ.4 documentation. Used to undertake spatial distance calculations, such as when 'domain.mask' is in geographic coordinate system. An example would be specifying Albers projected coordinates for Australia as: pcs.projargs="+init=epsg:3577" . (default = NULL, in which case the CRS of 'domain.mask' is used for distance calculations).}

\item{output.folder}{(string) A folder to save the outputs to. If none specified, no file is written.}

\item{output.name}{(string) A name to use in saving the outputs. Default: 'site_pairs_data_dens'.}

\item{verbose}{(boolean) Print messages to console. Default TRUE.}
}
\value{
Dataframe, site-pairs table, being first 6 columns of the GDM input table, with dissimilarities not calculated.
}
\description{
Sample site-pairs based on the number of nearby sites with composition data.
}
\examples{
output = sitepair_sample_density(My.site.env.data, n.pairs.target=10000, domain.mask=My.mask)

}
