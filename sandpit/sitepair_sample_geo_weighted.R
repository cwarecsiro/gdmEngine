#'@title Geographically weighted site-pair sampler
#'
#'@description Sample site-pairs based on the number of nearby sites with composition data. 
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.
#'@param n.pairs.target (integer) The number of site-pairs to select.
#'@param n.pairs.sample (integer) The number of pairs to assess simultaneously as part of the sampling process
#'@param domain.mask (raster layer) A raster layer specifying the analysis domain
#'@param a.used (float) The decay curve parameter (minimum y-value) for the number of times each site is used in selected pairs (default = 0.05)
#'@param b.used.factor (float) Multiplier for the decay curve parameter (x-value at curve inflection point) for the number of times each site is used in selected pairs. This factor is multiplied by the ratio of n.pairs.target:n.sites in site.env.data to obtain the b.used parameter (default = 2)
#'@param c.used (float) The decay curve parameter (slope of the curve) for the number of times each site is used in selected pairs. (default = 3) 
#'@param sigma.spair (float) The standard deviation of the isotropic smoothing kernel used in the 'density()' function from the spatstat package, which is used to determine the density of other sites around each site with compositional data. (default = 0.5)
#'@param a.spair (float) The decay curve parameter (minimum y-value) for the number of sites within the neighbourhood radius (default = 0.05)
#'@param b.spair.factor (float) Multiplier for the decay curve parameter (x-value at curve inflection point) for the density of other sites around each site. This factor is multiplied by the mean density of sites to obtain the b.spair parameter. (default = 1.0
#'@param c.spair (float) The decay curve parameter (slope of the curve) for the number of sites within the neighbourhood radius (default = 3)
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. Default: 'site_pairs_data_dens'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return Dataframe, site-pairs table, being first 6 columns of the GDM input table, with dissimilarities not calculated.
#'
#'@examples output = sitepair_sample_density(My.site.env.data, n.pairs.target=10000, domain.mask=My.mask)
#'
#'@importFrom raster pointDistance
#'@importFrom matrixStats rowMins rowMaxs
#'
#'@export
sitepair_sample_geo_weighted=function(site.env.data,
                                      n.pairs.target,
                                      bandwidth=1, # NEW - bandwidth, in same x/y units as 'domain.mask'
                                      #n.pairs.sample=NULL, 
                                      domain.mask,
                                      #a.used=0.05, 
                                      #b.used.factor=2, 
                                      #c.used=3, 
                                      #sigma.spair=0.5,
                                      #a.spair=0.05, 
                                      #b.spair.factor=1.0, 
                                      #c.spair=1, 
                                      output.folder = NULL,       
                                      output.name = "site_pairs_data_geowt",  
                                      verbose=FALSE)
{
## TO DO:
# - deal properly with lat long position of points vs metres distance
# - add some longer distance samples to the sample of sites around each sample point  
  
  
  # ARGUMENTS ################################
  # specify the distance from a sample point to use in skipping sample points
  b.skip <- 3 # i.e. 3 times the bandwidth
  # specify the distance between sample points, as a multiple of the bandwidth
  inter.sample.pt.b.factor <- 2 # i.e. 2 x bandwidth between sample points
  # specify the proportion of sites to add to the local sample that are drawn at random from the whole region
  prop.sites.background <- 0.1 # i.e. 10% of sites should be a random sample across the whole region
  ###########################################
  
  # First set up the net of 'regresion points' over the spatial grid
  
  d.ij <- 4
  w.ij <- exp(-0.5*((d.ij/bandwidth)^2))
  
  # Work out what the bandwidth is in metres (on average)
  x.d <- (domain.mask@extent@xmax + domain.mask@extent@xmin)/2
  y.d <- (domain.mask@extent@ymax + domain.mask@extent@ymin)/2
  b1.m <- raster::pointDistance(p1=c(x.d,y.d), p2=c((x.d+bandwidth),y.d), lonlat=T)
  b2.m <- raster::pointDistance(p1=c(x.d,y.d), p2=c(x.d,(y.d+bandwidth)), lonlat=T)
  bandwidth.m <- (b1.m + b2.m)/2
    
  ###### Set up the geographic sample points ############################################################
  # start of the geographic sample net
  start.x <- domain.mask@extent@xmin + runif(1,0,(bandwidth/2))
  start.y <- domain.mask@extent@ymin + runif(1,0,(bandwidth/2))  
  # The position of all x's & y's
  n.x <- ceiling((domain.mask@extent@xmax - domain.mask@extent@xmin) / (inter.sample.pt.b.factor*bandwidth))
  n.y <- ceiling((domain.mask@extent@ymax - domain.mask@extent@ymin) / (inter.sample.pt.b.factor*bandwidth))  
  sample.pts <- expand.grid(x.pos = seq(from=start.x, by=(inter.sample.pt.b.factor*bandwidth), length.out = n.x),
                            y.pos = seq(from=start.y, by=(inter.sample.pt.b.factor*bandwidth), length.out = n.y))
  ########################################################################################################

  ###### Select pairs based on weighting at each sample point ############################################
  
  # work out how many random samples to take at each sample point (assuming average weight for any grid cell = 1)
  n.samples <- ceiling((n.pairs.target/nrow(site.env.data))*2)
  # set up a catcher for the sampled pairs
  train.pairs<-NULL
  # Loop through the sample points
  for(i.pt in 1:nrow(sample.pts))
    {
    # determine the distance of each point with composition data from the sample point
    dist.to.pt <- raster::pointDistance(p1=sample.pts[i.pt,], p2=site.env.data[,c(3:4)], lonlat=T) 
    site.wt <- exp(-0.5*((dist.to.pt/bandwidth.m)^2)) 
    # if no sites are  within 3 bandwidths, skip this point
    if(max(site.wt) < (exp(-0.5*(((b.skip*bandwidth.m)/bandwidth.m)^2))) )
      {next}
    in.sites<-NULL
    for(i.sam in 1:n.samples)
      {
      # Use the site weights to randomly select sites
      site.sample <- rbinom(n=length(site.wt), size=1, prob=site.wt) 
      in.sites<-c(in.sites, which(site.sample>0))
      }#end for i.sam
    # If we have sampled more than one site, then process the sites into pairs and add them to the master
    # list of sampled site-pairs
    if(length(in.sites)>1)
      {
      # Sample additional long-distance sites randomly from the full set
      n.background.sites <- floor(length(in.sites) * prop.sites.background)
      if(n.background.sites>0)
        {
        add.sites<-sample.int(nrow(site.env.data), n.background.sites, replace=TRUE) # Note, probably better to sample off the other sample points, as this approach will be biased to heavily sampled areas
        in.sites<-c(in.sites,add.sites)
        } # end if n.background.sites>0
      # permute the order of sites in in.sites
      in.sites<-sample(in.sites, replace = FALSE)
      # Form the sites into pairs (remove the last site if there's an odd number)
      if(!(length(in.sites) %% 2 == 0 ))
        {in.sites<-in.sites[-(length(in.sites))]} 
      n.pairs.sample<-length(in.sites)/2
      ij.pairs<-cbind(in.sites[c(1:n.pairs.sample)], in.sites[c((n.pairs.sample+1):length(in.sites))])
      # rearrange indices for each site pair so that the smalles comes first
      temp.i<-matrixStats::rowMins(ij.pairs)
      temp.j<-matrixStats::rowMaxs(ij.pairs)
      ij.pairs<-cbind(temp.i,temp.j)
      # remove any pairs composed of the same site
      ij.pairs <- ij.pairs[(ij.pairs[,1] != ij.pairs[,2]),]
      # omit duplicate site pairs within this sample (note: we wouldn't do this step under a bootstrapping approach)
      ij.pairs<-unique(ij.pairs)
      # and omit site pairs that have already been selected
      train.pairs<-rbind(train.pairs, ij.pairs)
      train.pairs<-unique(train.pairs)
      }# end if length(in.sites)>1
    }# end for i.pt
  
  ## Processing the selected pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  # Provide a warning if not enough pairs were selected
  if(nrow(train.pairs) < n.pairs.target)
    {
    warning("less pairs selected than desired", call. = FALSE)
    }
  # If we've sampled too many pairs, randomly drop the appropriate amount
  # NOTE - perhaps do this based on densit of sites (dump more where we have more samples)
  if(nrow(train.pairs)>n.pairs.target)
    {train.pairs <- train.pairs[sample.int(nrow(train.pairs), n.pairs.target, replace=FALSE),]}
  # Prepare the start of a GDM input table for the pairs selected
  Pairs.table <- data.frame(distance	= 0,
                            weights = 1,
                            s1.xCoord = site.env.data$decimalLongitude[train.pairs[,1]],
                            s1.yCoord = site.env.data$decimalLatitude[train.pairs[,1]],
                            s2.xCoord = site.env.data$decimalLongitude[train.pairs[,2]],
                            s2.yCoord = site.env.data$decimalLatitude[train.pairs[,2]]) 
  # return the selected pairs
  return(Pairs.table)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  
    

}#end sitepair_sample_geo_weighted()  