#'@title Geographically weighted site-pair sampler
#'
#'@description Sample site-pairs based on the number of nearby sites with composition data. 
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.
#'@param n.pairs.target (integer) The number of site-pairs to select.
#'@param bandwidth (float) The bandwidth to use in the 'geowt' (geographically weighted) sample function. Units are the same x/y units as 'domain.mask' or pcs.projargs (if specified). (default = NULL, in which case bandwidth is 5% of the x-axis extent)
#'@param b.skip (float) The minimum distance (as a factor of the specified bandwidth) of any data from a geographic 'sample point'. Where all data are greater than 'b.skip' x bandwidth away from a sample point, that sample point will be not used. (default = 3)
#'@param inter.sample.pt.b.factor (float) The distance between sample points, as a factor to be multiplied by the bandwidth. (default=1, in which case the distance between sample points is equal to the bandwidth)
#'@param prop.sites.background (float) The proportion of sites relative to the geographically weighted sample that will be drawn at random from the whole region (default = 0.1 (i.e. 10%))
#'@param domain.mask (raster layer) A raster layer specifying the analysis domain
#'@param pcs.projargs (character) A character string of projection arguments; the arguments must be entered exactly as in the PROJ.4 documentation. Used to undertake spatial distance calculations, such as when 'domain.mask' is in geographic coordinate system. An example would be specifying Albers projected coordinates for Australia as: pcs.projargs="+init=epsg:3577" . (default = NULL, in which case the CRS of 'domain.mask' is used for distance calculations).
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. Default: 'site_pairs_data_dens'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return Dataframe, site-pairs table, being first 6 columns of the GDM input table, with dissimilarities not calculated.
#'
#'@examples output = sitepair_sample_density(My.site.env.data, n.pairs.target=10000, domain.mask=My.mask)
#'
#'@importFrom raster projectExtent
#'@importFrom matrixStats rowMins rowMaxs
#'@importFrom sp CRS
#'
#'@export
sitepair_sample_geo_weighted=function(site.env.data,
                                      n.pairs.target,
                                      bandwidth=NULL,              
                                      b.skip=3, 
                                      inter.sample.pt.b.factor=1,  
                                      prop.sites.background=0.1,   
                                      domain.mask,
                                      pcs.projargs=NULL,          
                                      output.folder = NULL,       
                                      output.name = "site_pairs_data_geowt",  
                                      verbose=FALSE)
{
## TO DO:
# - deal properly with lat long position of points vs metres distance
# - add some longer distance samples to the sample of sites around each sample point  
  
  

  ###### Set up the geographic sample points ############################################################
  if(is.null(pcs.projargs))
    {
    # if bandwidth isn't specified, set it at 5% of the grid x-axis extent
    if(is.null(bandwidth))
      {bandwidth <- ((domain.mask@extent@xmax - domain.mask@extent@xmin)/20)}
    # start of the geographic sample net
    start.x <- domain.mask@extent@xmin + runif(1,0,(bandwidth/2))
    start.y <- domain.mask@extent@ymin + runif(1,0,(bandwidth/2))  
    ## The number of x's & y's
    # n.x <- ceiling((domain.mask@extent@xmax - domain.mask@extent@xmin) / (inter.sample.pt.b.factor*bandwidth))
    # n.y <- ceiling((domain.mask@extent@ymax - domain.mask@extent@ymin) / (inter.sample.pt.b.factor*bandwidth))  
    # The number of x's & y's [hexagonal net]
    n.x1 <- ceiling((domain.mask@extent@xmax - start.x) / (inter.sample.pt.b.factor*bandwidth))
    n.x2 <- ceiling((domain.mask@extent@xmax - (start.x + ((inter.sample.pt.b.factor*bandwidth)/2))) / (inter.sample.pt.b.factor*bandwidth))
    y.interval <- sqrt(((inter.sample.pt.b.factor*bandwidth)^2)-(((inter.sample.pt.b.factor*bandwidth)/2)^2))
    n.y <- ceiling((domain.mask@extent@ymax - start.y) / y.interval)
    }else{
    # Find the extent of the domain, in projected coordinates
    pcs.domain.ext <- raster::projectExtent(domain.mask,
                                            crs=sp::CRS(pcs.projargs))
    # if bandwidth isn't specified, set it at 5% of the grid x-axis extent
    if(is.null(bandwidth))
      {bandwidth <- ((pcs.domain.ext@extent@xmax - pcs.domain.ext@extent@xmin)/20)}
    # start of the geographic sample net
    start.x <- pcs.domain.ext@extent@xmin + runif(1,0,(bandwidth/2))
    start.y <- pcs.domain.ext@extent@ymin + runif(1,0,(bandwidth/2))  
    ## The number of x's & y's [square net]
    # n.x <- ceiling((pcs.domain.ext@extent@xmax - pcs.domain.ext@extent@xmin) / (inter.sample.pt.b.factor*bandwidth))
    # n.y <- ceiling((pcs.domain.ext@extent@ymax - pcs.domain.ext@extent@ymin) / (inter.sample.pt.b.factor*bandwidth)) 
    # The number of x's & y's [hexagonal net]
    n.x1 <- ceiling((pcs.domain.ext@extent@xmax - start.x) / (inter.sample.pt.b.factor*bandwidth))
    n.x2 <- ceiling((pcs.domain.ext@extent@xmax - (start.x + ((inter.sample.pt.b.factor*bandwidth)/2))) / (inter.sample.pt.b.factor*bandwidth))
    y.interval <- sqrt(((inter.sample.pt.b.factor*bandwidth)^2)-(((inter.sample.pt.b.factor*bandwidth)/2)^2))
    n.y <- ceiling((pcs.domain.ext@extent@ymax - start.y) / y.interval)
    } # end if... else... is.null(pcs.projargs)
  ## The position of all x's & y's [square net]
  # sample.pts <- expand.grid(x.pos = seq(from=start.x, by=(inter.sample.pt.b.factor*bandwidth), length.out = n.x),
  #                          y.pos = seq(from=start.y, by=(inter.sample.pt.b.factor*bandwidth), length.out = n.y))  
  # The position of all x's & y's [hexagonal net]
  sample.pts.1 <- expand.grid(x.pos = seq(from=start.x, by=(inter.sample.pt.b.factor*bandwidth), length.out = n.x1),
                              y.pos = seq(from=start.y, by=(y.interval*2), length.out = ceiling(n.y/2)))  
  sample.pts.2 <- expand.grid(x.pos = seq(from=(start.x + ((inter.sample.pt.b.factor*bandwidth)/2)), by=(inter.sample.pt.b.factor*bandwidth), length.out = n.x2),
                              y.pos = seq(from=(start.y + y.interval), by=(y.interval*2), length.out = floor(n.y/2)))   
  sample.pts <- rbind(sample.pts.1,sample.pts.2)
  ########################################################################################################  

  ###### Select pairs based on weighting at each sample point ############################################
  # work out how many random samples to take at each sample point (assuming average weight for any grid cell = 1)
  n.samples <- ceiling((n.pairs.target/nrow(site.env.data))*1.5) #multiply by 1.5 for safety buffer
  # set up a catcher for the sampled pairs
  train.pairs<-NULL
  i.sam <- 1
  i.pair.tally <- 0
  while(i.pair.tally < n.pairs.target)
  {
  # Loop through the sample points
  for(i.pt in 1:nrow(sample.pts))
    {
    # determine the distance of each point with composition data from the sample point
    #### NEW ####
    dist.to.pt <- pts.euc.distance(x1=sample.pts[i.pt,1], y1=sample.pts[i.pt,2], x2=site.env.data[,3], y2=site.env.data[,4])
    #### OLD ####    dist.to.pt <- raster::pointDistance(p1=sample.pts[i.pt,], p2=site.env.data[,c(3:4)], lonlat=T) 
    site.wt <- exp(-0.5*((dist.to.pt/bandwidth)^2)) 
    # if no sites are  within 3 bandwidths, skip this point
    if(max(site.wt) < (exp(-0.5*(((b.skip*bandwidth)/bandwidth)^2))) ) #  w.ij <- exp(-0.5*((d.ij/bandwidth)^2))
      {next}
    in.sites<-NULL
  #for(i.sam in 1:n.samples) #try putting the multiple sampling on the outer loop
  #  {
      # Use the site weights to randomly select sites
      site.sample <- rbinom(n=length(site.wt), size=1, prob=site.wt) 
      in.sites<-c(in.sites, which(site.sample>0))
  #  }#end for i.sam
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
  i.pair.tally <- nrow(train.pairs)
  if(i.sam > n.samples)
    {break}
  i.sam <- i.sam + 1  
  }# end while i.pair.tally < n.pairs.target
  
  ## Processing the selected pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  # Provide a warning if not enough pairs were selected
  if(nrow(train.pairs) < n.pairs.target)
    {
    warning("less pairs selected than desired", call. = FALSE)
    }
  # If we've sampled too many pairs, randomly drop the appropriate amount
  # NOTE - perhaps do this based on density of sites (dump more where we have more samples)
  if(nrow(train.pairs)>n.pairs.target)
    {train.pairs <- train.pairs[sample.int(nrow(train.pairs), n.pairs.target, replace=FALSE),]}
  # Prepare the start of a GDM input table for the pairs selected
  Pairs.table <- data.frame(distance	= 0,
                            weights = 1,
                            s1.xCoord = site.env.data$xCoord[train.pairs[,1]],
                            s1.yCoord = site.env.data$yCoord[train.pairs[,1]],
                            s2.xCoord = site.env.data$xCoord[train.pairs[,2]],
                            s2.yCoord = site.env.data$yCoord[train.pairs[,2]],
                            s1.decimalLongitude = site.env.data$decimalLongitude[train.pairs[,1]],
                            s1.decimalLatitude = site.env.data$decimalLatitude[train.pairs[,1]],
                            s2.decimalLongitude = site.env.data$decimalLongitude[train.pairs[,2]],
                            s2.decimalLatitude = site.env.data$decimalLatitude[train.pairs[,2]]) 
  # return the selected pairs
  return(Pairs.table)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  
    

}#end sitepair_sample_geo_weighted()  