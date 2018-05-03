#'@title Sitepair sample assessor
#'
#'@description Assess sitepairs sampled for GDM model fitting. For the selected sample method, a variety of useful stats are generated.
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.
#'@param composition.data (dataframe) A dataframe holding the final species composition data: the records of all species in all grid cells to be used for modelling.
#'@param n.pairs.train (integer) The number of site-pairs to use in training the GDM. If not specified, the default is to use 10 percent of the total number of pairs possible from the training sites. (default = NULL)
#'@param n.pairs.per.site (integer) If 'n.pairs.train' is not specified, this parameter specifies the target of the average number of pairs each site is included in, to calculate n.pairs.train. (default = NULL) 
#'@param prop.pairs.train (float) If 'n.pairs.train' is not specified, and 'n.pairs.per.site' is not specified, a proportion of all possible sitepairs is specified by this argument and used. (default = 0.05 (i.e. 5% of possible number of sitepairs))
#'@param n.crossvalid.tests (integer) The number of cross-validation sets to use in sampling site-pairs. (default = 10)
#'@param sample.method (string) The site-pair sample method to use. Options are 'random', 'geodist' (geographic distance), 'envdist' (environmental distance), 'geodens' (geographic density), 'geowt' (geographically weighted). (default = 'random')
#'@param b.used.factor (float) Multiplier for the decay curve parameter (x-value at curve inflection point) for the number of times each site is used in selected pairs. This factor is multiplied by the ratio of n.pairs.target:n.sites in site.env.data to obtain the b.used parameter (default = 2)
#'@param b.dpair.factor (float) For sample method 'geodist'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the distance between sites in a pair. This factor is multiplied by the mean pairwise distance to obtain the b.dpairparameter. (default = 0.5)
#'@param b.epair.factor (float) For sample method 'envdist'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the distance between sites in a pair. This factor is multiplied by the mean environmental distance to obtain the b.epair parameter (default = 1)
#'@param sigma.spair (float) The standard deviation of the isotropic smoothing kernel used in the 'density()' function from the spatstat package, which is used to determine the density of other sites around each site with compositional data. (default = NULL, in which case the value is set to 5% of the total x-axis extent of 'domain.mask')
#'@param spair.factor (float) For sample method 'geodens'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the density of other sites around each site. This factor is multiplied by the mean density of sites to obtain the b.spair parameter. (default = 1.0)
#'@param domain.mask (raster layer) A raster layer specifying the analysis domain
#'@param pcs.projargs (character) A character string of projection arguments; the arguments must be entered exactly as in the PROJ.4 documentation. Used to undertake spatial distance calculations, for example when 'domain.mask' is in geographic coordinate system. For Australian Albers, pcs.projargs="+init=epsg:3577". (default = NULL, in which case the CRS of 'domain.mask' is used for distance calculations).
#'@param bandwidth.geowt (float) The bandwidth to use in the 'geowt' (geographically weighted) sample function. (default = NULL, in which case bandwidth is 5% of the x-axis extent)
#'@param bandwidth.skip (float) The minimum distance (as a factor of the specified bandwidth) of any data from a geographic 'sample point'. Where all data are greater than 'b.skip' x bandwidth away from a sample point, that sample point will be not used. (default = NULL)
#'@param bandwidth.DistFact (float) The distance between sample points, as a factor to be multiplied by the bandwidth. (default = NULL)
#'@param geowt.RndProp (float) The proportion of sites relative to the geographically weighted sample that will be drawn at random from the whole region (default = NULL)
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. (default = 'gdm_builder_output')
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return List, including useful assessment stats.
#'
#'@examples output = sitepair_sample_assessor(site.env.data=Site.Env.Data, composition.data=Selected.records, n.pairs.train=100000, domain.mask=Aus.domain.mask, pcs.projargs="+init=epsg:3577" ,output.folder = getwd(), output.name = 'MySitepairAssessment')
#'
#'@importFrom DescTools Gini
#'@importFrom sp CRS SpatialPoints spTransform
#'@importFrom raster projectExtent
#'
#'@export
sitepair_sample_assessor <- function(site.env.data, 
                                     composition.data,
                                     n.pairs.train = NULL,
                                     n.pairs.per.site = NULL,
                                     prop.pairs.train = 0.05,
                                     n.crossvalid.tests = 10,
                                     sample.method = 'random',
                                     b.used.factor=2,
                                     b.dpair.factor=0.5,
                                     b.epair.factor=1,
                                     sigma.spair=NULL,
                                     spair.factor=1,
                                     domain.mask=NULL,
                                     pcs.projargs=NULL, 
                                     bandwidth.geowt=NULL,
                                     bandwidth.skip=NULL,
                                     bandwidth.DistFact=NULL,
                                     geowt.RndProp=NULL,
                                     output.folder = NULL,       
                                     output.name = "Sitepair_sample_assessor_output",  
                                     verbose=TRUE,
                                     ...) 
{
  start.time <- proc.time()
  ## SETUP ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Establish working parameters for site x env data
  n.cols.start <- 4
  n.vars <- ncol(site.env.data) - n.cols.start 
  # Determine how many sites we are using for the training and testing sets
  n.sites.train <- nrow(site.env.data)
  # If the number of pairs to use in modelling is not specified, use 'prop.pairs.train' of the available pairs (10% as a default)
  if(is.null(n.pairs.train))
    {
    n.pairs.total <- ((n.sites.train^2)-n.sites.train)/2
    if(!is.null(n.pairs.per.site))
      {
      n.pairs.ideal <- 0.5 * n.pairs.per.site * n.sites.train
      n.pairs.train <- min(n.pairs.ideal, (n.pairs.total*0.5)) # maxing the possible pairs at half the full selection   
      }else{
      n.pairs.train <- floor(n.pairs.total*prop.pairs.train)
      } # end else
    }# end if is.null(n.pairs.model)
  # ensure specis names are a factor in composition.data
  if(!is.factor(composition.data$scientificName))
  {composition.data$scientificName<-as.factor(composition.data$scientificName)}
  # If we're doing spatial calculations in projected coordinate system, convert x & y
  # from site.env.data to pcs & add to site.env.data. Otherwise, just duplicate x & y.
  if(is.null(pcs.projargs))
  {
    xCoord <- site.env.data$decimalLongitude
    yCoord <- site.env.data$decimalLatitude
  }else{
    # create a spatial points object from long / lat
    site.pts.gcs <- sp::SpatialPoints(coords=site.env.data[,c(3,4)],
                                      proj4string=domain.mask@crs)
    # project the spatial points into the specified pcs  
    site.pts.pcs <- sp::spTransform(site.pts.gcs,
                                    CRSobj=sp::CRS(pcs.projargs))  
    xCoord <- site.pts.pcs$decimalLongitude
    yCoord <- site.pts.pcs$decimalLatitude  
  } # end if... else... is.null(pcs.projargs)
  # add the x & y to site.env.data
  site.env.data <- cbind(site.env.data[,c(1:2)], xCoord, yCoord, site.env.data[,c(3:ncol(site.env.data))])
  n.cols.start <- 6
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  
  ## CREATE CROSS_VALIDATION SETS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  # Store the site-pair data in a list
  train.lst <- list()
  # loop over the number of cross-validation tests, and create the site-pair table (first 6 cols)
  for(i.test in 1:n.crossvalid.tests)
  { 
    ## SPLIT DATA FOR CROSS-VALIDATION (TRAINING AND TESTING SETS) ##
    Train.Site.Env.Data <- site.env.data
    ## SUBSAMPLE SITE-PAIRS   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ##
    if(sample.method == 'random')
    {
      Pairs.Table.Train <- sitepair_sample_random(site.env.data = Train.Site.Env.Data,
                                                  n.pairs.target = n.pairs.train)
    }#end if sample.method == 'random'
    if(sample.method == 'geodist')
    {
      Pairs.Table.Train <- sitepair_sample_geographic(site.env.data = Train.Site.Env.Data,
                                                      n.pairs.target = n.pairs.train,
                                                      b.used.factor=b.used.factor,
                                                      b.dpair.factor=b.dpair.factor)
    }#end if sample.method == 'geodist'
    if(sample.method == 'envdist')
    {
      Pairs.Table.Train <- sitepair_sample_environment(site.env.data = Train.Site.Env.Data,
                                                       n.pairs.target = n.pairs.train,
                                                       b.used.factor=b.used.factor,
                                                       b.epair.factor=b.epair.factor)
    }#end if sample.method == 'envdist'
    if(sample.method == 'geodens')
    {
      Pairs.Table.Train <- sitepair_sample_density(site.env.data = Train.Site.Env.Data,
                                                   n.pairs.target = n.pairs.train,
                                                   domain.mask=domain.mask,
                                                   pcs.projargs=pcs.projargs,
                                                   b.used.factor=b.used.factor,
                                                   sigma.spair=sigma.spair,
                                                   b.spair.factor=spair.factor)
      
    }#end if sample.method == 'geodens'
    if(sample.method == 'geowt')
      {
      Pairs.Table.Train <- sitepair_sample_geo_weighted(site.env.data = Train.Site.Env.Data,
                                                        n.pairs.target = n.pairs.train,
                                                        bandwidth = bandwidth.geowt,
                                                        b.skip = bandwidth.skip, #
                                                        inter.sample.pt.b.factor = bandwidth.DistFact,  #
                                                        prop.sites.background = geowt.RndProp,   #
                                                        domain.mask = domain.mask,
                                                        pcs.projargs = pcs.projargs)    
      
      }# end if(sample.method == 'geowt')
    ##  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ##
    ## CALCULATE DISSIMILARITIES ##
    Pairs.Table.Train <- calculate_dissimilarities(pairs.table = Pairs.Table.Train, 
                                                   composition.data = composition.data,
                                                   verbose=FALSE) 
    ## ADD SITE NAMES ##
    Pairs.Table.Train$s1.site.ID <- paste(Pairs.Table.Train$s1.decimalLongitude, Pairs.Table.Train$s1.decimalLatitude, sep = '_')
    Pairs.Table.Train$s2.site.ID <- paste(Pairs.Table.Train$s2.decimalLongitude, Pairs.Table.Train$s2.decimalLatitude, sep = '_')
    ## PUT THE RAINING TABLES IN THE LISTS ##
    train.name <- paste('PairsTableTrain_',i.test, sep='')
    train.lst[[train.name]] <- Pairs.Table.Train
  } # end for i.test
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  ## ASSESS THE SAMPLED SITE-PAIRS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Prepare to catch the outputs of the assessment
  dissim.evenness<-rep(0,times=n.crossvalid.tests)
  dissim.summary<-matrix(0,nrow=n.crossvalid.tests, ncol=6)
  colnames(dissim.summary)<-c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
  sites.geo.evenness<-rep(0,times=n.crossvalid.tests)
  sitepairs.geo.evenness<-rep(0,times=n.crossvalid.tests)
  pairs.geo.distance<-dissim.summary
  sites.env.evenness<-rep(0,times=n.crossvalid.tests)
  sitepairs.env.evenness<-rep(0,times=n.crossvalid.tests)
  pairs.env.distance<-dissim.summary
  sites.ntimes.used<-dissim.summary
  
  # Loop through the samples, and catch the relevant attributes of the site-pairs
  for(i.test in 1:n.crossvalid.tests)  
    { 
  ## Grab the test and train site-pair data 
    Pairs.Table.Train <- train.lst[[i.test]]
  
  ## Dissimilarities ##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dissim.hist <- hist(Pairs.Table.Train$distance,
                        breaks = seq(from=0, to=1, by=0.025),
                        plot=TRUE)
    dissim.evenness[i.test] <- 1 - DescTools::Gini(dissim.hist$counts)#quantifies inequality, so high values are less even, zero = perfect evenness
    dissim.summary[i.test,] <- summary(Pairs.Table.Train$distance)

  ## Geographic distribution (of sites) ##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # break the grid up into 10km boxes and quantify evenness in the distribution of sites used in pairs   
    if(is.null(pcs.projargs))
      {
      bandwidth<-0.1 #degrees
      start.x <- domain.mask@extent@xmin
      start.y <- domain.mask@extent@ymin
      n.x <- ceiling((domain.mask@extent@xmax - start.x) / (bandwidth))
      n.y <- ceiling((domain.mask@extent@ymax - start.y) / (bandwidth))
    }else{
      bandwidth <- 10000 #metres
      pcs.domain.ext <- raster::projectExtent(domain.mask,
                                              crs=sp::CRS(pcs.projargs))
      start.x <- pcs.domain.ext@extent@xmin
      start.y <- pcs.domain.ext@extent@ymin
      n.x <- ceiling((pcs.domain.ext@extent@xmax - pcs.domain.ext@extent@xmin) / (bandwidth))
      n.y <- ceiling((pcs.domain.ext@extent@ymax - pcs.domain.ext@extent@ymin) / (bandwidth)) 
    } # end if... else... is.null(pcs.projargs)
  # Create a coarse spatial grid
  coarse.grid<-matrix(0,nrow=n.y,ncol=n.x)
  # count how many sites in each coarse grid cell
  for(i.site in 1:nrow(site.env.data))
    {
    this.col <- ceiling((site.env.data$xCoord[i.site] - start.x) / bandwidth) 
    this.row <- ceiling((site.env.data$yCoord[i.site] - start.y) / bandwidth) 
    coarse.grid[this.row, this.col] <- coarse.grid[this.row, this.col] + 1
    }# end for i.site
  # turn zeros to NA's & create a new grid from this for the site pairs
  coarse.grid[coarse.grid==0] <- NA
  coarse.grid.pairs <- coarse.grid
  coarse.grid.pairs[!(is.na(coarse.grid.pairs))] <- 0
  # count how many sites used in site pairs are in each coarse grid cell
  for(i.pair in 1:nrow(Pairs.Table.Train))
    {
    # site one in the pair
    this.col <- ceiling((Pairs.Table.Train$s1.xCoord[i.pair] - start.x) / bandwidth) 
    this.row <- ceiling((Pairs.Table.Train$s1.yCoord[i.pair] - start.y) / bandwidth) 
    coarse.grid.pairs[this.row, this.col] <- coarse.grid.pairs[this.row, this.col] + 1
    # site two in the pair
    this.col <- ceiling((Pairs.Table.Train$s2.xCoord[i.pair] - start.x) / bandwidth) 
    this.row <- ceiling((Pairs.Table.Train$s2.yCoord[i.pair] - start.y) / bandwidth) 
    coarse.grid.pairs[this.row, this.col] <- coarse.grid.pairs[this.row, this.col] + 1
    }# end for i.site
  # Now calculate the evenness of distrubution of sites used over the spatial grid
  coarse.grid.vec <- as.vector(coarse.grid)
  coarse.grid.vec <- coarse.grid.vec[!(is.na(coarse.grid.vec))]
  coarse.grid.pairs.vec <- as.vector(coarse.grid.pairs)
  coarse.grid.pairs.vec <- coarse.grid.pairs.vec[!(is.na(coarse.grid.pairs.vec))]  
  sites.geo.evenness[i.test] <- 1 - DescTools::Gini(coarse.grid.vec)# Gini quantifies inequality, so high values are less even, zero = perfect evenness
  sitepairs.geo.evenness[i.test] <- 1 - DescTools::Gini(coarse.grid.pairs.vec)# Gini quantifies inequality, so high values are less even, zero = perfect evenness
      
  ## Geographic distances ##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    pairs.distance <- pts.euc.distance(x1=Pairs.Table.Train$s1.xCoord, y1=Pairs.Table.Train$s1.yCoord, x2=Pairs.Table.Train$s2.xCoord, y2=Pairs.Table.Train$s2.yCoord)
    pairs.geo.distance[i.test,] <- summary(pairs.distance)
    
  ## Environmental distribution (of sites) ##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    env.colnames=c('PTA','TXX')
    # Work out which cols in site.env.data will be used to determine env distance
    oldw <- getOption("warn")
    options(warn = -1)
    env.cols <- which(colnames(site.env.data) == env.colnames)
    options(warn = oldw)
    std.site.env.vars<-as.matrix(site.env.data[,c(3,4,env.cols)])
    #std.site.env.vars[,c(3,4)]<-scale(as.matrix(std.site.env.vars[,c(3,4)])) # Don't think we need to standardise
    n.x <- 100
    n.y <- 100
    start.x <- min(std.site.env.vars[,3]) - 0.0001 # take a marginal amount off start.x to ensure min is within bounds
    start.y <- min(std.site.env.vars[,4]) - 0.0001 # take a marginal amount off start.y to ensure min is within bounds
    bandwidth.x <- (max(std.site.env.vars[,3]) + 0.0001 - start.x) / n.x 
    bandwidth.y <- (max(std.site.env.vars[,4]) + 0.0001 - start.y) / n.y
    # Create a coarse environmental grid
    coarse.grid<-matrix(0,nrow=n.y,ncol=n.x)
    # count how many sites in each coarse grid cell
    for(i.site in 1:nrow(site.env.data))
      {
      this.col <- ceiling((std.site.env.vars[i.site,3] - start.x) / bandwidth.x) 
      this.row <- ceiling((std.site.env.vars[i.site,4] - start.y) / bandwidth.y) 
      coarse.grid[this.row, this.col] <- coarse.grid[this.row, this.col] + 1
      }# end for i.site
    # to plot: 
    #     library(plotly)
    #     plot_ly(z=coarse.grid,type='heatmap')    
    coarse.grid[coarse.grid==0] <- NA
    coarse.grid.pairs <- coarse.grid
    coarse.grid.pairs[!(is.na(coarse.grid.pairs))] <- 0
    # count how many sites used in site pairs are in each coarse grid cell of env space
    # add the env data to each site in the pair
    s1.predictor1 <- site.env.data[match(as.character(Pairs.Table.Train$s1.site.ID), as.character(site.env.data$xy)), (env.cols[1])]
    s2.predictor1 <- site.env.data[match(as.character(Pairs.Table.Train$s2.site.ID), as.character(site.env.data$xy)), (env.cols[1])]
    s1.predictor2 <- site.env.data[match(as.character(Pairs.Table.Train$s1.site.ID), as.character(site.env.data$xy)), (env.cols[2])]
    s2.predictor2 <- site.env.data[match(as.character(Pairs.Table.Train$s2.site.ID), as.character(site.env.data$xy)), (env.cols[2])]
    # link the data for the same predictors together
    predictor1<-c(s1.predictor1, s2.predictor1)  
    predictor2<-c(s1.predictor2, s2.predictor2)    
    # Now loop through the data, aggregating the counts for each part of env space
    for(i.pair in 1:length(predictor1))
      {
      this.col <- ceiling((predictor1[i.pair] - start.x) / bandwidth.x) 
      this.row <- ceiling((predictor2[i.pair] - start.y) / bandwidth.y) 
      coarse.grid.pairs[this.row, this.col] <- coarse.grid.pairs[this.row, this.col] + 1
      }# end for i.pair
    # to plot: 
    #     library(plotly)
    #     plot_ly(z=coarse.grid.pairs,type='heatmap')
    # Now calculate the evenness of distrubution of sites used over the spatial grid
    coarse.grid.vec <- as.vector(coarse.grid)
    coarse.grid.vec <- coarse.grid.vec[!(is.na(coarse.grid.vec))]
    coarse.grid.pairs.vec <- as.vector(coarse.grid.pairs)
    coarse.grid.pairs.vec <- coarse.grid.pairs.vec[!(is.na(coarse.grid.pairs.vec))]  
    sites.env.evenness[i.test] <- 1 - DescTools::Gini(coarse.grid.vec)# Gini quantifies inequality, so high values are less even, zero = perfect evenness
    sitepairs.env.evenness[i.test] <- 1 - DescTools::Gini(coarse.grid.pairs.vec)# Gini quantifies inequality, so high values are less even, zero = perfect evenness
  
  ## Environmental distances ##  
    # Generate a new matrix with standardised values (mean = 0, variance=1) for the selected environment variables
    # Note that this uses the '' data generated above, and changes it through standardisation
    std.site.env.vars<-site.env.data[,c(1,env.cols)]
    std.site.env.vars[,c(2,3)] <- scale(as.matrix(std.site.env.vars[,c(2,3)]))
    # extract the row indices in "std.site.env.vars" for each site in each pair in the pairs table
    s1.row <- match(as.character(Pairs.Table.Train$s1.site.ID), as.character(std.site.env.vars$xy))
    s2.row <- match(as.character(Pairs.Table.Train$s2.site.ID), as.character(std.site.env.vars$xy))
    ij.pairs <- cbind(s1.row,s2.row)
    std.site.env.vars <- as.matrix(std.site.env.vars[,c(2,3)])
    pairs.distance <- PairsDist(std.site.env.vars, ij.pairs)
    pairs.env.distance[i.test,] <- summary(pairs.distance)

  
  ## Number of times used
    # Make a list of all the sites used, and find out how many times for each
    sites.used <- c(Pairs.Table.Train$s1.site.ID, Pairs.Table.Train$s2.site.ID)
    sites.times.used <- table(sites.used)
    sites.times.used <- as.data.frame(sites.times.used)
    sites.times.used$sites.used <- as.character(sites.times.used$sites.used)
    # and add any sites that weren't used, with zeros
    if(nrow(sites.times.used) < nrow(site.env.data))
      {
      sites.used.names <- sites.times.used$sites.used
      out.sites <- as.character(site.env.data$xy[!(as.character(site.env.data$xy) %in% sites.times.used$sites.used)]) 
      out.sites <- data.frame('sites.used' = as.character(out.sites),
                              'Freq' = rep(0, times=length(out.sites)))
      out.sites$sites.used <- as.character(out.sites$sites.used)
      sites.times.used <- rbind(sites.times.used, out.sites)
      } # end if length(sites.times.used) < nrow(site.env.data)
    # Now summarise the number of times sites have been used
    sites.ntimes.used[i.test,]<-summary(sites.times.used$Freq)
    
  } # end for each i.test
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  r.time <- proc.time() - start.time
  
  
  # Create an output list
  Sitepair_sample_assessor_args=data.frame('argument'=as.character(c('n.pairs.train',
                                                           'n.crossvalid.tests',
                                                           'sample.method',
                                                           'b.used.factor',
                                                           'b.dpair.factor',
                                                           'b.epair.factor',
                                                           'sigma.spair',
                                                           'b.spair.factor',
                                                           'pcs.projargs', 
                                                           'bandwidth.geowt',
                                                           'output.folder',       
                                                           'output.name')),
                                     'value'=c(as.character(rep(NA,times=12))))
  Sitepair_sample_assessor_args$value <- as.character(Sitepair_sample_assessor_args$value)
  if(!is.null(n.pairs.train)) {Sitepair_sample_assessor_args$value[1] <- as.character(n.pairs.train)}
  Sitepair_sample_assessor_args$value[2] <- as.character(n.crossvalid.tests)
  Sitepair_sample_assessor_args$value[3] <- as.character(sample.method)
  Sitepair_sample_assessor_args$value[4] <- as.character(b.used.factor)
  Sitepair_sample_assessor_args$value[5] <- as.character(b.dpair.factor)
  Sitepair_sample_assessor_args$value[6] <- as.character(b.epair.factor)
  if(!is.null(sigma.spair)) {Sitepair_sample_assessor_args$value[7] <- as.character(sigma.spair)}
  Sitepair_sample_assessor_args$value[8] <- as.character(spair.factor)
  if(!is.null(pcs.projargs)) {Sitepair_sample_assessor_args$value[9] <- as.character(pcs.projargs)} 
  if(!is.null(bandwidth.geowt)) {Sitepair_sample_assessor_args$value[10] <- as.character(bandwidth.geowt)}
  if(!is.null(output.folder)) {Sitepair_sample_assessor_args$value[11] <- as.character(output.folder)}      
  Sitepair_sample_assessor_args$value[12] <- as.character(output.name)
  
  Sitepair_assessor_results = list(Inputs = match.call(),
                             Arguments = Sitepair_sample_assessor_args,
                             ProcessingTime = r.time,
                             DissimilarityEvenness = dissim.evenness,
                             DissimilaritySummary = dissim.summary,
                             SitesGeoEvenness = sites.geo.evenness,
                             SitepairsGeoEvenness = sitepairs.geo.evenness,
                             SitepairsGeoDistanceSummary = pairs.geo.distance,
                             SitesEnvEvenness = sites.env.evenness,
                             SitepairsEnvEvenness = sitepairs.env.evenness,
                             SitepairsEnvDistanceSummary = pairs.env.distance,
                             nTimesSitesUsedInPairs = sites.ntimes.used)
  
  ## Write the output list to file if specified ##
  if(!is.null(output.folder))
  {
    if(!dir.exists(output.folder))
    {
      dir.create(output.folder)
    }# end if !dir.exists
    out.path <- file.path(output.folder,paste0(output.name,"_",Sys.Date(),".Rdata")) 
    save(Sitepair_assessor_results, file=out.path)
  }#end if !is.null(output.folder) 
  ## Return the output list to the console ##
  return(Sitepair_assessor_results)
  
} # end sitepair_sample_assessor()  
