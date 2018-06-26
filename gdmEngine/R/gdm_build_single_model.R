#'@title GDM build single model
#'
#'@description Fit a GDM for a dataset, using a specified set of predictor variables, and specified site-pair sampling settings. Produces a final GDM object (average parameters across site-pair samples) and a single sample input table (from one of the site-pair samples).
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data. Sites are rows, Col 1 = 'xy' (character of x & y coordinates joined with underscore between), col2 = ignored (i.e. anything), col3 = decimalLongitude, col4 = decimalLatitude, col5+ = predictor data.
#'@param composition.data (dataframe) A dataframe holding the final species composition data: the records of all species in all grid cells to be used for modelling. Col1 = 'scientificName', col2 = 'decimalLongitude', col3 = decimalLongitude
#'@param predictor.names (character vector) A list of the names of the predictor variables to use in the model. These must match column names in 'site.env.data'.
#'@param geo (boolean) Whether geographic distance should be considered as a predictor when deriving the GDM. (default = TRUE)
#'@param train.proportion (float) The proportion of sites in 'site.env.data' to use in training the GDM, with the remaining proportion used to test the model. (default = 0.8)
#'@param n.pairs.train (integer) The number of site-pairs to use in training the GDM. If not specified, the default is to use 10 percent of the total number of pairs possible from the training sites. (default = NULL)
#'@param n.pairs.test (integer) The number of site-pairs to use in testing the GDM. If not specified, the default is to use the same ratio of site-pairs to availabel sites as was used to train the model. (default = NULL)
#'@param n.crossvalid.tests (integer) The number of cross-validation sets to use in deriving the GDM. (default = 10)
#'@param sample.method (string) The site-pair sample method to use. Options are 'random', 'geodist' (geographic distance), 'envdist' (environmental distance), 'geodens' (geographic density), 'geowt' (geographically weighted). (default = 'random')
#'@param b.used.factor (float) Multiplier for the decay curve parameter (x-value at curve inflection point) for the number of times each site is used in selected pairs. This factor is multiplied by the ratio of n.pairs.target:n.sites in site.env.data to obtain the b.used parameter (default = 2)
#'@param b.dpair.factor (float) For sample method 'geodist'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the distance between sites in a pair. This factor is multiplied by the mean pairwise distance to obtain the b.dpairparameter. (default = 0.5)
#'@param b.epair.factor (float) For sample method 'envdist'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the distance between sites in a pair. This factor is multiplied by the mean environmental distance to obtain the b.epair parameter (default = 1)
#'@param sigma.spair (float) The standard deviation of the isotropic smoothing kernel used in the 'density()' function from the spatstat package, which is used to determine the density of other sites around each site with compositional data. (default = NULL, in which case the value is set to 5 percent of the total x-axis extent of 'domain.mask')
#'@param spair.factor (float) For sample method 'geodens'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the density of other sites around each site. This factor is multiplied by the mean density of sites to obtain the b.spair parameter. (default = 1.0)
#'@param domain.mask (raster layer) A raster layer specifying the analysis domain
#'@param pcs.projargs (character) A character string of projection arguments; the arguments must be entered exactly as in the PROJ.4 documentation. Used to undertake spatial distance calculations, for example when 'domain.mask' is in geographic coordinate system. For Australian Albers, pcs.projargs="+init=epsg:3577". (default = NULL, in which case the CRS of 'domain.mask' is used for distance calculations).
#'@param bandwidth.geowt (float) The bandwidth to use in the 'geowt' (geographically weighted) sample function. (default = NULL, in which case bandwidth is 5 percent of the x-axis extent)
#'@param bandwidth.skip (float) The minimum distance (as a factor of the specified bandwidth) of any data from a geographic 'sample point'. Where all data are greater than 'b.skip' x bandwidth away from a sample point, that sample point will be not used. (default = NULL)
#'@param bandwidth.DistFact (float) The distance between sample points, as a factor to be multiplied by the bandwidth. (default = NULL)
#'@param geowt.RndProp (float) The proportion of sites relative to the geographically weighted sample that will be drawn at random from the whole region (default = NULL)
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. (default = 'gdm_builder_output')
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return List, including GDM model object and cross validation stats.
#'
#'@examples output = gdm_builder(My.site.env.data, My.composition.data, output.folder = 'C:/Users/processed_data', output.name = 'My.sitepair.data')
#'
#'@importFrom gdm gdm
#'@importFrom gdm predict.gdm
#'@importFrom DescTools Gini
#'@importFrom sp CRS SpatialPoints
#'
#'@export
gdm_build_single_model <- function(site.env.data, 
                                   composition.data,
                                   predictor.names,
                                   geo=TRUE,
                                   train.proportion = 0.8,
                                   n.pairs.train = NULL,
                                   n.pairs.test = NULL,
                                   n.crossvalid.tests = 10,
                                   sample.method = 'random',
                                   b.used.factor=2,
                                   b.dpair.factor=0.5,
                                   b.epair.factor=1,
                                   sigma.spair=1,
                                   spair.factor=NULL,
                                   domain.mask=NULL,
                                   pcs.projargs=NULL, 
                                   bandwidth.geowt=NULL,
                                   bandwidth.skip=NULL,
                                   bandwidth.DistFact=NULL,
                                   geowt.RndProp=NULL,
                                   output.folder = NULL,       
                                   output.name = "gdm_single_model_output",  
                                   verbose=TRUE,
                                   ...) 
{
  ## NOTE - THIS FUNCTION ESSENTIALLY ASSUMES YOU ARE DOING SOME SITE-PAIR SUBSAMPLING. REQUIRES RE-CODING TO
  ## DEAL WITH CASES WHERE THERE IS NO SITE-PAIR SUBSAMPLING
  ## To fix: 
  ##         1) find out why models are better in backward elimination than final (? Test this further. Can't see it yet)
  
  
  start.time <- proc.time()
  ## SETUP ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Determine how many sites we are using for the training and testing sets
  n.sites.train <- floor(train.proportion * nrow(site.env.data))
  n.sites.test <- nrow(site.env.data) - n.sites.train
  # If the number of pairs to use in modelling is not specified, use 10 percent the available pairs as a default
  if(is.null(n.pairs.train))
    {
    n.pairs.total <- ((n.sites.train^2)-n.sites.train)/2
    n.pairs.train <- floor(n.pairs.total*0.1)
    }# end if is.null(n.pairs.model)
  # If the number of pairs to use in testing the model is not specified, use the same proportion of sites
  # to pairs as used for the model training
  if(is.null(n.pairs.test))
    {
    pairs.sites.ratio.train <- n.pairs.train/n.sites.train
    n.pairs.total <- ((n.sites.test^2)-n.sites.test)/2
    n.pairs.test <- min(floor(n.sites.test*pairs.sites.ratio.train),n.pairs.total)
    }# end if is.null(n.pairs.model)
  # Establish working parameters for site x env data
  n.cols.start <- 4
  n.vars <- length(predictor.names) 
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
  # add the x & y to site.env.data, including only the specified predictors
  site.env.data <- cbind(site.env.data[,c(1:2)], xCoord, yCoord, site.env.data[,c(3,4)], site.env.data[,which(colnames(site.env.data) %in% predictor.names)])
  n.cols.start <- 6
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~## 
 
  
  ptm <- proc.time()
  ## CREATE CROSS_VALIDATION SETS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  # Store the site-pair data in a list
  train.lst <- list()
  test.lst <- list()
  test.lst.rnd <- list()
  # loop over the number of cross-validation tests, and create the site-pair table (first 6 cols)
  for(i.test in 1:n.crossvalid.tests)
  { 
    ## SPLIT DATA FOR CROSS-VALIDATION (TRAINING AND TESTING SETS) ##
    train.indices <- sample(seq_len(nrow(site.env.data)), size = n.sites.train)
    Train.Site.Env.Data <- site.env.data[train.indices, ]
    Test.Site.Env.Data <- site.env.data[-train.indices, ]
    ## SUBSAMPLE SITE-PAIRS   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  ##
    if(sample.method == 'random')
    {
      Pairs.Table.Train <- sitepair_sample_random(site.env.data = Train.Site.Env.Data,
                                                  n.pairs.target = n.pairs.train)
      Pairs.Table.Test <- sitepair_sample_random(site.env.data = Test.Site.Env.Data,
                                                 n.pairs.target = n.pairs.test)
    }#end if sample.method == 'random'
    if(sample.method == 'geodist')
    {
      Pairs.Table.Train <- sitepair_sample_geographic(site.env.data = Train.Site.Env.Data,
                                                      n.pairs.target = n.pairs.train,
                                                      b.used.factor=b.used.factor,
                                                      b.dpair.factor=b.dpair.factor)
      Pairs.Table.Test <- sitepair_sample_geographic(site.env.data = Test.Site.Env.Data,
                                                     n.pairs.target = n.pairs.test,
                                                     b.used.factor=b.used.factor,
                                                     b.dpair.factor=b.dpair.factor)
    }#end if sample.method == 'geodist'
    if(sample.method == 'envdist')
    {
      Pairs.Table.Train <- sitepair_sample_environment(site.env.data = Train.Site.Env.Data,
                                                       n.pairs.target = n.pairs.train,
                                                       b.used.factor=b.used.factor,
                                                       b.epair.factor=b.epair.factor)
      Pairs.Table.Test <- sitepair_sample_environment(site.env.data = Test.Site.Env.Data,
                                                      n.pairs.target = n.pairs.test,
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
      Pairs.Table.Test <- sitepair_sample_density(site.env.data = Test.Site.Env.Data,
                                                  n.pairs.target = n.pairs.test,
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
      
      Pairs.Table.Test <- sitepair_sample_geo_weighted(site.env.data = Test.Site.Env.Data,
                                                        n.pairs.target = n.pairs.test,
                                                        bandwidth = bandwidth.geowt,
                                                        b.skip = bandwidth.skip, #
                                                        inter.sample.pt.b.factor = bandwidth.DistFact,  #
                                                        prop.sites.background = geowt.RndProp,   #
                                                        domain.mask = domain.mask,
                                                        pcs.projargs = pcs.projargs) 
    }# end if(sample.method == 'geowt')
    # And always have a purely random set of site-pairs for model testing as well
    Pairs.Table.Test.Rnd <- sitepair_sample_random(site.env.data = Test.Site.Env.Data,
                                                   n.pairs.target = n.pairs.test)
    ##  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ##
    ## CALCULATE DISSIMILARITIES ##
    Pairs.Table.Train <- calculate_dissimilarities(pairs.table = Pairs.Table.Train, 
                                                   composition.data = composition.data,
                                                   verbose=FALSE) # Time consuming
    Pairs.Table.Test <- calculate_dissimilarities(pairs.table = Pairs.Table.Test, 
                                                  composition.data = composition.data,
                                                  verbose=FALSE) # Time consuming
    Pairs.Table.Test.Rnd <- calculate_dissimilarities(pairs.table = Pairs.Table.Test.Rnd, 
                                                      composition.data = composition.data,
                                                      verbose=FALSE) # Time consuming
    ## ADD SITE NAMES ##
    Pairs.Table.Train$s1.site.ID <- paste(Pairs.Table.Train$s1.decimalLongitude, Pairs.Table.Train$s1.decimalLatitude, sep = '_')
    Pairs.Table.Train$s2.site.ID <- paste(Pairs.Table.Train$s2.decimalLongitude, Pairs.Table.Train$s2.decimalLatitude, sep = '_')
    Pairs.Table.Test$s1.site.ID <- paste(Pairs.Table.Test$s1.decimalLongitude, Pairs.Table.Test$s1.decimalLatitude, sep = '_')
    Pairs.Table.Test$s2.site.ID <- paste(Pairs.Table.Test$s2.decimalLongitude, Pairs.Table.Test$s2.decimalLatitude, sep = '_')
    Pairs.Table.Test.Rnd$s1.site.ID <- paste(Pairs.Table.Test.Rnd$s1.decimalLongitude, Pairs.Table.Test.Rnd$s1.decimalLatitude, sep = '_')
    Pairs.Table.Test.Rnd$s2.site.ID <- paste(Pairs.Table.Test.Rnd$s2.decimalLongitude, Pairs.Table.Test.Rnd$s2.decimalLatitude, sep = '_')
    ## PUT THE TEST AND TRAINING TABLES IN THE LISTS ##
    train.name <- paste('PairsTableTrain_',i.test, sep='')
    test.name <- paste('PairsTableTest_',i.test, sep='')
    test.name.rnd <- paste('PairsTableTestRnd_',i.test, sep='')
    train.lst[[train.name]] <- Pairs.Table.Train
    test.lst[[test.name]] <- Pairs.Table.Test
    test.lst.rnd[[test.name.rnd]] <- Pairs.Table.Test.Rnd
  } # end for i.test
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  proc.time() - ptm


  ptm <- proc.time()
  ## Assess the performance of the specified model through cross-validation ~~~~~~~~~~~~~~~##
  in.vars <- c(1:length(predictor.names))
  # get the starting number of predictors
  if(geo){
    n.preds <- length(in.vars) + 1
    }else{
    n.preds <- length(in.vars)
    }# END if(geo)
  # Set up GDM input tables for each of the cross-validation sets#####################
  var.names <- c(colnames(site.env.data[,c(n.cols.start + in.vars)])) 
  # Create a catcher for the final model
  final.mod.MAE.set <- rep(0, times=n.crossvalid.tests)
  final.mod.RMSE.set <- rep(0, times=n.crossvalid.tests)
  final.mod.equRMSE.set <- rep(0, times=n.crossvalid.tests)
  final.mod.D2.set <- rep(0, times=n.crossvalid.tests)
  final.mod.MAE.rnd.set <- rep(0, times=n.crossvalid.tests)
  final.mod.RMSE.rnd.set <- rep(0, times=n.crossvalid.tests)
  final.mod.equRMSE.rnd.set <-rep(0, times=n.crossvalid.tests)
  final.mod.D2.rnd.set <- rep(0, times=n.crossvalid.tests)
  ## Loop through the crossvalidation sets
  for(i.test in 1:n.crossvalid.tests)  
    {
    # Grab the test and train site-pair data 
    Pairs.Table.Train <- train.lst[[i.test]]
    Pairs.Table.Test <- test.lst[[i.test]]
    Pairs.Table.Test.Rnd <- test.lst.rnd[[i.test]]
    # Format these dataframes into GDM input tables
    Training.table.In <- Pairs.Table.Train[,c(1:6)]
    Testing.table.In <- Pairs.Table.Test[,c(1:6)]
    Testing.Rnd.table.In <- Pairs.Table.Test.Rnd[,c(1:6)]
    # Prepare predictor variables table
    Training.GDM.input.vars <- matrix(0, nrow=nrow(Training.table.In), ncol=(length(in.vars)*2))
    Testing.GDM.input.vars <- matrix(0, nrow=nrow(Testing.table.In), ncol=(length(in.vars)*2))
    Testing.Rnd.GDM.input.vars <- matrix(0, nrow=nrow(Testing.Rnd.table.In), ncol=(length(in.vars)*2)) 
    colnames(Training.GDM.input.vars) <- c(paste0("s1.",var.names),paste0("s2.",var.names))
    colnames(Testing.GDM.input.vars) <- c(paste0("s1.",var.names),paste0("s2.",var.names))
    colnames(Testing.Rnd.GDM.input.vars) <- c(paste0("s1.",var.names),paste0("s2.",var.names))
    # find the row indices in env.data for the sites in the pairs table
    s1.row.indices.train <- match(as.character(Pairs.Table.Train$s1.site.ID), as.character(site.env.data$xy))
    s2.row.indices.train <- match(as.character(Pairs.Table.Train$s2.site.ID), as.character(site.env.data$xy))    
    s1.row.indices.test <- match(as.character(Pairs.Table.Test$s1.site.ID), as.character(site.env.data$xy))
    s2.row.indices.test <- match(as.character(Pairs.Table.Test$s2.site.ID), as.character(site.env.data$xy))    
    s1.row.indices.test.rnd <- match(as.character(Pairs.Table.Test.Rnd$s1.site.ID), as.character(site.env.data$xy))
    s2.row.indices.test.rnd <- match(as.character(Pairs.Table.Test.Rnd$s2.site.ID), as.character(site.env.data$xy))
    # catch the env data for both sites in the pair
    for(i.var in 1:length(in.vars))
    {
      # TRAINING
      Training.GDM.input.vars[,i.var] <- site.env.data[s1.row.indices.train, (n.cols.start+i.var)]
      Training.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[s2.row.indices.train, (n.cols.start+i.var)]
      # TESTING
      Testing.GDM.input.vars[,i.var] <- site.env.data[s1.row.indices.test, (n.cols.start+i.var)]
      Testing.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[s2.row.indices.test, (n.cols.start+i.var)]
      # RANDOM TESTING
      Testing.Rnd.GDM.input.vars[,i.var] <- site.env.data[s1.row.indices.test.rnd, (n.cols.start+i.var)]
      Testing.Rnd.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[s2.row.indices.test.rnd, (n.cols.start+i.var)]
    } # end for i.var
    # Join the variables to the site-pair table
    Training.table.In <- cbind(Training.table.In, Training.GDM.input.vars)
    Testing.table.In <- cbind(Testing.table.In, Testing.GDM.input.vars)
    Testing.Rnd.table.In <- cbind(Testing.Rnd.table.In, Testing.Rnd.GDM.input.vars)
    # remove the rows with no data
    Training.table.In<-Training.table.In[complete.cases(Training.table.In),]
    Testing.table.In<-Testing.table.In[complete.cases(Testing.table.In),]
    Testing.Rnd.table.In<-Testing.Rnd.table.In[complete.cases(Testing.Rnd.table.In),]
    # For the applied sub-sampling scheme
    validation.results.smp<- gdm_SingleCrossValidation(Training.table.In, 
                                                       Testing.table.In,
                                                       geo=geo)
    final.mod.MAE.set[i.test] <- validation.results.smp$Mean.Absolute.Error
    final.mod.RMSE.set[i.test] <- validation.results.smp$Root.Mean.Squre.Error
    final.mod.equRMSE.set[i.test] <- validation.results.smp$Equalised.RMSE
    final.mod.D2.set[i.test] <- validation.results.smp$Test.Deviance.Explained
    # For a purely random sample
    validation.results.rnd<- gdm_SingleCrossValidation(Training.table.In, 
                                                       Testing.Rnd.table.In,
                                                       geo=geo)
    final.mod.MAE.rnd.set[i.test] <- validation.results.rnd$Mean.Absolute.Error
    final.mod.RMSE.rnd.set[i.test] <- validation.results.rnd$Root.Mean.Squre.Error
    final.mod.equRMSE.rnd.set[i.test] <- validation.results.rnd$Equalised.RMSE
    final.mod.D2.rnd.set[i.test] <- validation.results.rnd$Test.Deviance.Explained
  }# end for i.test
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  proc.time() - ptm

  
    
  ptm <- proc.time()  
  ## Now we have a final set of predictors, fit a full model sampling site-pairs from
  ## the full set of sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Set up catching matrices for the model parameters
  n.params<-n.preds*3 # ASSUMING 3 KNOTS PER VARIABLE AT THE MOMENT
  intercept.set<-rep(0, times=n.crossvalid.tests)
  deviance.explained.set <- rep(0, times=n.crossvalid.tests)
  final.mod.obs.dissim.evenness <- rep(0, times=n.crossvalid.tests)
  coefficients.set<-matrix(0, nrow=n.crossvalid.tests, ncol=n.params)
  knots.set<-matrix(0, nrow=n.crossvalid.tests, ncol=n.params)
  final.mod.dissimilarity<-matrix(0, nrow=n.crossvalid.tests, ncol=6, dimnames=list(paste0("sample",c(1:n.crossvalid.tests)),c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")))
  # Loop through the samples of site pairs, fit GDMs, catch statistics
  for(i.test in 1:n.crossvalid.tests)
  { 
    ## SUBSAMPLE SITE-PAIRS  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ##
    if(sample.method == 'random')
    {
      Pairs.Table.Train <- sitepair_sample_random(site.env.data = site.env.data,
                                                  n.pairs.target = n.pairs.train)
    }#end if sample.method == 'random'
    if(sample.method == 'geodist')
    {
      Pairs.Table.Train <- sitepair_sample_geographic(site.env.data = site.env.data,
                                                      n.pairs.target = n.pairs.train,
                                                      b.used.factor=b.used.factor,
                                                      b.dpair.factor=b.dpair.factor)
    }#end if sample.method == 'geodist'
    if(sample.method == 'envdist')
    {
      Pairs.Table.Train <- sitepair_sample_environment(site.env.data = site.env.data,
                                                       n.pairs.target = n.pairs.train,
                                                       b.used.factor=b.used.factor,
                                                       b.epair.factor=b.epair.factor)
    }#end if sample.method == 'envdist'
    if(sample.method == 'geodens')
    {
      Pairs.Table.Train <- sitepair_sample_density(site.env.data = site.env.data,
                                                   n.pairs.target = n.pairs.train,
                                                   domain.mask=domain.mask,
                                                   pcs.projargs=pcs.projargs,
                                                   b.used.factor=b.used.factor,
                                                   sigma.spair=sigma.spair,
                                                   b.spair.factor=spair.factor)
    }#end if sample.method == 'geodens'
    if(sample.method == 'geowt')
    {
      Pairs.Table.Train <- sitepair_sample_geo_weighted(site.env.data = site.env.data,
                                                        n.pairs.target = n.pairs.train,
                                                        bandwidth = bandwidth.geowt,
                                                        b.skip = bandwidth.skip, #
                                                        inter.sample.pt.b.factor = bandwidth.DistFact,  #
                                                        prop.sites.background = geowt.RndProp,   #
                                                        domain.mask = domain.mask,
                                                        pcs.projargs = pcs.projargs)   
      
    }# end if(sample.method == 'geowt')
    ##  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ##
    Pairs.Table.Train <- calculate_dissimilarities(pairs.table = Pairs.Table.Train, 
                                                   composition.data = composition.data,
                                                   verbose=FALSE) 
    Pairs.Table.Train$s1.site.ID <- paste(Pairs.Table.Train$s1.decimalLongitude, Pairs.Table.Train$s1.decimalLatitude, sep = '_')
    Pairs.Table.Train$s2.site.ID <- paste(Pairs.Table.Train$s2.decimalLongitude, Pairs.Table.Train$s2.decimalLatitude, sep = '_')
    # Format these dataframes into GDM input tables
    Training.table.In <- Pairs.Table.Train[,c(1:6)]
    # Prepare predictor variables table
    Training.GDM.input.vars <- matrix(0, nrow=nrow(Training.table.In), ncol=(length(in.vars)*2))
    colnames(Training.GDM.input.vars) <- c(paste0("s1.",var.names),paste0("s2.",var.names))
    # catch the env data for both sites in the pair
    for(i.var in 1:length(in.vars))
    {
      Training.GDM.input.vars[,i.var] <- site.env.data[match(as.character(Pairs.Table.Train$s1.site.ID), as.character(site.env.data$xy)), (n.cols.start+i.var)]
      Training.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[match(as.character(Pairs.Table.Train$s2.site.ID), as.character(site.env.data$xy)), (n.cols.start+i.var)]
    } # end for i.var
    # Join the variables to the site-pair table
    Training.table.In <- cbind(Training.table.In, Training.GDM.input.vars)
    Training.table.In<-Training.table.In[complete.cases(Training.table.In),]
    # Fit a GDM [without warnings]
    oldw <- getOption("warn")
    options(warn = -1)
    Final.mod <- gdm(Training.table.In,
                     geo=geo)
    options(warn = oldw)
    # Catch the parameters/stats
    intercept.set[i.test] <- Final.mod$intercept
    deviance.explained.set[i.test] <- Final.mod$explained
    coefficients.set[i.test,] <- Final.mod$coefficients
    knots.set[i.test,] <- Final.mod$knots
    final.mod.dissimilarity[i.test,]<-summary(Training.table.In$distance)
    dissim.hist <- hist(Training.table.In$distance, breaks=seq(from=0, to=1, by=0.025),plot=FALSE)
    final.mod.obs.dissim.evenness[i.test] <- 1 - (DescTools::Gini(dissim.hist$counts))

## NEW - RUN THE SIGNIFICANCE TEST ON THIS SAMPLE     
    # gdm_ext_sigtest(dllpath="//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/GDM_EXT_For_Karel/GDM4Rext.dll",
    #                 wdpath = data.processing.folder,
    #                 datatable = "//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/MyGDMinputTable.csv", # GDM input table saved to .csv
    #                 outname = "sig_test",
    #                 iterations = 100,
    #                 do_geo = geo)
## END NEW    
    
  } # end for i.test
  #replace the final model object data with the mean values across models
  Final.mod$intercept<-mean(intercept.set)
  Final.mod$explained<-mean(deviance.explained.set)
  Final.mod$coefficients<-colMeans(coefficients.set)
  Final.mod$knots<-colMeans(knots.set)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~## 
  proc.time() - ptm
  r.time <- proc.time() - start.time
 
#TEMP - WRITE OUT THE LAST GDM INPUT TABLE TO FILE, TO ENABLE RUNNING THE SIGNIFICANCE TEST  
table.path <- file.path(output.folder,paste0(output.name,"_GDMtable_",Sys.Date(),".csv")) 
write.csv(Training.table.In, file=table.path, row.names = FALSE)   
#END TEMP  
  
  ## Now format and return outputs of the model builder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Create an output list
  GDM_Builder_arguments=data.frame(argument=as.character(c('geo',
                                                           'train.proportion',
                                                           'n.pairs.train',
                                                           'n.pairs.test',
                                                           'n.crossvalid.tests',
                                                           'sample.method',
                                                           'b.used.factor',
                                                           'b.dpair.factor',
                                                           'b.epair.factor',
                                                           'sigma.spair',
                                                           'b.spair.factor',
                                                           'domain.mask',
                                                           'pcs.projargs', 
                                                           'bandwidth.geowt',
                                                           'bandwidth.skip',
                                                           'bandwidth.DistFact',
                                                           'geowt.RndProp',
                                                           'output.folder',       
                                                           'output.name')),
                                   value=as.character(rep(NA,times=19)))
  GDM_Builder_arguments$value <- as.character(GDM_Builder_arguments$value)
  GDM_Builder_arguments$value[1] <- as.character(geo)
  GDM_Builder_arguments$value[2] <- as.character(train.proportion)
  GDM_Builder_arguments$value[3] <- as.character(n.pairs.train)
  GDM_Builder_arguments$value[4] <- as.character(n.pairs.test)
  GDM_Builder_arguments$value[5] <- as.character(n.crossvalid.tests)
  GDM_Builder_arguments$value[6] <- as.character(sample.method)
  GDM_Builder_arguments$value[7] <- as.character(b.used.factor)
  GDM_Builder_arguments$value[8] <- as.character(b.dpair.factor)
  GDM_Builder_arguments$value[9] <- as.character(b.epair.factor)
  GDM_Builder_arguments$value[10] <- as.character(sigma.spair)
  if(!is.null(spair.factor)) {GDM_Builder_arguments$value[11] <- as.character(spair.factor)}
  if(!is.null(domain.mask)) {GDM_Builder_arguments$value[12] <- as.character(domain.mask@file@name)}
  if(!is.null(pcs.projargs)) {GDM_Builder_arguments$value[13] <- as.character(pcs.projargs)}
  if(!is.null(bandwidth.geowt)) {GDM_Builder_arguments$value[14] <- as.character(bandwidth.geowt)}
  if(!is.null(bandwidth.skip)) {GDM_Builder_arguments$value[15] <- as.character(bandwidth.skip)}
  if(!is.null(bandwidth.DistFact)) {GDM_Builder_arguments$value[16] <- as.character(bandwidth.DistFact)}
  if(!is.null(geowt.RndProp)) {GDM_Builder_arguments$value[17] <- as.character(geowt.RndProp)}
  if(!is.null(output.folder))GDM_Builder_arguments$value[18] <- as.character(output.folder)       
  GDM_Builder_arguments$value[19] <- as.character(output.name)
  
  
  GDM_Final_Model = list(Inputs = match.call(),
                         Arguments = GDM_Builder_arguments,
                         ProcessingTime = r.time,
                         Deviance.Explained = deviance.explained.set,
                         Intercept = intercept.set,
                         Dissimilarities.Evenness = final.mod.obs.dissim.evenness,
                         Mean.Absolute.Error = final.mod.MAE.set,
                         Root.Mean.Squre.Error = final.mod.RMSE.set,
                         Equalised.RMSE = final.mod.equRMSE.set,
                         Deviance.Exp.CrossVal = final.mod.D2.set,
                         rnd.Mean.Absolute.Error = final.mod.MAE.rnd.set,
                         rnd.Root.Mean.Squre.Error = final.mod.RMSE.rnd.set,
                         rnd.Equalised.RMSE = final.mod.equRMSE.rnd.set,
                         rnd.Deviance.Exp.CrossVal = final.mod.D2.rnd.set,
                         dissimilarity.summary = final.mod.dissimilarity,
                         Mean.Final.GDM = Final.mod)
  
  ## Write the output list to file if specified ##
  if(!is.null(output.folder))
  {
    if(!dir.exists(output.folder))
    {
      dir.create(output.folder)
    }# end if !dir.exists
    out.path <- file.path(output.folder,paste0(output.name,"_",Sys.Date(),".Rdata")) 
    save(GDM_Final_Model, file=out.path)
  }#end if !is.null(output.folder) 
  ## Return the output list to the console ##
  return(GDM_Final_Model)
  
} # end gdm_build_single_model() 
