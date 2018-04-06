#'@title GDM builder
#'
#'@description Derive a GDM for a dataset, including automated variable selection based on cross-validation predictive performance. Produces a final GDM object (average parameters across site-pair samples) and associated cross-validation test metrics.
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.
#'@param composition.data (dataframe) A dataframe holding the final species composition data: the records of all species in all grid cells to be used for modelling.
#'@param geo (boolean) Whether geographic distance should be considered as a predictor when deriving the GDM. (default = TRUE)
#'@param train.proportion (float) The proportion of sites in 'site.env.data' to use in training the GDM, with the remaining proportion used to test the model. (default = 0.8)
#'@param n.pairs.train (integer) The number of site-pairs to use in training the GDM. If not specified, the default is to use 10 percent of the total number of pairs possible from the training sites. (default = NULL)
#'@param n.pairs.test (integer) The number of site-pairs to use in testing the GDM. If not specified, the default is to use the same ratio of site-pairs to availabel sites as was used to train the model. (default = NULL)
#'@param n.crossvalid.tests (integer) The number of cross-validation sets to use in deriving the GDM. (default = 10)
#'@param correlation.threshold (float) The maximum correlation (Pearson's R) permitted between candidate predictor variables, as derived from the sites in 'site.env.data'. (default = 0.7)
#'@param selection.metric (string) The model test metric to use in backward elimination of variables for model selection. Options are 'D2' (deviance explained), 'RMSE' (root mean square error), or 'equalised RMSE'. (default = 'D2')
#'@param sample.method (string) The site-pair sample method to use. Options are 'random', 'geodist' (geographic distance), 'envdist' (environmental distance), 'geodens' (geographic density). (default = 'random')
#'@param Indiv.Dev.Explained.Min (float) The minimum amount of deviance explained when potential predictors are assessed individually. (default = 1.0)
#'@param n.predictors.min (integer) The target (or minimum) number of predictor variables in the final model. (default = 8)
#'@param b.used.factor (float) Multiplier for the decay curve parameter (x-value at curve inflection point) for the number of times each site is used in selected pairs. This factor is multiplied by the ratio of n.pairs.target:n.sites in site.env.data to obtain the b.used parameter (default = 2)
#'@param b.dpair.factor (float) For sample method 'geodist'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the distance between sites in a pair. This factor is multiplied by the mean pairwise distance to obtain the b.dpairparameter. (default = 0.5)
#'@param b.epair.factor (float) For sample method 'envdist'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the distance between sites in a pair. This factor is multiplied by the mean environmental distance to obtain the b.epair parameter (default = 1)
#'@param sigma.spair (float) The standard deviation of the isotropic smoothing kernel used in the 'density()' function from the spatstat package, which is used to determine the density of other sites around each site with compositional data. (default = NULL, in which case the value is set to 5% of the total x-axis extent of 'domain.mask')
#'@param spair.factor (float) For sample method 'geodens'. Multiplier for the decay curve parameter (x-value at curve inflection point) for the density of other sites around each site. This factor is multiplied by the mean density of sites to obtain the b.spair parameter. (default = 1.0)
#'@param domain.mask (raster layer) A raster layer specifying the analysis domain
#'@param pcs.projargs (character) A character string of projection arguments; the arguments must be entered exactly as in the PROJ.4 documentation. Used to undertake spatial distance calculations, for example when 'domain.mask' is in geographic coordinate system. For Australian Albers, pcs.projargs="+init=epsg:3577". (default = NULL, in which case the CRS of 'domain.mask' is used for distance calculations).
#'@param bandwidth.geowt (float) The bandwidth to use in the 'geowt' (geographically weighted) sample function. (default = NULL, in which case bandwidth is 5% of the x-axis extent)
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. (default = 'gdm_builder_output')
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return List, including GDM model object and cross validation stats.
#'
#'@examples output = gdm_builder(My.site.env.data, My.composition.data, output.folder = 'C:/Users/processed_data', output.name = 'My.sitepair.data')
#'
#'@importFrom reshape2 dcast
#'@importFrom betapart beta.pair
#'@importFrom gdm gdm
#'@importFrom gdm predict.gdm
#'@importFrom raster pointDistance
#'@importFrom DescTools Gini
#'@importFrom sp CRS
#'
#'@export
gdm_builder <- function(site.env.data, 
                        composition.data,
                        geo=TRUE,
                        train.proportion = 0.8,
                        n.pairs.train = NULL,
                        n.pairs.test = NULL,
                        n.crossvalid.tests = 10,
                        correlation.threshold = 0.7,
                        selection.metric = 'D2',
                        sample.method = 'random',
                        Indiv.Dev.Explained.Min = 1.0,
                        n.predictors.min = 8,
                        b.used.factor=2,
                        b.dpair.factor=0.5,
                        b.epair.factor=1,
                        sigma.spair=NULL,
                        spair.factor=1,
                        domain.mask=NULL,
                        pcs.projargs=NULL, 
                        bandwidth.geowt=NULL,
                        output.folder = NULL,       
                        output.name = "gdm_builder_output",  
                        verbose=TRUE,
                        ...) 
{
  ## NOTE - THIS FUNCTION ESSENTIALLY ASSUMES YOU ARE DOING SOME SITE-PAIR SUBSAMPLING. REQUIRES RE-CODING TO
  ## DEAL WITH CASES WHERE THERE IS NO SITE-PAIR SUBSAMPLING
  start.time <- proc.time()
  
  ## SETUP ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Determine how many sites we are using for the training and testing sets
  n.sites.train <- floor(train.proportion * nrow(site.env.data))
  n.sites.test <- nrow(site.env.data) - n.sites.train
  # If the number of pairs to use in modelling is not specified, use 10% the available pairs as a default
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
  n.vars <- ncol(site.env.data) - n.cols.start 
  # Codify the metric to be used for model selection
  if(!selection.metric=='RMSE'){
    if(selection.metric=='D2'){
      test.col<-4 # deviance explained on test data
      }else{
      test.col<-3 # equalised RMSE
      }# end else
    }else{
    test.col<-2 # RMSE
    } # end if !selection.metric...
  # ensure specis names are a factor in composition.data
  if(!is.factor(composition.data$scientificName))
    {composition.data$scientificName<-as.factor(composition.data$scientificName)}
#### NEW ####
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
#### END NEW ####  
    
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
                                                        domain.mask = domain.mask,
                                                        pcs.projargs = pcs.projargs)    
      Pairs.Table.Test <- sitepair_sample_geo_weighted(site.env.data = Test.Site.Env.Data,
                                                        n.pairs.target = n.pairs.test,
                                                        bandwidth = bandwidth.geowt, 
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
  ## Fit a GDM to each variable independently ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Create the output catcher
  ind.var.test.stats.sampairs <- matrix(0, nrow = (n.vars+1), ncol=4)
  row.names(ind.var.test.stats.sampairs)<- c(colnames(site.env.data[,c((n.cols.start+1):ncol(site.env.data))]),"Geographic")
  colnames(ind.var.test.stats.sampairs) <- c("Mean.Absolute.Error", "Root.Mean.Squre.Error", "Equalised.RMSE","Deviance.Explained")
  # Create a copy of site.env.data and remove all outliers (make them NA),  as they blow up GDM in single-variable fitting
  site.env.data.nout <- site.env.data
  for(i.var in 1:n.vars)
    {
    site.env.data.nout[site.env.data.nout[,(n.cols.start+i.var)] %in% boxplot.stats(site.env.data.nout[,(n.cols.start+i.var)])$out, (n.cols.start+i.var)] <- NA
    }#end for i.var  
    
  # Loop through cross-validation sets  
  for(i.test in 1:n.crossvalid.tests)
    { 
    # Grab the test and train site-pair data 
    Pairs.Table.Train <- train.lst[[i.test]]
    Pairs.Table.Test <- test.lst[[i.test]]
    # loop through each predictor variable
    for(i.var in 1:n.vars)
      {  
      # Catch the env data for both sites in the pair
      # TRAINING
      s1.predictor <- site.env.data.nout[match(as.character(Pairs.Table.Train$s1.site.ID), as.character(site.env.data.nout$xy)), (n.cols.start+i.var)]
      s2.predictor <- site.env.data.nout[match(as.character(Pairs.Table.Train$s2.site.ID), as.character(site.env.data.nout$xy)), (n.cols.start+i.var)]
      Training.table.In <- cbind(Pairs.Table.Train[,c(1:6)], s1.predictor, s2.predictor)
      # TESTING
      s1.predictor <- site.env.data[match(as.character(Pairs.Table.Test$s1.site.ID), as.character(site.env.data.nout$xy)), (n.cols.start+i.var)]
      s2.predictor <- site.env.data[match(as.character(Pairs.Table.Test$s2.site.ID), as.character(site.env.data.nout$xy)), (n.cols.start+i.var)]
      Testing.table.In <- cbind(Pairs.Table.Test[,c(1:6)], s1.predictor, s2.predictor)
      # Run the cross-validation - for strategically sampled test data 
      validation.results<- gdm_SingleCrossValidation(Training.table.In, 
                                                     Testing.table.In)    
      ind.var.test.stats.sampairs[i.var,1] <- as.numeric(ind.var.test.stats.sampairs[i.var,1]) + validation.results$Mean.Absolute.Error
      ind.var.test.stats.sampairs[i.var,2] <- as.numeric(ind.var.test.stats.sampairs[i.var,2]) + validation.results$Root.Mean.Squre.Error
      ind.var.test.stats.sampairs[i.var,3] <- as.numeric(ind.var.test.stats.sampairs[i.var,3]) + validation.results$Equalised.RMSE
      ind.var.test.stats.sampairs[i.var,4] <- as.numeric(ind.var.test.stats.sampairs[i.var,4]) + validation.results$Deviance.Explained
      }# end for i.var
    # Now test geographic distance 
    Training.table.In <- Pairs.Table.Train[,c(1:6)]
    Testing.table.In <- Pairs.Table.Test[,c(1:6)]
    validation.results<- gdm_SingleCrossValidation(Training.table.In, 
                                                   Testing.table.In,
                                                   geo=TRUE)    
    ind.var.test.stats.sampairs[(n.vars+1),1] <- as.numeric(ind.var.test.stats.sampairs[(n.vars+1),1]) + validation.results$Mean.Absolute.Error
    ind.var.test.stats.sampairs[(n.vars+1),2] <- as.numeric(ind.var.test.stats.sampairs[(n.vars+1),2]) + validation.results$Root.Mean.Squre.Error
    ind.var.test.stats.sampairs[(n.vars+1),3] <- as.numeric(ind.var.test.stats.sampairs[(n.vars+1),3]) + validation.results$Equalised.RMSE
    ind.var.test.stats.sampairs[(n.vars+1),4] <- as.numeric(ind.var.test.stats.sampairs[(n.vars+1),4]) + validation.results$Test.Deviance.Explained
    } # end for i.test  
  # Now calculate the means
  ind.var.test.stats.sampairs <- ind.var.test.stats.sampairs / n.crossvalid.tests    
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  proc.time() - ptm
  
  ptm <- proc.time()
  ## Select candidate variables for modelling~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # First assess correlation between variables
  var.cor<-cor(site.env.data[,c((n.cols.start+1):ncol(site.env.data))], use="pairwise.complete.obs",method="pearson")
  # And now combine the model individual variable explanatory power with the level of 
  # corellation between variables to select a set for use in an initial multivariate model.
  # Set up the variable assessment
  var.impt <- ind.var.test.stats.sampairs[,test.col]
  var.impt[ind.var.test.stats.sampairs[,4] < Indiv.Dev.Explained.Min] <- 9999  
  var.impt <- var.impt[-length(var.impt)] #remove geographic distance
  in.var.lst <- c()
  i.var<-1
  while(i.var<n.vars)
    {
    nxt.bst.env<-as.numeric(which.min(var.impt))
    if(as.logical(var.impt[nxt.bst.env]>=9999))
      {break}
    # check if this variable is too correlated to variables already selected
    if(length(in.var.lst)>0)
      {
      cor.to.in.vars<-var.cor[c(in.var.lst),nxt.bst.env]
      if(abs(max(cor.to.in.vars)) > correlation.threshold)
        {
        var.impt[nxt.bst.env] <- 9999
        }# end if(max(cor.to.in.vars) > correlation.threshold)
      if(abs(max(cor.to.in.vars)) <= correlation.threshold)
        {
        in.var.lst<-c(in.var.lst,nxt.bst.env)
        var.impt[nxt.bst.env] <- 9999       
      } # end if(max(cor.to.in.vars) <= correlation.threshold) 
      } # end if length(in.var.lst)>0
    if(length(in.var.lst)<1)
      {
      # catch the variable index and set its error to extreme (9999)
      in.var.lst<-c(in.var.lst,nxt.bst.env)
      var.impt[nxt.bst.env] <- 9999  
    } # end if length(in.var.lst)<1    
    i.var <- i.var + 1
  } # end while i.var<n.env.variables
  # Check if geographic distance should be used too, based on it's individual deviance explained
  if(geo){
    if(ind.var.test.stats.sampairs[nrow(ind.var.test.stats.sampairs),4] < Indiv.Dev.Explained.Min)
      {
      geo<-FALSE
      } # end if(ind.var.test.stats.sampairs...
    }# end if(geo) 
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  proc.time() - ptm
  
  ptm <- proc.time()
##NEW  ## Now run a Backward elimination variable selection routine based ~~~~~~~~~~~~~~~##
  ## on performance under cross-validation. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  in.vars <- in.var.lst
  in.vars.in <- rep(1,times=length(in.vars))
  # get the starting number of predictors
  if(geo){
    n.preds <- length(in.vars) + 1
    }else{
    n.preds <- length(in.vars)
    }
  # Check we have more than the minimum number of predictors to choose from. Otherwise just use all the predictors  
  if(n.preds > n.predictors.min)
    {  
    drop.sequence <-c(n.preds:n.predictors.min)
    }else{
    drop.sequence <- n.preds 
    }
  # Set up GDM input tables for each of the cross-validation sets#####################
  var.names <- c(colnames(site.env.data[,c(n.cols.start + in.vars)])) 
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
      Training.GDM.input.vars[,i.var] <- site.env.data[s1.row.indices.train, (n.cols.start+in.vars[i.var])]
      Training.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[s2.row.indices.train, (n.cols.start+in.vars[i.var])]
      # TESTING
      Testing.GDM.input.vars[,i.var] <- site.env.data[s1.row.indices.test, (n.cols.start+in.vars[i.var])]
      Testing.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[s2.row.indices.test, (n.cols.start+in.vars[i.var])]
      # RANDOM TESTING
      Testing.Rnd.GDM.input.vars[,i.var] <- site.env.data[s1.row.indices.test.rnd, (n.cols.start+in.vars[i.var])]
      Testing.Rnd.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[s2.row.indices.test.rnd, (n.cols.start+in.vars[i.var])]
      } # end for i.var
    # Join the variables to the site-pair table
    Training.table.In <- cbind(Training.table.In, Training.GDM.input.vars)
    Testing.table.In <- cbind(Testing.table.In, Testing.GDM.input.vars)
    Testing.Rnd.table.In <- cbind(Testing.Rnd.table.In, Testing.Rnd.GDM.input.vars)
    ## PUT THE TEST AND TRAINING TABLES IN THE LISTS ## This replaces the first 6 cols with the full GDM tables
    train.name <- paste('PairsTableTrain_',i.test, sep='')
    test.name <- paste('PairsTableTest_',i.test, sep='')
    test.name.rnd <- paste('PairsTableTestRnd_',i.test, sep='')
    train.lst[[train.name]] <- Training.table.In
    test.lst[[test.name]] <- Testing.table.In
    test.lst.rnd[[test.name.rnd]] <- Testing.Rnd.table.In
  }# end for i.test
  ####### DATA IS PREPPED, NOW BACKWARD VAR SELECTION
  # Set up the output catcher for the backward elimination
  if(geo){
    drop.names <- c(var.names,"Geographic")
    }else{
    drop.names <- c(var.names)
    }# end else geo
  drop.stats.D2 <- matrix(nrow=length(drop.names),ncol=length(drop.sequence),dimnames=list(drop.names,paste0(drop.sequence,"_variables")) )
  drop.stats.RMSE <- matrix(nrow=length(drop.names),ncol=length(drop.sequence),dimnames=list(drop.names,paste0(drop.sequence,"_variables")) )
  drop.stats.eRMSE <- matrix(nrow=length(drop.names),ncol=length(drop.sequence),dimnames=list(drop.names,paste0(drop.sequence,"_variables")) )
  # Loop backwards, dropping the worst predictor
  out.col<-1
  for(i.drp in drop.sequence) 
    {
    var.index.vars.in <- which(in.vars.in > 0)
    n.vars.in <- sum(in.vars.in)
    in.vars.cols <- c(1:6, (6+var.index.vars.in), (6+length(in.vars.in)+var.index.vars.in) )
    drop.stats.D2[var.index.vars.in,out.col] <- 0
    drop.stats.RMSE[var.index.vars.in,out.col] <- 0  
    drop.stats.eRMSE[var.index.vars.in,out.col] <- 0
    if(geo)
      {
      drop.stats.D2[nrow(drop.stats.D2),out.col] <- 0
      drop.stats.RMSE[nrow(drop.stats.RMSE),out.col] <- 0
      drop.stats.eRMSE[nrow(drop.stats.eRMSE),out.col] <- 0
      }#end if geo
    # Loop through the cross-validation data sets
    for(i.test in 1:n.crossvalid.tests)  
      {
      # Grab the test and train site-pair data 
      Training.table.In <- train.lst[[i.test]]
      Testing.table.In <- test.lst[[i.test]]
      # Remove previously omitted variables
      Training.table.In <- Training.table.In[,in.vars.cols]
      Testing.table.In <- Testing.table.In[,in.vars.cols] 
      # remove the rows with no data
      Training.table.In<-Training.table.In[complete.cases(Training.table.In),]
      Testing.table.In<-Testing.table.In[complete.cases(Testing.table.In),]
      # Drop each variable, and see how the error changes
      for(i.var in 1:n.vars.in)
        {
        #drop i.var
        p.Training.table.In <- Training.table.In[,-c((6+i.var),(6+n.vars.in+i.var))]
        p.Testing.table.In <- Testing.table.In[,-c((6+i.var),(6+n.vars.in+i.var))]
        # Test the reduced model against the testing data
        validation.results <- gdm_SingleCrossValidation(p.Training.table.In, 
                                                        p.Testing.table.In,
                                                        geo=geo)
        drop.stats.D2[var.index.vars.in[i.var],out.col] <- drop.stats.D2[var.index.vars.in[i.var],out.col] + validation.results$Test.Deviance.Explained
        drop.stats.RMSE[var.index.vars.in[i.var],out.col] <- drop.stats.RMSE[var.index.vars.in[i.var],out.col] + validation.results$Root.Mean.Squre.Error
        drop.stats.eRMSE[var.index.vars.in[i.var],out.col] <- drop.stats.eRMSE[var.index.vars.in[i.var],out.col] + validation.results$Equalised.RMSE
        } #end for i.var
      # If geographic distance is still in, drop it to see the effect
      if(geo)
        {
        validation.results <- gdm_SingleCrossValidation(Training.table.In, 
                                                        Testing.table.In,
                                                        geo=FALSE)
        drop.stats.D2[nrow(drop.stats.D2),out.col] <- drop.stats.D2[nrow(drop.stats.D2),out.col] + validation.results$Test.Deviance.Explained
        drop.stats.RMSE[nrow(drop.stats.RMSE),out.col] <- drop.stats.RMSE[nrow(drop.stats.RMSE),out.col] + validation.results$Root.Mean.Squre.Error
        drop.stats.eRMSE[nrow(drop.stats.eRMSE),out.col] <- drop.stats.eRMSE[nrow(drop.stats.eRMSE),out.col] + validation.results$Equalised.RMSE
        }#end if(geo)
    } # end for each i.test
    # Now average the results across cross-validation sets
    drop.stats.D2[!is.na(drop.stats.D2[,out.col]),out.col] <- drop.stats.D2[!is.na(drop.stats.D2[,out.col]),out.col] / n.crossvalid.tests
    drop.stats.RMSE[!is.na(drop.stats.RMSE[,out.col]),out.col] <- drop.stats.RMSE[!is.na(drop.stats.RMSE[,out.col]),out.col] / n.crossvalid.tests
    drop.stats.eRMSE[!is.na(drop.stats.eRMSE[,out.col]),out.col] <- drop.stats.eRMSE[!is.na(drop.stats.eRMSE[,out.col]),out.col] / n.crossvalid.tests
    # Work out which variable to drop, and drop it
    if(n.preds>n.predictors.min)
      {
      if(test.col == 2)
        {drop.var <- which.min(drop.stats.RMSE[,out.col])}
      if(test.col == 3)
        {drop.var <- which.min(drop.stats.eRMSE[,out.col])}
      if(test.col == 4)
        {drop.var <- which.max(drop.stats.D2[,out.col])}
      if(geo){
        if(drop.var<n.preds){ # then we are dropping an env predictor
          in.vars.in[drop.var] <- 0
          n.preds<-n.preds-1
        }else{ # then we must be dropping geo
          geo<-FALSE
          n.preds<-n.preds-1
        }
      }else{
        in.vars.in[drop.var] <- 0
        n.preds<-n.preds-1
      } # end else
    }# end if n.preds>n.predictors.min
    out.col<-out.col+1
  }# end for i.drp
  ###################################################################################
  ## Now we have a final model, run cross-validation with random test set         ##
  if(n.preds <= n.predictors.min) # -- FINAL -- FINAL -- FINAL -- FINAL --
    {
    # Create a catcher for the final model
    final.mod.MAE.set <- rep(0, times=n.crossvalid.tests)
    final.mod.RMSE.set <- rep(0, times=n.crossvalid.tests)
    final.mod.equRMSE.set <- rep(0, times=n.crossvalid.tests)
    final.mod.D2.set <- rep(0, times=n.crossvalid.tests)
    final.mod.MAE.rnd.set <- rep(0, times=n.crossvalid.tests)
    final.mod.RMSE.rnd.set <- rep(0, times=n.crossvalid.tests)
    final.mod.equRMSE.rnd.set <-rep(0, times=n.crossvalid.tests)
    final.mod.D2.rnd.set <- rep(0, times=n.crossvalid.tests)
    for(i.test in 1:n.crossvalid.tests)  
      {
      # Grab the test and train site-pair data 
      Training.table.In <- train.lst[[i.test]]
      Testing.table.In <- test.lst[[i.test]]      
      Testing.Rnd.table.In <- test.lst.rnd[[i.test]]
      # Remove previously omitted variables
      Training.table.In <- Training.table.In[,in.vars.cols]
      Testing.table.In <- Testing.table.In[,in.vars.cols] 
      Testing.Rnd.table.In <- Testing.Rnd.table.In[,in.vars.cols] 
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
      } # end for each i.test
   } # end if n.preds <= n.predictors.min  
  ###################################################################################
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  proc.time() - ptm
  
  ptm <- proc.time()  
  ## Now we have a final set of predictors, fit a full model sampling site-pairs from
  ## the full set of sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Set up catching matrices for the model parameters
  in.vars <- in.vars[which(in.vars.in > 0)]
  var.names <- var.names[which(in.vars.in > 0)]
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
      Training.GDM.input.vars[,i.var] <- site.env.data[match(as.character(Pairs.Table.Train$s1.site.ID), as.character(site.env.data$xy)), (n.cols.start+in.vars[i.var])]
      Training.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[match(as.character(Pairs.Table.Train$s2.site.ID), as.character(site.env.data$xy)), (n.cols.start+in.vars[i.var])]
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
    } # end for i.test
  #replace the final model object data with the mean values across models
  Final.mod$intercept<-mean(intercept.set)
  Final.mod$explained<-mean(deviance.explained.set)
  Final.mod$coefficients<-colMeans(coefficients.set)
  Final.mod$knots<-colMeans(knots.set)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~## 
  proc.time() - ptm
  r.time <- proc.time() - start.time
  
  ## Now format and return outputs of the model builder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Create an output list
  GDM_Builder_arguments=data.frame(argument=as.character(c('geo',
                                                           'train.proportion',
                                                           'n.pairs.train',
                                                           'n.pairs.test',
                                                           'n.crossvalid.tests',
                                                           'correlation.threshold',
                                                           'selection.metric',
                                                           'sample.method',
                                                           'Indiv.Dev.Explained.Min',
                                                           'n.predictors.min',
                                                           'b.used.factor',
                                                           'b.dpair.factor',
                                                           'b.epair.factor',
                                                           'sigma.spair',
                                                           'b.spair.factor',
                                                           'output.folder',       
                                                           'output.name')),
                                   value=as.character(c(geo,
                                                        train.proportion,
                                                        n.pairs.train,
                                                        n.pairs.test,
                                                        n.crossvalid.tests,
                                                        correlation.threshold,
                                                        selection.metric,
                                                        sample.method,
                                                        Indiv.Dev.Explained.Min,
                                                        n.predictors.min,
                                                        b.used.factor,
                                                        b.dpair.factor,
                                                        b.epair.factor,
                                                        sigma.spair,
                                                        spair.factor,
                                                        output.folder,       
                                                        output.name)))
  
  GDM_Builder_results = list(Inputs = match.call(),
                             Arguments = GDM_Builder_arguments,
                             ProcessingTime = r.time,
                             Backward.Elim.D2 = drop.stats.D2,
                             Backward.Elim.RMSE = drop.stats.RMSE, 
                             Backward.Elim.eRMSE = drop.stats.eRMSE, 
                             Predictors = Final.mod$predictors,
                             Deviance.Explained = deviance.explained.set,
                             Intercept = intercept.set,
                             Dissimilarities.Evenness = final.mod.obs.dissim.evenness,
                             Mean.Absolute.Error = final.mod.MAE.set,
                             Root.Mean.Squre.Error = final.mod.RMSE.set,
                             Equalised.RMSE = final.mod.equRMSE.set,
                             Deviance.Exp = final.mod.D2.set,
                             rnd.Mean.Absolute.Error = final.mod.MAE.rnd.set,
                             rnd.Root.Mean.Squre.Error = final.mod.RMSE.rnd.set,
                             rnd.Equalised.RMSE = final.mod.equRMSE.rnd.set,
                             rnd.Deviance.Exp = final.mod.D2.rnd.set,
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
    save(GDM_Builder_results, file=out.path)
    }#end if !is.null(output.folder) 
  ## Return the output list to the console ##
  return(GDM_Builder_results)
  
} # end gdm_builder() 
    