#'@title GDM builder
#'
#'@description Derive a GDM for a dataset, including automated variable selection based on cross-validation predictive performance. Produces a final GDM object (average parameters across site-pair samples) and associated cross-validation test metrics.
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.
#'@param composition.data (dataframe) A dataframe holding the final species composition data: the records of all species in all grid cells to be used for modelling.
#'@param geo (boolean) Whether geographic distance should be considered as a predictor when deriving the GDM. (default = TRUE)
#'@param train.proportion (float) The proportion of sites in 'site.env.data' to use in training the GDM, with the remaining proportion used to test the model. (default = 0.8)
#'@param n.pairs.train (integer) The number of site-pairs to use in training the GDM. If not specified, the default is to use 10% of the total number of pairs possible from the training sites. (default = NULL)
#'@param n.pairs.test (integer) The number of site-pairs to use in testing the GDM. If not specified, the default is to use the same ratio of site-pairs to availabel sites as was used to train the model. (default = NULL)
#'@param n.crossvalid.tests (integer) The number of cross-validation sets to use in deriving the GDM. (default = 10)
#'@param correlation.threshold (float) The maximum correlation (Pearson's R) permitted between candidate predictor variables, as derived from the sites in 'site.env.data'. (default = 0.7)
#'@param selection.metric (float) The model test metric to use in backward elimination of variables for model selection. Options are 'RMSE' or 'equalised RMSE'. (default = 'RMSE')
#'@param Indiv.Dev.Explained.Min (float) The minimum amount of deviance explained when potential predictors are assessed individually. (default = 1.0)
#'@param n.predictors.min (integer) The target (or minimum) number of predictor variables in the final model. (default = 8)
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. (default = 'gdm_builder_output')
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@returns List, including GDM model object and cross validation stats.
#'
#'@examples output = gdm_builder(My.site.env.data, My.composition.data, output.folder = 'C:/Users/processed_data', output.name = 'My.sitepair.data')
#'
#'@importFrom reshape2 dcast
#'@importFrom betapart beta.pair
#'@importFrom gdm gdm gdm.predict
#'@importFrom raster pointDistance
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
                        selection.metric = 'RMSE',
                        Indiv.Dev.Explained.Min = 1.0,
                        n.predictors.min = 8,
                        output.folder = NULL,       
                        output.name = "gdm_builder_output",  
                        verbose=TRUE) 
{
  ## NOTE - THIS FUNCTION ESSENTIALLY ASSUMES YOU ARE DOING SOME SITE-PAIR SUBSAMPLING. REQUIRES RE-CODING TO
  ## DEAL WITH CASES WHERE THERE IS NO SITE-PAIR SUBSAMPLING
  
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
    test.col<-3 # equalised RMSE
    }else{
    test.col<-2 # RMSE
    } # end if !selection.metric...
  
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
    ## SUBSAMPLE SITE-PAIRS   *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ##
    # Random .... (or alternative)
    Pairs.Table.Train <- sitepair_sample_random(site.env.data = Train.Site.Env.Data,
                                                n.pairs.target = n.pairs.train)
    Pairs.Table.Test <- sitepair_sample_random(site.env.data = Test.Site.Env.Data,
                                               n.pairs.target = n.pairs.test)
    # And always have a purely random set of site-pairs for model testing as well
    Pairs.Table.Test.Rnd <- sitepair_sample_random(site.env.data = Test.Site.Env.Data,
                                                   n.pairs.target = n.pairs.test)
    ##  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ##
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
    Pairs.Table.Train$s1.site.ID <- paste(Pairs.Table.Train$s1.xCoord, Pairs.Table.Train$s1.yCoord, sep = '_')
    Pairs.Table.Train$s2.site.ID <- paste(Pairs.Table.Train$s2.xCoord, Pairs.Table.Train$s2.yCoord, sep = '_')
    Pairs.Table.Test$s1.site.ID <- paste(Pairs.Table.Test$s1.xCoord, Pairs.Table.Test$s1.yCoord, sep = '_')
    Pairs.Table.Test$s2.site.ID <- paste(Pairs.Table.Test$s2.xCoord, Pairs.Table.Test$s2.yCoord, sep = '_')
    Pairs.Table.Test.Rnd$s1.site.ID <- paste(Pairs.Table.Test.Rnd$s1.xCoord, Pairs.Table.Test.Rnd$s1.yCoord, sep = '_')
    Pairs.Table.Test.Rnd$s2.site.ID <- paste(Pairs.Table.Test.Rnd$s2.xCoord, Pairs.Table.Test.Rnd$s2.yCoord, sep = '_')
    ## PUT THE TEST AND TRAINING TABLES IN THE LISTS ##
    train.name <- paste('PairsTableTrain_',i.test, sep='')
    test.name <- paste('PairsTableTest_',i.test, sep='')
    test.name.rnd <- paste('PairsTableTestRnd_',i.test, sep='')
    train.lst[[train.name]] <- Pairs.Table.Train
    test.lst[[test.name]] <- Pairs.Table.Test
    test.lst.rnd[[test.name.rnd]] <- Pairs.Table.Test.Rnd
    } # end for i.test
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  ## Fit a GDM to each variable independently ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Create the output catcher
  ind.var.test.stats.sampairs <- matrix(0, nrow = (n.vars+1), ncol=4)
  row.names(ind.var.test.stats.sampairs)<- c(colnames(site.env.data[,c((n.cols.start+1):ncol(site.env.data))]),"Geographic")
  colnames(ind.var.test.stats.sampairs) <- c("Mean.Absolute.Error", "Root.Mean.Squre.Error", "Equalised.RMSE","Deviance.Explained")
  # copy this output file for the randomly selected testing pairs [?? Not sure we need to do this at this stage]
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
      s1.predictor <- site.env.data[match(as.character(Pairs.Table.Train$s1.site.ID), as.character(site.env.data$xy)), (n.cols.start+i.var)]
      s2.predictor <- site.env.data[match(as.character(Pairs.Table.Train$s2.site.ID), as.character(site.env.data$xy)), (n.cols.start+i.var)]
      Training.table.In <- cbind(Pairs.Table.Train[,c(1:6)], s1.predictor, s2.predictor)
      # TESTING
      s1.predictor <- site.env.data[match(as.character(Pairs.Table.Test$s1.site.ID), as.character(site.env.data$xy)), (n.cols.start+i.var)]
      s2.predictor <- site.env.data[match(as.character(Pairs.Table.Test$s2.site.ID), as.character(site.env.data$xy)), (n.cols.start+i.var)]
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
    ind.var.test.stats.sampairs[(n.vars+1),4] <- as.numeric(ind.var.test.stats.sampairs[(n.vars+1),4]) + validation.results$Deviance.Explained
    } # end for i.test  
  # Now calculate the means
  ind.var.test.stats.sampairs <- ind.var.test.stats.sampairs / n.crossvalid.tests    
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
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

  ## Now run a Backward elimination variable selection routine based ~~~~~~~~~~~~~~~##
  ## on performance under cross-validation. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  in.vars <- in.var.lst
  final.mod.MAE.set <- rep(0, times=n.crossvalid.tests)
  final.mod.RMSE.set <- rep(0, times=n.crossvalid.tests)
  final.mod.equRMSE.set <- rep(0, times=n.crossvalid.tests)
  final.mod.MAE.rnd.set <- rep(0, times=n.crossvalid.tests)
  final.mod.RMSE.rnd.set <- rep(0, times=n.crossvalid.tests)
  final.mod.equRMSE.rnd.set <-rep(0, times=n.crossvalid.tests)
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
  # Loop backwards, dropping the worst predictor
  #for(i.drp in n.preds:n.predictors.min)
  for(i.drp in drop.sequence) 
    {
    var.names <- c(colnames(site.env.data[,c(n.cols.start + in.vars)])) 
    # Set up the output catcher
    drop.var.test.stats <- matrix(0, nrow=n.preds, ncol=3)
    colnames(drop.var.test.stats) <- c("deviance.explained", "Root.Mean.Squre.Error", "Equalised.RMSE")
    if(geo){
      row.names(drop.var.test.stats)<- c(paste0("drop.",var.names),"drop.Geographic")
      }else{
      row.names(drop.var.test.stats)<- c(paste0("drop.",var.names))
      }# end else geo
    # Loop through the cross-validation data sets
    for(i.test in 1:n.crossvalid.tests)  
      {
      # Grab the test and train site-pair data 
      Pairs.Table.Train <- train.lst[[i.test]]
      Pairs.Table.Test <- test.lst[[i.test]]
      # Format these dataframes into GDM input tables
      Training.table.In <- Pairs.Table.Train[,c(1:6)]
      Testing.table.In <- Pairs.Table.Test[,c(1:6)]
      # Prepare predictor variables table
      Training.GDM.input.vars <- matrix(0, nrow=nrow(Training.table.In), ncol=(length(in.vars)*2))
      Testing.GDM.input.vars <- matrix(0, nrow=nrow(Testing.table.In), ncol=(length(in.vars)*2))
      colnames(Training.GDM.input.vars) <- c(paste0("s1.",var.names),paste0("s2.",var.names))
      colnames(Testing.GDM.input.vars) <- c(paste0("s1.",var.names),paste0("s2.",var.names))
      # catch the env data for both sites in the pair
      for(i.var in 1:length(in.vars))
        {
        # TRAINING
        Training.GDM.input.vars[,i.var] <- site.env.data[match(as.character(Pairs.Table.Train$s1.site.ID), as.character(site.env.data$xy)), (n.cols.start+in.vars[i.var])]
        Training.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[match(as.character(Pairs.Table.Train$s2.site.ID), as.character(site.env.data$xy)), (n.cols.start+in.vars[i.var])]
        # TESTING
        Testing.GDM.input.vars[,i.var] <- site.env.data[match(as.character(Pairs.Table.Test$s1.site.ID), as.character(site.env.data$xy)), (n.cols.start+in.vars[i.var])]
        Testing.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[match(as.character(Pairs.Table.Test$s2.site.ID), as.character(site.env.data$xy)), (n.cols.start+in.vars[i.var])]
        } # end for i.var
      # Join the variables to the site-pair table
      Training.table.In <- cbind(Training.table.In, Training.GDM.input.vars)
      Testing.table.In <- cbind(Testing.table.In, Testing.GDM.input.vars)
      # remove the rows with no data
      Training.table.In<-Training.table.In[complete.cases(Training.table.In),]
      Testing.table.In<-Testing.table.In[complete.cases(Testing.table.In),]
      # Drop each variable, and see how the error changes
      for(i.var in 1:length(in.vars))
        {
        #drop i.var
        p.Training.table.In <- Training.table.In[,-c((6+i.var),(6+length(in.vars)+i.var))]
        p.Testing.table.In <- Testing.table.In[,-c((6+i.var),(6+length(in.vars)+i.var))]
        # Test the reduced model against the testing data
        validation.results <- gdm_SingleCrossValidation(p.Training.table.In, 
                                                          p.Testing.table.In,
                                                          geo=geo)
        drop.var.test.stats[i.var,1] <- drop.var.test.stats[i.var,1] + validation.results$Deviance.Explained
        drop.var.test.stats[i.var,2] <- drop.var.test.stats[i.var,2] + validation.results$Root.Mean.Squre.Error
        drop.var.test.stats[i.var,3] <- drop.var.test.stats[i.var,3] + validation.results$Equalised.RMSE
        } #end for i.var
      # If geographic distance is still in, drop it to see the effect
      if(geo)
        {
        validation.results <- gdm_SingleCrossValidation(Training.table.In, 
                                                          Testing.table.In,
                                                          geo=FALSE)
        drop.var.test.stats[n.preds,1] <- drop.var.test.stats[n.preds,1] + validation.results$Deviance.Explained
        drop.var.test.stats[n.preds,2] <- drop.var.test.stats[n.preds,2] + validation.results$Root.Mean.Squre.Error
        drop.var.test.stats[n.preds,3] <- drop.var.test.stats[n.preds,3] + validation.results$Equalised.RMSE
        }#end if(geo)
      # If this is the final variable set, also fit a model to the full set of variables (for cross-valid stats)
      if(n.preds <= n.predictors.min) # -- FINAL -- FINAL -- FINAL -- FINAL --
        {
        # For the applied sub-sampling scheme
        validation.results.smp<- gdm_SingleCrossValidation(Training.table.In, 
                                                             Testing.table.In,
                                                             geo=geo)
        final.mod.MAE.set[i.test] <- validation.results.smp$Mean.Absolute.Error
        final.mod.RMSE.set[i.test] <- validation.results.smp$Root.Mean.Squre.Error
        final.mod.equRMSE.set[i.test] <- validation.results.smp$Equalised.RMSE
        # For a completely random sample
        Pairs.Table.Test <- test.lst.rnd[[i.test]]
        Testing.table.In <- Pairs.Table.Test[,c(1:6)]
        Testing.GDM.input.vars <- matrix(0, nrow=nrow(Testing.table.In), ncol=(length(in.vars)*2))
        colnames(Testing.GDM.input.vars) <- c(paste0("s1.",var.names),paste0("s2.",var.names))
        for(i.var in 1:length(in.vars))
          {
          Testing.GDM.input.vars[,i.var] <- site.env.data[match(as.character(Pairs.Table.Test$s1.site.ID), as.character(site.env.data$xy)), (n.cols.start+in.vars[i.var])]
          Testing.GDM.input.vars[,(length(in.vars)+i.var)] <- site.env.data[match(as.character(Pairs.Table.Test$s2.site.ID), as.character(site.env.data$xy)), (n.cols.start+in.vars[i.var])]
          } # end for i.var
        Testing.table.In <- cbind(Testing.table.In, Testing.GDM.input.vars)
        Testing.table.In<-Testing.table.In[complete.cases(Testing.table.In),]
        validation.results.rnd<- gdm_SingleCrossValidation(Training.table.In, 
                                                             Testing.table.In,
                                                             geo=geo)
        final.mod.MAE.rnd.set[i.test] <- validation.results.rnd$Mean.Absolute.Error
        final.mod.RMSE.rnd.set[i.test] <- validation.results.rnd$Root.Mean.Squre.Error
        final.mod.equRMSE.rnd.set[i.test] <- validation.results.rnd$Equalised.RMSE
        } #end if length(in.vars) == n.predictors.min # -- FINAL -- FINAL -- FINAL -- 
      } # end for each i.test
    # Now average the results across cross-validation sets
    drop.var.test.stats <- drop.var.test.stats / n.crossvalid.tests
    # Work out which variable to drop, and drop it
    if(n.preds>n.predictors.min)
      {
      if(geo){
        drop.var <- which.min(drop.var.test.stats[,test.col])
        if(drop.var<n.preds){ # then we are dropping an env predictor
          in.vars<-in.vars[-(which.min(drop.var.test.stats[,test.col]))]
          n.preds<-n.preds-1
        }else{ # then we must be dropping geo
          geo<-FALSE
          n.preds<-n.preds-1
          }
      }else{
        in.vars<-in.vars[-(which.min(drop.var.test.stats[,test.col]))]
        n.preds<-n.preds-1
        } # end else
      }# end if n.preds>n.predictors.min
    }# end for i.drp
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  

  ## Now we have a final set of predictors, fit a full model sampling site-pairs from
  ## the full set of sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Set up catching matrices for the model parameters
  n.params<-n.preds*3 # ASSUMING 3 KNOTS PER VARIABLE AT THE MOMENT
  intercept.set<-rep(0, times=n.crossvalid.tests)
  deviance.explained.set <- rep(0, times=n.crossvalid.tests)
  coefficients.set<-matrix(0, nrow=n.crossvalid.tests, ncol=n.params)
  knots.set<-matrix(0, nrow=n.crossvalid.tests, ncol=n.params)
  # Loop through the samples of site pairs, fit GDMs, catch statistics
  for(i.test in 1:n.crossvalid.tests)
    { 
    ## SUBSAMPLE SITE-PAIRS  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ##
    # Random .... (or alternative)
    Pairs.Table.Train <- sitepair_sample_random(site.env.data = site.env.data,
                                                n.pairs.target = n.pairs.train)
    ##  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  ##
    Pairs.Table.Train <- calculate_dissimilarities(pairs.table = Pairs.Table.Train, 
                                                   composition.data = composition.data,
                                                   verbose=FALSE) # Time consuming
    Pairs.Table.Train$s1.site.ID <- paste(Pairs.Table.Train$s1.xCoord, Pairs.Table.Train$s1.yCoord, sep = '_')
    Pairs.Table.Train$s2.site.ID <- paste(Pairs.Table.Train$s2.xCoord, Pairs.Table.Train$s2.yCoord, sep = '_')
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
    } # end for i.test
  #replace the final model object data with the mean values across models
  Final.mod$intercept<-mean(intercept.set)
  Final.mod$explained<-mean(deviance.explained.set)
  Final.mod$coefficients<-colMeans(coefficients.set)
  Final.mod$knots<-colMeans(knots.set)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~## 
  
  ## Now format and return outputs of the model builder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Create an output list
  GDM_Builder_results = list(Inputs = match.call(),
                             Predictors = Final.mod$predictors,
                             Deviance.Explained = deviance.explained.set,
                             Intercept = intercept.set,
                             Mean.Absolute.Error = final.mod.MAE.set,
                             Root.Mean.Squre.Error = final.mod.RMSE.set,
                             Equalised.RMSE = final.mod.equRMSE.set,
                             rnd.Mean.Absolute.Error = final.mod.MAE.rnd.set,
                             rnd.Root.Mean.Squre.Error = final.mod.RMSE.rnd.set,
                             rnd.Equalised.RMSE = final.mod.equRMSE.rnd.set,
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
    