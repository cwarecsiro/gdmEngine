#'@title GDM builder
#'
#'@description Derive a GDM for a dataset, including automated variable selection based on cross-validation predictive performance. Produces 
#'
#'@param site.env.data (dataframe) A dataframe holding the location and environmental conditions for each site (grid cell) for which we have species composition data.
#'@param composition.data (dataframe) A dataframe holding the final species composition data: the records of all species in all grid cells to be used for modelling.
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. Default: 'pairs_table_dissim'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@returns List, including GDM model object and cross validation stats.
#'
#'@examples output = gdm_builder(My.site.env.data, My.composition.data, output.folder = 'C:/Users/processed_data', output.name = 'My.sitepair.data')
#'
#'@importFrom reshape2 dcast
#'@importFrom betapart beta.pair
#'
#'@export
gdm_builder <- function(site.env.data, 
                        composition.data,
                        train.proportion = 0.8,
                        n.pairs.train = NULL,
                        n.pairs.test = NULL,
                        n.crossvalid.tests = 10,
                        output.folder = NULL,       
                        output.name = "gdm_derivation_output",  
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
  
  
  ## CREATE CROSS_VALIDATION SETS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##  
  # Store the site-pair data in a list
  train.lst <- list()
  test.lst <- list()
  test.lst.rnd <- list()
  # loop over the number of cross-validation tests, and create the site-pair table (first 6 cols)
  for(i.test in 1:n.crossvalid.tests)
    { 
    ## SPLIT DATA FOR CROSS-VALIDATION (TRAINING AND TESTING SETS) ---------------------------##
    train.indices <- sample(seq_len(nrow(site.env.data)), size = n.sites.train)
    Train.Site.Env.Data <- site.env.data[train.indices, ]
    Test.Site.Env.Data <- site.env.data[-train.indices, ]
  
    ## SUBSAMPLE SITE-PAIRS ------------------------------------------------------------------##
    # Random .... (or alternative)
    Pairs.Table.Train <- sitepair_sample_random(site.env.data = Train.Site.Env.Data,
                                                n.pairs.target = n.pairs.train)
    Pairs.Table.Test <- sitepair_sample_random(site.env.data = Test.Site.Env.Data,
                                               n.pairs.target = n.pairs.test)
    # And always have a purely random set of site-pairs for model testing as well
    Pairs.Table.Test.Rnd <- sitepair_sample_random(site.env.data = Test.Site.Env.Data,
                                                   n.pairs.target = n.pairs.test)
    
    ## CALCULATE DISSIMILARITIES -------------------------------------------------------------##
    Pairs.Table.Train <- calculate_dissimilarities(pairs.table = Pairs.Table.Train, 
                                                   composition.data = composition.data,
                                                   verbose=FALSE) # Time consuming
    Pairs.Table.Test <- calculate_dissimilarities(pairs.table = Pairs.Table.Test, 
                                                  composition.data = composition.data,
                                                  verbose=FALSE) # Time consuming
    Pairs.Table.Test.Rnd <- calculate_dissimilarities(pairs.table = Pairs.Table.Test.Rnd, 
                                                      composition.data = composition.data,
                                                      verbose=FALSE) # Time consuming
    
    ## ADD SITE NAMES -----------------------------------------------------------------------##
    Pairs.Table.Train$s1.site.ID <- paste(Pairs.Table.Train$s1.xCoord, Pairs.Table.Train$s1.yCoord, sep = '_')
    Pairs.Table.Train$s2.site.ID <- paste(Pairs.Table.Train$s2.xCoord, Pairs.Table.Train$s2.yCoord, sep = '_')
    Pairs.Table.Test$s1.site.ID <- paste(Pairs.Table.Test$s1.xCoord, Pairs.Table.Test$s1.yCoord, sep = '_')
    Pairs.Table.Test$s2.site.ID <- paste(Pairs.Table.Test$s2.xCoord, Pairs.Table.Test$s2.yCoord, sep = '_')
    Pairs.Table.Test.Rnd$s1.site.ID <- paste(Pairs.Table.Test.Rnd$s1.xCoord, Pairs.Table.Test.Rnd$s1.yCoord, sep = '_')
    Pairs.Table.Test.Rnd$s2.site.ID <- paste(Pairs.Table.Test.Rnd$s2.xCoord, Pairs.Table.Test.Rnd$s2.yCoord, sep = '_')

    ## PUT THE TEST AND TRAINING TABLES IN THE LISTS -----------------------------------------##
    train.name <- paste('PairsTableTrain_',i.test, sep='')
    test.name <- paste('PairsTableTest_',i.test, sep='')
    test.name.rnd <- paste('PairsTableTestRnd_',i.test, sep='')
    
    train.lst[[train.name]] <- Pairs.Table.Train
    test.lst[[test.name]] <- Pairs.Table.Test
    test.lst.rnd[[test.name.rnd]] <- Pairs.Table.Test.Rnd
    } # end for i.test
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  ## Fit a GDM to each variable independently ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  # Establish working parameters
  n.cols.start <- 4
  n.vars <- ncol(site.env.data) - n.cols.start 
  # Create the output catcher
  ind.var.test.stats.sampairs <- matrix(0, nrow = n.vars, ncol=3)
  row.names(ind.var.test.stats.sampairs)<- c(colnames(site.env.data[,c((n.cols.start+1):ncol(site.env.data))]))
  colnames(ind.var.test.stats.sampairs) <- c("Mean.Absolute.Error", "Root.Mean.Squre.Error", "Equalised.RMSE")
  # copy this output file for the randomly selected testing pairs [?? Not sure we need to do this at this stage]
#  ind.var.test.stats.randpairs <- ind.var.test.stats.sampairs
  # Loop through cross-validation sets  
  for(i.test in 1:n.crossvalid.tests)
    { 
    # Grab the test and train site-pair data 
    Pairs.Table.Train <- train.lst[[i.test]]
    Pairs.Table.Test <- test.lst[[i.test]]
#    Pairs.Table.Test.Rnd <- test.lst.rnd[[i.test]]
    # loop through each predictor variable
    for(i.var in 1:n.env.variables)
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
      validation.results<- gdm.SingleCrossValidation(Training.table.In, 
                                                     Testing.table.In,
                                                     geo=FALSE,
                                                     splines=NULL,
                                                     knots=NULL)    
      ind.var.test.stats.sampairs[i.var,1] <- as.numeric(ind.var.test.stats.sampairs[i.var,1]) + validation.results$Mean.Absolute.Error
      ind.var.test.stats.sampairs[i.var,2] <- as.numeric(ind.var.test.stats.sampairs[i.var,2]) + validation.results$Root.Mean.Squre.Error
      ind.var.test.stats.sampairs[i.var,3] <- as.numeric(ind.var.test.stats.sampairs[i.var,3]) + validation.results$Equalised.RMSE
      }# end for i.var
    } # end for i.test  
  # Now calculate the means
  ind.var.test.stats.sampairs <- ind.var.test.stats.sampairs / n.crossvalid.tests    
  
  
  
  
    
} # end gdm_builder() 
    