#'@title Single GDM cross-validation test
#'
#'@description Undertake a cross-validation assessment of a GDM, using a single training and testing dataset.
#'
#'@param GDM_Table_Training_In (dataframe) A dataframe holding the GDM input table for model fitting.
#'@param GDM_Table_Testing_In (dataframe) A dataframe holding the GDM input table for model testing, having identical column names to 'GDM_Table_Training_In' but using different site-pairs.
#'@param geo (boolean) Geographic distance to be used in model fitting (default = FALSE) .
#'@param splines (vector) An optional vector of the number of I-spline basis functions to be used for each predictor in fitting the model.
#'@param knots (vector) An optional vector of knots in units of the predictor variables to be used in the fitting process. 
#'
#'@returns List, providing cross-validation statistics.
#'
#'@examples output = gdm_SingleCrossValidation(My.training.table, My.testing.table, geo=TRUE)
#'
#'@importFrom gdm gdm
#'@importFrom gdm predict.gdm
#'
#'@export
gdm_SingleCrossValidation=function(GDM_Table_Training_In, 
                                   GDM_Table_Testing_In,
                                   geo=FALSE,
                                   splines=NULL,
                                   knots=NULL)
{
  # omit rows with NAs
  GDM_Table_Training <- GDM_Table_Training_In[complete.cases(GDM_Table_Training_In),]
  GDM_Table_Testing <- GDM_Table_Testing_In[complete.cases(GDM_Table_Testing_In),]
  
  # suppress the warnings GDM often throws up
  oldw <- getOption("warn")
  options(warn = -1)
  
  # Fit the gdm on the training set of sites
  train.mod <- gdm(GDM_Table_Training, 
                   geo=geo,
                   splines=splines,
                   knots=knots)
  # now predict the dissimilarity for the test sites (pairs)
  pred.test <- predict(train.mod,
                       GDM_Table_Testing)
  
  # reset the warnings
  options(warn = oldw)
  
  # calculate an array of test statistics using observed & predicted dissimilarities in the
  # testing data
  # Error
  error = GDM_Table_Testing[,1] - pred.test
  # Mean Error
  ME = mean(error)
  # Mean Absolute Error
  MAE = mean(abs(error))
  # Root Mean Square Error
  RMSE = sqrt(mean(error^2))
  # correlation, obs vs exp
  ObsPred.R <- cor(pred.test , GDM_Table_Testing[,1] , method = "pearson")
  # something new breaking error up over the range of observed values
  obs.breaks<-seq(from=0, to=1, by=0.05)
  obs.value.err <- data.frame("obs.dissim"=seq(from=0.025, to=0.975, by=0.05),
                              "obs.count"=rep(0, times=20),
                              "pred.RMSE"=rep(0, times=20))
  obs.err<-cbind(GDM_Table_Testing[,1], error)
  for(i.cat in 1:nrow(obs.value.err))
  {
    # find the observed values within the specified range of dissimilarities
    if(i.cat == 1){
      cat.obs.err<-obs.err[(obs.err[,1]>=obs.breaks[i.cat] & obs.err[,1]<=obs.breaks[(i.cat+1)]),2]
    }else{
      cat.obs.err<-obs.err[(obs.err[,1]>obs.breaks[i.cat] & obs.err[,1]<=obs.breaks[(i.cat+1)]),2]
    }
    # check if there are any observations in this category
    if(length(cat.obs.err)>0) {
      obs.value.err[i.cat,2]<-length(cat.obs.err)
      obs.value.err[i.cat,3]<-sqrt(mean(cat.obs.err^2))
    }else{
      obs.value.err[i.cat,3]<-NA  
    }
  } #end for i.cat
  # summarise the prediction error across dissimilarity bands
  equ.RMSE<-mean(obs.value.err[,3], na.rm=TRUE)
  # write the outputs of the function
  list(Deviance.Explained = train.mod$explained,
       Mean.Error = ME,
       Mean.Absolute.Error = MAE,
       Root.Mean.Squre.Error = RMSE,
       Obs.Pred.Correlation = ObsPred.R,
       Equalised.RMSE = equ.RMSE,
       Error.by.Observed.Value = obs.value.err)
} # end gdm_SingleCrossValidation function
