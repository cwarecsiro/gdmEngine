#'@title Fit geographic distance to GDM residuals
#'
#'@description Fits geographic distance to the residuals of a GDM that has not used geographic distance as a predictor.
#'
#'@param GDM_ModelObject_NoDist (list) The list (model object) that is prduced when fitting a GDM using the gdm package.
#'@param GDM_InputTable_NoDist (dataframe) A dataframe holding the GDM input table used to create 'GDM_ModelObject_NoDist'.
#'@param lonlat (boolean) whether the x/y coordinates in 'GDM_InputTable_NoDist' are longitude/latitude. If True, these are converted to km. default = TRUE.
#'
#'@returns A GDM model object following fitting of geographic distance to residuals.
#'
#'@examples output = gdm_DistanceResidualFit(My.GDM, My.GDM.table)
#'
#'@importFrom raster pointDistance
#'
#'@export
gdm_DistanceResidualFit <- function(GDM_ModelObject_NoDist,
                                    GDM_InputTable_NoDist,
                                    lonlat=TRUE,
                                    )
  {
  # Determine the geographic distance between site-pairs  
  if(lonlat){
    XYDist = pointDistance(GDM_InputTable_NoDist[,3:4], GDM_InputTable_NoDist[,5:6],lonlat=T)/1000 
    } else {
    XYDist = pointDistance(GDM_InputTable_NoDist[,3:4], GDM_InputTable_NoDist[,5:6],lonlat=F)
    }
  
  # Prepare the data for fitting the distance
  dataTab = data.frame("Match" = GDM_ModelObject_NoDist$observed,
                       "Ecological" = GDM_ModelObject_NoDist$ecological,
                       "EcologicalINT" = GDM_ModelObject_NoDist$ecological+GDM_ModelObject_NoDist$intercept,
                       "Distance" = XYDist,
                       "fitted" = GDM_ModelObject_NoDist$predicted)
  
  # The first function is not-crucial, but can help to provide starting
  # coefficients for the next function to optimise. [ & suppress warnings]
  oldw <- getOption("warn")
  options(warn = -1)
  fit.A = glm(Match ~ Distance,
              family = binomial(link=negexp()),
              data = dataTab,
              method = 'nnls_fit')
  fit.B = glm(Match ~ offset(Ecological) + Distance,
              family = binomial(link=negexp()),
              data = dataTab,
              etastart = dataTab$Ecological,
              mustart = dataTab$fitted,
              start = coefficients(fit.A),
              control = list(maxit=1000),
              method = 'nnls_fit')
  options(warn = oldw)  
 

  
  
  
   
  } # end gdm_DistanceResidualFit()  







#summary(fit.B)
dist.coef = fit.B$coefficients[2]
# Thus this coefficient can then be added in when making predictions later on.
print( paste0("Dissimilarity increases by ",round(dist.coef*100,3),
              " for every 100km a pair is seperated.") )
## [1] "Dissimilarity increases by 0.022 for every 100km a pair is seperated."
print(paste("Spatial GDM model fitted. Model explains approximately",
            round(spatial.GDM$explained,1)," % of the deviance.",sep=""))
## [1] "Spatial GDM model fitted. Model explains approximately30.9 % of the deviance."
