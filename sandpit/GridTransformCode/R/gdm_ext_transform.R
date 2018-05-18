#'@title Generate GDM transformed grids rapidly and efficiently
#'
#'@description Derive a GDM for a dataset, including automated variable selection based on cross-validation predictive performance. Produces a final GDM object (average parameters across site-pair samples) and associated cross-validation test metrics.
#'
#'@param model A gdm model object created with gdm.fit
#'@param inlist A string containing the path to a .csv file containing the paths to the ESRI Binary Export Grids that will be transformed. Use one path per line and no header in the file list. Grids should be in the same order as the predictors in the gdm model.
#'@param outdir A string containing the path to the output directory to write the transformed grids.
#'@param extrap_type a string containing one of the four following extrapolation methods: "Default10", "WholeGrad", "Conservative", "None".
#'
#'@details gdm_ext.transform utilises the quantiles and coefficients in a gdm model to transform a set of ESRI Binary Export grids. The grids may extend beyond the quantiles defined by the model and these extrapolations may be dealt with using one of the extrap_type methods defined above.
#' "Default10" is the default method and it uses a linear extrapolation of the 0%-10% and 90%-100% slopes to set the extrapolated values.
#' "WholeGrad" uses a linear extrapolation of the 0%-100% slope to set the extrapolated values.
#' "Conservative" makes a choice of the minimum slope at either end of the continuum between the Default10 and WholeGrad methods.
#' "None" uses the minium and maximum transformed values for any gridcells that are to be extrapolated.
#' A grid with an "_ERR" file extenstion is also created detailing any of the grids cells that were extrapolated and also a grid called ABS_ERR_SUM is created with the sums of the extrapolated gridcells. Both these outputs may be useful for diagnostic purposes.
#'
#'@return gdm_ext_transform returns null. 
#'
#'@examples ## 
#'  inpaths <- "d://gridpaths//myInputGrids.csv"
#'  outdir <- "d://outputs"
#'  Gdm_ext_transform(model, inpaths, outdir)
#'  ##
#'
#' @references Ferrier, S., Manion, G., Elith, J. and Richardson, K. (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. Diversity and Distributions 13: 252-264
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib gdmEngine
#'@export
gdm_ext_transform <- function(model, 
                              inlist, 
                              outdir, 
                              extrap_type=c("Default10", "WholeGrad", "Conservative", "None")) 
{
    ##
    ##
    ## model MUST contain the following fields...
    ## geo,preds,splines,quantiles,coefficients
    ##
    ##
    
    # setup the extrapolation type here
    if (missing(extrap_type))
      etype <- 0 #Default10
    else
      switch(extrap_type,
             Default10=etype <- 0,
             WholeGrad=etype <- 1,
             Conservative=etype <- 2,
             None=etype <- 3)
    
    
    # dyn.load(dllpath)
    # 
    # 
    # z <- .C( "CreateTransformsFromR", 
    #          inlist, 
    #          outdir,
    #          as.integer(etype),
    #          as.integer(model$geo),
    #          as.integer(model$preds),
    #          as.integer(model$splines),
    #          as.double(model$quantiles),
    #          as.double(model$coefficients))
    # 
    # dyn.unload(dllpath)
    
    out.value <-  rcpp_CreateTransformsFromR(inlist, #char **ppInputGridTable,
                                             outdir, #char **ppOutputGridsWorkDir,
                                             as.integer(etype), #int *Extrap_Method,
                                             as.integer(model$geo), #bool *DoGeo, 
                                             as.integer(model$preds), #int *NumPredictors,
                                             as.integer(model$splines), #int *pSplines,
                                             as.double(model$quantiles), #double *pQuantiles, 
                                             as.double(model$coefficients) ) #double *pCoeffs)
    
    
    return(out.value)
  } # end gdm_ext_transform()


