#'@title Transform grids
#'
#'@description The R function to call in order to transform big grids. This in turn calls a Rcpp function, that does the magic
#'
#'@param gdm.model (gdm model object) A GDM model object for which transformed grids of each predictor will be generated.
#'@param env.grids.stk (raster stack) A raster stack holding grids for all the predictors to be transformed.
#'@param extrap.method (string) The method used to derive transformed values when the raw values are outside the range used to train the model. Options are: "Clamped" (transformed values clamped at the extreme value in the training data); "End10" (transformed values extrapolated based on the slope of the spline function in the final (or initial) 10% range of values used to train the model), "WholeGrad" (transformed values extrapolated based on the slope of the spline function across the whole gradient of values used to train the model), "Conservative" (the lowest spline slope out of End10 and WholeGrad), (default = "Conservative").
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. Default: 'pairs_table_dissim'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return Transformed grids are written to the specified output filepaths.
#'
#'@examples TransformGrids(My.model, My.env.stk, output.folder = 'C:/Users/processed_data', output.name = 'My.Trans')
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib gdmEngine
#'@export 
TransformGrids <- function(gdm.model,
                           env.grids.stk,
                           extrap.method = 'Conservative',
                           output.folder = NULL,       
                           output.name = "tran",  
                           verbose=TRUE) 
{
  ## Do some checking of the provided inputs #######################################
  # Does the raster stack exist?
  if(!class(env.grids.stk) == "RasterStack")
    {stop("A raster stack is required for 'env.grids.stk'") }
  if(!file.exists(output.folder))
    {stop("The specified output.folder doesn't exist") }
  
  ##################################################################################  
  
  ## Get the number of env predictors in the GDM
  npreds <- length(gdm.model$predictors)
  if(gdm.model$geo){
    geo.in <- as.integer(1)
    npreds <- length(gdm.model$predictors)-1
    pred.names <- gdm.model$predictors[-1]
    trans.files <- character(length = (npreds+2))
  }else{
    geo.in <- as.integer(0)
    npreds <- length(gdm.model$predictors)
    pred.names <- gdm.model$predictors
    trans.files <- character(length = npreds)
  } # end if gdm.model$geo
  
  ## Get the filepaths for the predictors and the trans grids
  env.files <- character(length = npreds)
  for(i.prd in 1:npreds)
    {
    layer.index <- which(names(env.grids.stk) == pred.names[i.prd])
    if(length(layer.index)>0){
      env.files[i.prd] <- env.grids.stk@layers[[layer.index]]@file@name
      if(substr(env.files[i.prd], nchar(env.files[i.prd])-3, nchar(env.files[i.prd])) != ".flt")
        {stop("env grids need to be in .flt format")}
      if(!file.exists(env.files[i.prd]))
        {stop(paste0("Cannot find the .flt file for ",pred.names[i.prd]))}
      trans.files[i.prd] <- file.path(output.folder, paste0(pred.names[i.prd], "_", output.name, ".flt") )
      }else{
      stop(paste0("Raster layer for ",pred.names[i.prd], "not found in env.grids.stk"))
      }
    }#end for i.prd
  
  ## And name the geo trans grids if we're creating them
  if(gdm.model$geo)
    {
    trans.files[(npreds+1)] <- file.path(output.folder, paste0("GeoX_", output.name, ".flt") )
    trans.files[(npreds+2)] <- file.path(output.folder, paste0("GeoY_", output.name, ".flt") )
    }
  
  ## establish which extrapolation method to use
    extrap.code <- 1 # assume the default ("Conservative"), even if misspecified
  if(extrap.method == "End10")
    {extrap.code <- 2}
  if(extrap.method == "WholeGrad")
    {extrap.code <- 3}
  if(extrap.method == "Clamped")
    {extrap.code <- 0}
  extrap.code <- as.integer(extrap.code)
  
  ## Finally, set up some key parameters
  n.splines.total <- as.integer(sum(gdm.model$splines))
  grid.rows <- as.integer(nrow(env.grids.stk))
  grid.cols <- as.integer(ncol(env.grids.stk))
  npreds <- as.integer(npreds)
  vec.splines <- as.integer(gdm.model$splines)
  vec.knots <- gdm.model$knots
  vec.coeffs <- gdm.model$coefficients
 
  out.val <- 1
 ## Call BigGridTransform() from BigGridTransform.cpp
  out.val <- BigGridTransform(grid.rows,
                              grid.cols,
                              npreds,
                              geo.in,
                              n.splines.total,
                              vec.splines,
                              vec.knots,
                              vec.coeffs,
                              extrap.code,
                              env.files,
                              trans.files)

  if(out.val < 1)
    {stop("Grid transformation failed  :( ")} 

  ##########################################################################  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
  # Now write out the data to file (if specified) and return the aggregated records
  # write the data to file, if an output folder is specified
  if(!is.null(output.folder))
  {
    if(!dir.exists(output.folder))
    {
      dir.create(output.folder)
    }# end if !dir.exists
    # write a log file describing how the data was created *************************************
    fileConn<-file(file.path(output.folder,paste0(output.name,"_",Sys.Date(),"_log_file.txt")),'w')
    writeLines("#######################################################################",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("### Transform Grids log file ",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines(paste0("### Created ",Sys.time()," using the TransformGrids() function."),con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    writeLines("",con = fileConn)
    writeLines(paste0("Output folder = ", output.folder),con = fileConn)
    writeLines("",con = fileConn)
    writeLines(paste0("--Model Details--"),con = fileConn)
    writeLines("",con = fileConn)
    writeLines(paste0("GDM created: ",gdm.model$creationdate),con = fileConn)
    writeLines(paste0("GDM deviance explained: ",gdm.model$explained, " %"),con = fileConn)
    writeLines(paste0("Predictors: "),con = fileConn)
    writeLines(paste0(gdm.model$predictors),con = fileConn)
    writeLines(paste0("Number of splines: "),con = fileConn)
    writeLines(paste0(gdm.model$splines),con = fileConn)
    writeLines(paste0("Knots: "),con = fileConn)
    writeLines(paste0(gdm.model$knots),con = fileConn)
    writeLines(paste0("Coefficients: "),con = fileConn)
    writeLines(paste0(gdm.model$coefficients),con = fileConn)
    writeLines(paste0("Input environment grids: "),con = fileConn)
    writeLines(paste0(env.files),con = fileConn)
    writeLines(paste0("Extrapolation method: ", extrap.method),con = fileConn)
    writeLines(paste0("Output transformed grids: "),con = fileConn)
    writeLines(paste0(trans.files),con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    close(fileConn) #**************************************************************************
  } # end if !is.null(output.folder)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
  
  # write some feedback to the terminal
  if(verbose)
    {
    msg1 = paste0('Grids have been transformed and written to ', output.folder)
    cat(msg1)      
    }# end if verbose
  

} # end calculate_dissimilarities
  
