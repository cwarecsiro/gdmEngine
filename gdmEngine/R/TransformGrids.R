#'@title Transform grids
#'
#'@description The R function to call in order to transform big grids. This in turn calls a Rcpp function, that does the magic
#'
#'@param pairs.table (dataframe) A dataframe holding the site-pairs. This is the first six columns of a gdm input table, plus 4 columns hold s1 & s2 long & lat.
#'@param composition.data (dataframe) A dataframe holding the final species composition data: the records of all species in all grid cells to be used for modelling.
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. Default: 'pairs_table_dissim'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return Dataframe.
#'
#'@examples output = calculate_dissimilarities(My.pairs.table, My.composition.data, output.folder = 'C:/Users/processed_data', output.name = 'My.sitepair.data')
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib gdmEngine
#'@export 
TransformGrids <- function(gdm.model,
                           env.grids.stk, 
                           output.folder = NULL,       
                           output.name = "tran",  
                           verbose=TRUE) 
{
  ## Do some checking of the provided inputs
  # Does the raster stack exist?
  if(!class(env.grids.stk) == "RasterStack")
    {stop("A raster stack is required for 'env.grids.stk'") }
  
  # Get the number of layers in the env grid stack
  n.grids <- as.integer(length(env.grids.stk@layers))

  # Get the list of filenames for the env grids
  env.files <- character(length = n.grids)
  for(i.grd in 1:n.grids)
    {env.files[i.grd] <- env.grids.stk@layers[[i.grd]]@file@name}
  # Create a list for the output filenames
  if(gdm.model$geo){
    trans.files <- character(length = (n.grids+2))
  }else{
    trans.files <- character(length = n.grids)
  }
  
  # Check the filepaths (input & output) to make sure they exist
  # Set outfolder to working directory, if not specified
  if(is.null(output.folder))
    {output.folder<-getwd()}
  if(!dir.exists(output.folder))
    {stop("specified output folder doesn't exist")}
  out.filepath <- file.path(output.folder.path,output.name)
  for(i.grd in 1:n.grids)
    {
    if(!file.exists(env.files[i.grd]))
      {stop("cannot find specified env grids on file")}
    if(substr(env.files[i.grd], nchar(env.files[i.grd])-3, nchar(env.files[i.grd])) != ".flt")
      {stop("env grids need to be in .flt format")}
    # Now make the output transformed grid filename
    trans.files[i.grd] <- file.path(output.folder)
    } # end for i.grd    
  
  # And if geographic distance is in the model, make some output filenames for its x & y grids
  # if(gdm.model$geo)
  #   {
  #   trans.files <-
  #   }

 k.one = 1
 k.two = 2
 
 #Call BigGridTransform() from BigGridTransform.cpp
  out.val <- BigGridTransform(k.one,
                              k.two,
                              n.grids,
                              env.files,
                              trans.files)
  
## Irrelevant code below  
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
    out.path <- file.path(output.folder,paste0(output.name,"_",Sys.Date(),".csv")) 
    write.csv(pairs.table.new, out.path, row.names=FALSE)
    # write a log file describing how the data was created *************************************
    fileConn<-file(file.path(output.folder,paste0(output.name,"_",Sys.Date(),"_log_file.txt")),'w')
    writeLines("#######################################################################",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("### Calculate dissimilarities log file ",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines(paste0("### Created ",Sys.time()," using the calculate_dissimilarities() function."),con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    writeLines("",con = fileConn)
    writeLines(paste0("Output data file = ", out.path),con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    close(fileConn) #**************************************************************************
  } # end if !is.null(output.folder)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
  
  # write some feedback to the terminal
  if(verbose)
  {
    msg1 = 'Returned object is a dataframe.'
    msg2 = paste('Dissimilarities have been calculated for ', k.one, ' site-pairs.')
    msg3 = paste(((sum(pairs.table.new$distance < 1)/k.one)*100) , '% of site-pairs have non-complete dissimilarity.')
    if(!is.null(output.folder)){
      msg4 = paste('These data have been also been written to ', out.path)
      cat(paste(msg1, msg2, msg3, msg4, sep = '\n'))      
    }else{
      cat(paste(msg1, msg2, msg3, sep = '\n'))
    }
  }# end if verbose
  
  # Return the freshly filled site-pair table  
  return(out.val)
  
} # end calculate_dissimilarities
  
