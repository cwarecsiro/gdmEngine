#'@title Extract environmental data
#'
#'@description Extracts environmental data for specified locations, from a set of specified spatial grids. Returns a dataframe (site x env). 
#'
#'@param ALA.composition.data (dataframe) A dataframe holding the final species composition data.
#'@param environment.stk (raster stack) A stack of rasters from which environment data will be extracted.
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string), A name to use in saving the outputs. Default: 'filtered_data'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return Dataframe.
#'
#'@examples output = extract_env_data(My.composition.data, My.env.stk, output.folder = 'C:/Users/processed_data', output.name = My.filtered.ALA.data, domain.mask = Aust.ras)
#'
#'@importFrom raster extract
#'
#'@export
extract_env_data = function(ALA.composition.data,             
                            environment.stk,
                            output.folder = NULL,       
                            output.name = "site_env_data",  
                            verbose=TRUE)
{
  # Create a grid-cell name 'xy'
  ALA.composition.data$xy <- paste(ALA.composition.data$decimalLongitude, ALA.composition.data$decimalLatitude, sep = '_')
  
  # Calculate species richness
  agg.check <- count(ALA.composition.data, c('xy', 'scientificName'))
  agg.check$freq <- ifelse(agg.check$freq > 1, 1, agg.check$freq)
  agg.check <- as.data.frame(table(agg.check[,c(1,3)]))
  agg.check <- agg.check[,-2]
  agg.check$xy <- as.character(agg.check$xy)
  colnames(agg.check) [2] <- "richness"

  # add the x/y values back into the working dataframe
  ll <- strsplit(agg.check[,1], '_')
  agg.check$decimalLongitude <- unlist(ll)[ c(TRUE,FALSE) ]
  agg.check$decimalLatitude <- unlist(ll)[ c(FALSE,TRUE) ]
  agg.check$decimalLongitude <- as.numeric(agg.check$decimalLongitude)
  agg.check$decimalLatitude <- as.numeric(agg.check$decimalLatitude)  
  
  # extract the environment data for the cells with records
  env.dat <- raster::extract(environment.stk, SpatialPoints(cbind(agg.check$decimalLongitude,agg.check$decimalLatitude)))
  
  # bind the environment data back
  agg.check <- cbind(agg.check, as.data.frame(env.dat))

  # Now write out the data to file (if specified) and return the aggregated records
  # write the data to file, if an output folder is specified
  if(!is.null(output.folder))
  {
    if(!dir.exists(output.folder))
    {
      dir.create(output.folder)
    }# end if !dir.exists
    out.path <- file.path(output.folder,paste0(output.name,"_",Sys.Date(),".csv")) 
    write.csv(agg.check, out.path, row.names=FALSE)
    # write a log file describing how the data was created *************************************
    fileConn<-file(file.path(output.folder,paste0(output.name,"_",Sys.Date(),"_log_file.txt")),'w')
    writeLines("#######################################################################",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("### Extract environment data log file ",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines(paste0("### Created ",Sys.time()," using the extract_env_data() function."),con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    writeLines("",con = fileConn)
    writeLines(paste0("Output data file = ", out.path),con = fileConn)
    writeLines(paste0("Number of grid cells for wich environment data were extracted = ", nrow(agg.check)),con = fileConn)
    writeLines(paste0("Number of environmental variables applied = ", ncol(env.dat)),con = fileConn)
    writeLines(paste0("Environmental variables extracted:" ),con = fileConn)
    for(i.env in 1:ncol(env.dat))
      {
      writeLines(paste0("Var ",i.env,": ", colnames(env.dat)[i.env]), con = fileConn)
      }# end for i.env
    writeLines("#######################################################################",con = fileConn)
    close(fileConn) #**************************************************************************
  } # end if !is.null(output.folder)
  
  # write some feedback to the terminal
  if(verbose)
  {
    msg1 = 'Returned object is a dataframe.'
    msg2 = paste(ncol(env.dat),' environment variables have been extracted for ', nrow(agg.check), ' grid cells.')
    if(!is.null(output.folder)){
      msg3 = paste('These data have been also been written to ', out.path)
      cat(paste(msg1, msg2, msg3, sep = '\n'))      
    }else{
      cat(paste(msg1, msg2, sep = '\n'))
    }
  }# end if verbose
  
  # And hand back the aggregated data
  return(agg.check)

} # end extract_env_data()
