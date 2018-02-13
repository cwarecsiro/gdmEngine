#'@title Select grid cells composition
#'
#'@description Select grid cells with sufficient records for modelling composition. 
#'
#'@param ALA.aggregated.data (dataframe) A dataframe of appropriate format held in R working environment.
#'@param domain.mask (raster layer) A raster layer specifying the analysis domain
#'@param min.richness.threshold (integer) A minimum threshold number of species for selected grid cells. Default: NULL.
#'@param max.richness.threshold (integer) A maximum threshold number of species for selected grid cells.  Default: NULL.
#'@param reference.radius.ncells (float) A radius (in grid cells) to use in assessing relative richness to maximum.  Default: NULL.
#'@param min.proportion.max.richness (float) The minimum threshold proportion of max richness in reference area. Default: 0.25.
#'@param agg.radius.ncells (float) A spatial radius, specified as number of grid cells, over which to aggregate records. If NULL, no radial aggregation is performed. Default: NULL.
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string), A name to use in saving the outputs. Default: 'selected_gridcell_composition'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@return Dataframe. Also writes this dataframe to file if output.folder is specified. 
#'
#'@examples output = aggregate_data(my.data, my.raster)
#'
#'@export
select_gridcells_composition = function(ALA.aggregated.data,
                                        domain.mask,
                                        min.richness.threshold = NULL,
                                        max.richness.threshold = NULL,
                                        reference.radius.ncells = NULL,
                                        min.proportion.max.richness = 0.25,
                                        output.folder = NULL,       
                                        output.name = "selected_gridcell_composition",
                                        verbose = TRUE)
{
  
  # First we need to determine the number of species observed in each grid cell
  # create a grid-cell name 'xy'
  ALA.aggregated.data$xy <- paste(ALA.aggregated.data$decimalLongitude, ALA.aggregated.data$decimalLatitude, sep = '_')
  agg.check <- count(ALA.aggregated.data, c('xy', 'scientificName'))
  agg.check$freq <- ifelse(agg.check$freq > 1, 1, agg.check$freq)
  agg.check <- as.data.frame(table(agg.check[,c(1,3)]))
  agg.check <- agg.check[,-2]
  agg.check$xy <- as.character(agg.check$xy)
  colnames(agg.check) [2] <- "richness"
  agg.check$select <- 1
  nrecs.initial <- nrow(ALA.aggregated.data) 
  
  # If we're applying a minimum richness threshold, remove those cells
  if(!is.null(min.richness.threshold))
    {
    agg.check$select[agg.check$richness < min.richness.threshold] <- 0
    }#end if !is.null(min.richness.threshold) 

  # If we're applying a minimum richness threshold, remove those cells
  if(!is.null(max.richness.threshold))
    {
    agg.check$select[agg.check$richness > max.richness.threshold] <- 0
    }#end if !is.null(max.richness.threshold) 
  
  # If we're applying the more complex min proportion of maximum richness in the 
  # surrounding area threshold, then do it laddie.
  if(!is.null(reference.radius.ncells))
    {
    # determine the radius in spatial units
    radius <- reference.radius.ncells * xres(domain.mask) 
    # add the x/y values back into the working dataframe
    ll <- strsplit(agg.check[,1], '_')
    agg.check$decimalLongitude <- unlist(ll)[ c(TRUE,FALSE) ]
    agg.check$decimalLatitude <- unlist(ll)[ c(FALSE,TRUE) ]
    agg.check$decimalLongitude <- as.numeric(agg.check$decimalLongitude)
    agg.check$decimalLatitude <- as.numeric(agg.check$decimalLatitude)
    # Now loop through each cell, find the max richness in the radius, run the
    # threshold test, and either leave it or remove it.
    # Reduce the dataframe to only those that have yet to be omitted.
    agg.check.rad <- agg.check[which(agg.check$select>0),]
    # Then loop
    for(i.row in 1:nrow(agg.check.rad))
      {
      x.cell <- agg.check.rad$decimalLongitude[i.row]
      y.cell <- agg.check.rad$decimalLatitude[i.row]  
      neighbours <- agg.check.rad
      neighbours$cell.dist <- sqrt((neighbours$decimalLongitude - x.cell)^2 + (neighbours$decimalLatitude - y.cell)^2)
      neighbours <- neighbours[which(neighbours$cell.dist <= radius),]
      rich.thrsh <- max(neighbours$richness) * min.proportion.max.richness 
      if(agg.check.rad$richness[i.row] < rich.thrsh)
        { 
        agg.check$select[which(agg.check$xy == agg.check.rad$xy[i.row])] <- 0
        }#end if(agg.check.rad$richness...)
      }# end for i.row
    }#end if !is.null(reference.radius.ncells) 
  
  # Now we've finished selecting cells, format up the records for output
  agg.check <- agg.check[,-which(names(agg.check) %in% c("decimalLatitude", "decimalLongitude"))]
  ALA.aggregated.data <- merge(ALA.aggregated.data, agg.check, by="xy", all.x=TRUE,  sort=FALSE)
  # remove the relevant records
  ALA.aggregated.data <- ALA.aggregated.data[which(ALA.aggregated.data$select > 0), ]
  # refine to the key columns
  ALA.aggregated.data <- ALA.aggregated.data[,which(names(ALA.aggregated.data) %in% c("decimalLongitude", "decimalLatitude", "scientificName"))]
  # remove duplicate records within a grid cell
  ALA.aggregated.data <- unique(ALA.aggregated.data)
  
  # Now write out the data to file (if specified) and return the aggregated records
  # write the data to file, if an output folder is specified
  if(!is.null(output.folder))
  {
    if(!dir.exists(output.folder))
    {
      dir.create(output.folder)
    }# end if !dir.exists
    out.path <- file.path(output.folder,paste0(output.name,"_",Sys.Date(),".csv")) 
    write.csv(ALA.aggregated.data, out.path, row.names=FALSE)
    # write a log file describing how the data was created *************************************
    fileConn<-file(file.path(output.folder,paste0(output.name,"_",Sys.Date(),"_log_file.txt")),'w')
    writeLines("#######################################################################",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("### ALA data grid cell selection log file ",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines(paste0("### Created ",Sys.time()," using the select_gridcells_composition() function."),con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    writeLines("",con = fileConn)
    writeLines(paste0("Output data file = ", out.path),con = fileConn)
    writeLines(paste0("Domain mask applied = ", domain.mask@file@name),con = fileConn)
    if(is.null(min.richness.threshold)){
      writeLines(paste0("No minimum richness threshold applied."),con = fileConn)
      }else{
      writeLines(paste0("Minimum richness threshold applied = ", min.richness.threshold," species."),con = fileConn)
      }
    if(is.null(max.richness.threshold)){
      writeLines(paste0("No maximum richness threshold applied."),con = fileConn)
      }else{
      writeLines(paste0("Maximum richness threshold applied = ", max.richness.threshold," species."),con = fileConn)
      }
    if(is.null(reference.radius.ncells)){
      writeLines(paste0("No threshold richness within neighbourhood applied "),con = fileConn)
      }else{
      writeLines(paste0("Minimum richness proportion of ", min.proportion.max.richness," of the maximum richness within radius of ", reference.radius.ncells," cells."),con = fileConn)  
      } # end if else is.null(agg.radius.ncells)
    writeLines(paste0("Number of records before grid cell selection = ", nrecs.initial),con = fileConn)
    writeLines(paste0("Number of records after grid cell selection = ", nrow(ALA.aggregated.data)),con = fileConn)
    writeLines(paste0("Number of grid cells before selection = ", nrow(agg.check)),con = fileConn)
    writeLines(paste0("Number of grid cells after selection = ", sum(agg.check$select)),con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    close(fileConn) #**************************************************************************
  } # end if !is.null(output.folder)
  
  # write some feedback to the terminal
  if(verbose)
    {
    msg1 = 'Returned object is a dataframe.'
    msg2 = paste(sum(agg.check$select),' grid cells have been selected, containing ', nrow(ALA.aggregated.data), ' species records.')
    if(!is.null(output.folder)){
      msg3 = paste('These data have been also been written to ', out.path)
      cat(paste(msg1, msg2, msg3, sep = '\n'))      
      }else{
      cat(paste(msg1, msg2, sep = '\n'))
      }
    }# end if verbose
  
  # And hand back the aggregated data
  return(ALA.aggregated.data)
  
} # end select_gridcells_composition()
  