#'@title Aggregate ALA data to grid cells
#'
#'@description Take filtered occurrence data and group their position to grid cell centres. 
#'
#'@param ALA.filtered.data (string) Filepath to directory where loadeds are saved. Must not be other .zip files stored in src.
#'@param domain.mask, (raster layer) A raster layer specifying the analysis domain
#'@param agg.radius.ncells (float) A spatial radius, specified as number of grid cells, over which to aggregate records. If NULL, no radial aggregation is performed. Default: NULL.
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string), A name to use in saving the outputs. Default: 'aggregated_taxa_data'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@returns Dataframe. Also writes this dataframe to file if output.folder is specified. 
#'
#'@examples output = aggregate_data(my.data, my.raster)
#'
#'@export
aggregate_ALA_data = function(ALA.filtered.data,
                              domain.mask,
                              agg.radius.ncells = NULL,
                              output.folder = NULL,       
                              output.name = "aggregated_taxa_data",
                              verbose = TRUE)
{

  
  # Aggregate the records to the centres of the grid cells they occur in
  xy <- cbind(ALA.filtered.data$decimalLongitude, ALA.filtered.data$decimalLatitude)
  origin <- c(domain.mask@extent@xmin, domain.mask@extent@ymin)
  grid_res <- res(domain.mask)
  out <- as.data.frame(t(apply(xy, 1, function(z) grid_res/2 + origin + grid_res*(floor((z - origin)/grid_res)))))
  #replace the long, lat with the grid cell centred version
  ALA.filtered.data$decimalLongitude<-out[,1]
  ALA.filtered.data$decimalLatitude<-out[,2]

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ## IF WE'RE AGGREGATING RECORDS ACROSS NEIGHBOURING GRID CELLS, DO SOME MORE COMPLEX AGGREGATION
  if(!is.null(agg.radius.ncells))
    {
    # determine the radius in spatial units
    radius <- agg.radius.ncells * xres(domain.mask)
    
    # create a grid-cell name 'xy'
    ALA.filtered.data$xy <- paste(ALA.filtered.data$decimalLongitude, ALA.filtered.data$decimalLatitude, sep = '_')
  
    # Create a new data frame with the number of records in each grid cell
    nrec <- count(ALA.filtered.data$xy)
    colnames(nrec)<- c("xy","n.records")
    nrec$xy <- as.character(nrec$xy)
    ll <- strsplit(nrec[,1], '_')
    nrec$decimalLongitude <- unlist(ll)[ c(TRUE,FALSE) ]
    nrec$decimalLatitude <- unlist(ll)[ c(FALSE,TRUE) ]
    nrec$decimalLongitude <- as.numeric(nrec$decimalLongitude)
    nrec$decimalLatitude <- as.numeric(nrec$decimalLatitude)
  
    # Now for each grid cell with records in it, find how many records there'd be if we included 
    # records form gid cells within the specified radius [note: this is quite slow, maybe a more efficient way to code this, or use rcpp ?]
    nrec$neighbour.recs <- 0
    for(i.row in 1:nrow(nrec))
      {
      x.cell <- nrec$decimalLongitude[i.row]
      y.cell <- nrec$decimalLatitude[i.row]  
      neighbours <- nrec
      neighbours$cell.dist <- sqrt((nrec$decimalLongitude - x.cell)^2 + (nrec$decimalLatitude - y.cell)^2)
      neighbours <- neighbours[which(neighbours$cell.dist <= radius),]
      nrec$neighbour.recs[i.row] <- sum(neighbours$n.records)
      }# end for i.row
    
    # Order the cells by most no. records in neighbourhood (descending)
    nrec.reduced <- nrec[order(nrec$neighbour.recs, nrec$n.records, decreasing=TRUE),]
    # put the cell with the most records in it's neighbourhood somewhere special
    best.record.summary <- nrec.reduced[1,]
    # remove all cells within the neighbourhood of the best cell (including the best cell)
    x.cell <- nrec.reduced$decimalLongitude[1]
    y.cell <- nrec.reduced$decimalLatitude[1]
    cell.dist <- sqrt((nrec.reduced$decimalLongitude - x.cell)^2 + (nrec.reduced$decimalLatitude - y.cell)^2)
    neigh.cells <- which(cell.dist <= radius)
    nrec.reduced <- nrec.reduced[-neigh.cells,]

    # Now start a while loop through the rest of the data, keeping only the cells with the most 
    # records in their neighbourhood
    while(nrow(nrec.reduced)>0)
      {
      best.record.summary <- rbind(best.record.summary, nrec.reduced[1,])
      x.cell <- nrec.reduced$decimalLongitude[1]
      y.cell <- nrec.reduced$decimalLatitude[1]
      cell.dist <- sqrt((nrec.reduced$decimalLongitude - x.cell)^2 + (nrec.reduced$decimalLatitude - y.cell)^2)
      neigh.cells <- which(cell.dist <= radius)
      nrec.reduced <- nrec.reduced[-neigh.cells,]
      }# end while
    
    # Now re-allocate longitude and latitude to the aggregated cell centres. For each cell in nrec, 
    # find the closest cell in the aggregated neighbourhood cells, and give the record that x,y 
    # if it is within the right radius (otherwise, leave it alone).
    for(i.row in 1:nrow(nrec))
      {
      x.cell <- nrec$decimalLongitude[i.row]
      y.cell <- nrec$decimalLatitude[i.row]  
      rec.dist <- sqrt((best.record.summary$decimalLongitude - x.cell)^2 + (best.record.summary$decimalLatitude - y.cell)^2)      
      if(min(rec.dist) < radius)
        {
        nrec$decimalLongitude[i.row] <- best.record.summary$decimalLongitude[which.min(rec.dist)]
        nrec$decimalLatitude[i.row] <- best.record.summary$decimalLatitude[which.min(rec.dist)]  
        }# end if min(rec.dist) < radius
      } # end for i.row
    
    # And finally, change the longitude & latitude in the records, to match the aggregated values
    ALA.filtered.data.agg <- ALA.filtered.data[,-which(names(ALA.filtered.data) %in% c("decimalLatitude", "decimalLongitude"))]
    nrec <- nrec[,which(names(nrec) %in% c("xy", "decimalLongitude", "decimalLatitude"))]
    ALA.filtered.data.agg <- merge(ALA.filtered.data.agg, nrec, by="xy", all.x=TRUE,  sort=FALSE)
    ALA.filtered.data <- ALA.filtered.data.agg[,-which(names(ALA.filtered.data.agg) %in% c("xy"))]
    } # end if !is.null(agg.radius.ncells)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
  # Get the number of grid cells that the data has been aggregated to
  cells.location <- ALA.filtered.data[,which(names(ALA.filtered.data) %in% c("decimalLongitude", "decimalLatitude"))]
  cells.location <- unique(cells.location)
  
  # Now write out the data to file (if specified) and return the aggregated records
  # write the data to file, if an output folder is specified
  if(!is.null(output.folder))
  {
    if(!dir.exists(output.folder))
    {
      dir.create(output.folder)
    }# end if !dir.exists
    out.path <- file.path(output.folder,paste0(output.name,"_",Sys.Date(),".csv")) 
    write.csv(ALA.filtered.data, out.path, row.names=FALSE)
    # write a log file describing how the data was created *************************************
    fileConn<-file(file.path(output.folder,paste0(output.name,"_",Sys.Date(),"_log_file.txt")),'w')
    writeLines("#######################################################################",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("### ALA data aggregation log file ",con = fileConn)
    writeLines("###",con = fileConn)
    writeLines(paste0("### Created ",Sys.time()," using the aggregate_ALA_data() function."),con = fileConn)
    writeLines("###",con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    writeLines("",con = fileConn)
    writeLines(paste0("Output data file = ", out.path),con = fileConn)
    writeLines(paste0("Domain mask applied = ", domain.mask@file@name),con = fileConn)
    if(is.null(agg.radius.ncells)){
      writeLines(paste0("No radial cell aggregation applied "),con = fileConn)
      }else{
      writeLines(paste0("Data aggregated to grid cells within radius of ", agg.radius.ncells," cells."),con = fileConn)  
      } # end if else is.null(agg.radius.ncells)
    writeLines(paste0("Number of records = ", nrow(ALA.filtered.data)),con = fileConn)
    writeLines(paste0("Number of grid cells records are aggregated to = ", nrow(cells.location)),con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    close(fileConn) #**************************************************************************
  } # end if !is.null(output.folder)
  
  # write some feedback to the terminal
  if(verbose)
    {
    msg1 = 'Returned object is a dataframe.'
    msg2 = paste('These data have been also been written to ', out.path)
    msg3 = paste(nrow(ALA.filtered.data), ' records have been aggregated to ', nrow(cells.location),' grid cells.')
    cat(paste(msg1, msg2, msg3, sep = '\n'))
    }# end if verbose
  
  # And hand back the aggregated data
  return(ALA.filtered.data)
  
} # end aggregate_ALA_data() 
