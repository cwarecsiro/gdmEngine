#'@title Calculate dissimilarities
#'
#'@description Calculates the compositional dissimilarity between all specified site-pairs. Fills the dissimilarities in the input dataframe 'pairs.table'. Note this uses an approach that trades processing time for memory efficiency, so it can handle very large numbers of site pairs on any machine.
#'
#'@param pairs.table (dataframe) A dataframe holding the site-pairs. This is the first six columns of a gdm input table.
#'@param composition.data (dataframe) A dataframe holding the final species composition data: the records of all species in all grid cells to be used for modelling.
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string) A name to use in saving the outputs. Default: 'pairs_table_dissim'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@returns Dataframe.
#'
#'@examples output = calculate_dissimilarities(My.pairs.table, My.composition.data, output.folder = 'C:/Users/processed_data', output.name = 'My.sitepair.data')
#'
#'@importFrom reshape2 dcast
#'@importFrom betapart beta.pair
#'
#'@export
calculate_dissimilarities <- function(pairs.table, 
                                      composition.data,
                                      output.folder = NULL,       
                                      output.name = "pairs_table_dissim",  
                                      verbose=TRUE) 
  {
    # Calculate obs number of species shared betwen all the pairs
    # Process a specified number of pairs at a time
    n.pairs.proc <- 40
    n.procs <- ceiling(nrow(pairs.table) / n.pairs.proc)
    # Add 'xy' site IDs to the composition data and the pairs table
    composition.data$Site.ID <- paste(composition.data$decimalLongitude, composition.data$decimalLatitude, sep = '_')
    pairs.table$s1.site.ID <- paste(pairs.table$s1.xCoord, pairs.table$s1.yCoord, sep = '_')
    pairs.table$s2.site.ID <- paste(pairs.table$s2.xCoord, pairs.table$s2.yCoord, sep = '_')
    # Create a list to catch the dissimilarities
    distance <- NULL    
    for(i.proc in 1:n.procs)
      {
      # set the start and end of the processing zone
      i.start.row <- ((i.proc - 1) * n.pairs.proc) + 1
      i.end.row <- (i.proc * n.pairs.proc)
      if(i.end.row > nrow(pairs.table))
        {i.end.row <-  nrow(pairs.table)}
      # Trim the composition table to only include sites in the selected pairs
      sites.proc <- c(as.character(pairs.table$s1.site.ID[c(i.start.row:i.end.row)]), as.character(pairs.table$s2.site.ID[c(i.start.row:i.end.row)]))
      sites.proc <- unique(sites.proc)
      composition.table.proc <- composition.data[as.character(composition.data$Site.ID) %in% sites.proc, ]
      composition.table.proc <- unique(composition.table.proc)
      # Convert to sites x spp table
      present <- rep(1, times=nrow(composition.table.proc))
      composition.table.proc <- cbind(composition.table.proc, present)
      composition.mat <- dcast(composition.table.proc, Site.ID ~ scientificName, sum, value.var='present')
      row.names(composition.mat) <- composition.mat[,1]
      composition.mat <- composition.mat[,-1]
      # Calculate compositional dissimilaities using the betapart package
      composition.dis <- beta.pair(composition.mat, index.family="sorensen")
      dissimilarity.data <- composition.dis$beta.sor
      # Convert dissimilarities to a matrix
      dissimilarity.data <- as.matrix(dissimilarity.data)
      # Get the dissimilarities for our site-pairs of interest from the dissimilarity matrix
      pairs.table.names <- as.matrix(pairs.table[c(i.start.row:i.end.row),c(7,8)])
      mat.row <- match(pairs.table.names[,1],rownames(dissimilarity.data))
      mat.col <- match(pairs.table.names[,2],colnames(dissimilarity.data))
      dissim.lst <- dissimilarity.data[cbind(mat.row,mat.col)]
      # Now join these dissimilarities to our running list
      distance <- c(distance, dissim.lst)
      } # end for i.proc

  # create a new dataframe to return, with the scaled dissimilarities
  pairs.table.new <- pairs.table[,-c(7,8)]
  pairs.table.new$distance <- distance

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
    writeLines(paste0("Number of site-pairs = ", nrow(pairs.table.new)),con = fileConn)
    writeLines(paste0("Dissimilarity metric = Sorensen's"),con = fileConn)
    writeLines(paste0("Average site-pair dissimilarity = ", mean(pairs.table.new$distance)),con = fileConn)
    writeLines(paste0("Proportion of site-pairs with non-complete dissimilarity = ", (sum(pairs.table.new$distance < 1)/nrow(pairs.table.new))),con = fileConn)
    writeLines("#######################################################################",con = fileConn)
    close(fileConn) #**************************************************************************
  } # end if !is.null(output.folder)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#
  
  # write some feedback to the terminal
  if(verbose)
    {
    msg1 = 'Returned object is a dataframe.'
    msg2 = paste('Dissimilarities have been calculated for ', nrow(pairs.table.new), ' site-pairs.')
    msg3 = paste(((sum(pairs.table.new$distance < 1)/nrow(pairs.table.new))*100) , '% of site-pairs have non-complete dissimilarity.')
    if(!is.null(output.folder)){
      msg4 = paste('These data have been also been written to ', out.path)
      cat(paste(msg1, msg2, msg3, msg4, sep = '\n'))      
      }else{
      cat(paste(msg1, msg2, msg3, sep = '\n'))
      }
    }# end if verbose
  
  # Return the freshly filled site-pair table  
  return(pairs.table.new)
    
  } # end Scale.Dissimilarity.NSW.BBA() function
  
