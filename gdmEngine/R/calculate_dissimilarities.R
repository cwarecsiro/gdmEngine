#'@title Calculate dissimilarities
#'
#'@description Calculates the compositional dissimilarity between all specified site-pairs. Fills the dissimilarities in the input dataframe 'pairs.table'. Note this uses an approach that trades processing time for memory efficiency, so it can handle very large numbers of site pairs on any machine.
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


# composition.data = Selected.records
#head(composition.data)

calculate_dissimilarities <- function(pairs.table, 
                                      composition.data,
                                      output.folder = NULL,       
                                      output.name = "pairs_table_dissim",  
                                      verbose=TRUE) 
{
  # First indexify the sites and species in composition.data
  composition.data$Site.ID <- as.factor(paste(composition.data$decimalLongitude, composition.data$decimalLatitude, sep = '_'))
  composition.data.site.indices<-levels(composition.data$Site.ID)
  
  # CW: account for scenario where scientific name is not factor
  if (class(composition.data$scientificName) == 'character'){
    
    # do it this way to avoid creating a column with class factor where it 
    # was previously char... just because factors are trouble.
    composition.data.spp.indices <- unique(composition.data$scientificName)
    
  } else {
    
    composition.data$scientificName = as.factor(composition.data$scientificName)
    composition.data.spp.indices <- levels(composition.data$scientificName)
    
  }
  
  #composition.data.spp.indices<-levels(composition.data$scientificName)
  composition.data$site.index <- match(composition.data$Site.ID, composition.data.site.indices)  
  composition.data$spp.index <- match(composition.data$scientificName, composition.data.spp.indices)
  
  # generalise
  site.spp.index <- as.matrix(composition.data[,c('site.index', 'spp.index')]) 
  
  # Then get the index for each site in the pairs.table
  pairs.table$s1.site.ID <- paste(pairs.table$s1.decimalLongitude, pairs.table$s1.decimalLatitude, sep = '_')
  pairs.table$s2.site.ID <- paste(pairs.table$s2.decimalLongitude, pairs.table$s2.decimalLatitude, sep = '_')
  pairs.table$S1.index <- match(pairs.table$s1.site.ID, composition.data.site.indices)
  pairs.table$S2.index <- match(pairs.table$s2.site.ID, composition.data.site.indices)
  #pairs.site.index <- as.matrix(pairs.table[,c(13,14)])
  # generalise
  pairs.site.index <- as.matrix(pairs.table[,c('S1.index', 'S2.index')])
  
  # Determine the richness of each site 
  site.richness <- tabulate(composition.data$site.index)
  max.richness <- as.integer(max(site.richness))
  
  # Now run some rcpp code to format the data & calculate dissimilarities for the selected pairs
  distance <- PairsDissim(site.spp.index, 
                          pairs.site.index,
                          site.richness,
                          max.richness)
  
  # create a new dataframe to return, with the scaled dissimilarities
  #pairs.table.new <- pairs.table[,-c(11:14)] # If we want to remove the xy site names
  # generalise:
  pairs.table.new <- pairs.table[,-c('s1.site.ID', 's2.site.ID')] # If we want to remove the xy site names
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
  
} # end calculate_dissimilarities


##-------------------------------------------------------------------------------------------------------------##
# #'@title Compositional dissimilarity calculation for specified pairs of sites
# #'
# #'@description Calculate the dissimilarity between specified pairs of sites. 
# #'
# #'@param site_spp (matrix) A matrix of the indices for species (col2) occurring in each site (col1), with each row an occurrence
# #'@param pair_rows (matrix) A matrix of the site index for site 1 (col 1) and site 2 (col 2) in each pair of sites (rows)
# #'@param site_rich (vector) A vector giving the total number of species in each site (element)
# #'@param max_richness (integer)
# #'
# #'@returns Vector, the Sorensen dissimilarity between the specified pairs.
# #'
# #'@examples output = PairsDissim(site.spp, ij.pairs, site.rich, max.rich)
# #'@importFrom Rcpp cppFunction
# #'
# #'@export
# cppFunction('NumericVector PairsDissim(IntegerMatrix site_spp, IntegerMatrix pair_rows, IntegerVector site_rich, int max_richness) {
#             
#   int n_sites = site_rich.size();
#   int n_pairs = pair_rows.nrow();
#   int n_records = site_spp.nrow();  
#   IntegerMatrix comp(n_sites,max_richness);        
#   IntegerVector upto_index(n_sites);
#   NumericVector out(n_pairs);
# 
# 
#   for(int i_site=0; i_site < n_sites; i_site++) {
#     upto_index[i_site] = 0;
#     }
#   for(int i_rec = 0; i_rec < n_records; i_rec++) {
#     int site_index = site_spp(i_rec,0) - 1;
#     comp(site_index, upto_index[site_index]) = site_spp(i_rec,1);
#     upto_index[site_index] += 1;
#     }
# 
#   for(int i = 0; i < n_pairs; i++) {
#     int site_one_index = pair_rows(i,0) - 1;
#     int site_two_index = pair_rows(i,1) - 1;
#     float n_spp_common = 0;
#     for(int i_spp_one=0; i_spp_one<site_rich[site_one_index]; i_spp_one++){
#       for(int i_spp_two=0; i_spp_two<site_rich[site_two_index]; i_spp_two++){
#         if(comp(site_one_index,i_spp_one) == comp(site_two_index,i_spp_two)){
#           n_spp_common += 1;
#           }
#         }      
#       }  
#     float sum_rich = site_rich[site_one_index] + site_rich[site_two_index];
#     out[i] = (1 - ((2 * n_spp_common) / (sum_rich)));
#     }
#   return out;
#   }')
# ##-------------------------------------------------------------------------------------------------------------##
