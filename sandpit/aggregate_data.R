#'@title Aggregate ALA data to grid cells
#'
#'@description Take filtered occurrence data and group their position to grid cell centres. 
#'
#'@param src (string) Filepath to directory where loadeds are saved. Must not be other .zip files stored in src.
#'@param output.folder (string) A folder to save the outputs to. If none specified, no file is written.
#'@param output.name (string), A name to use in saving the outputs. Default: 'merged_taxa_data'.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@returns List. $data is a data.frame with the combined records. $log is a data.frame with the taxa names searched, 
#'records returned, and any errors.
#'
#'@examples output = merge_downloads('C:/Users/raw_files')
#'
#'@export
merge_downloads = function(ALA.filtered.data, 
                           output.folder = NULL,       
                           output.name = "aggregated_taxa_data",
                           verbose = TRUE)
{
  
} # end merge_downloads() 