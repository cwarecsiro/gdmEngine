
#'@title Simple counter to monitor a loop
#'
#'@description Counter that doesn't clog the console but enables progress to be monitored
#'@param msg (string) Default: 'Doing loop.' Whatever needs to be printed...
#'@param vec (vector) A vector of things to iterate through with each i in vec being printed as progress is made. Default: NULL. 
#'@param iterable An object that is iterated over in a loop (i.e. the i in the above example)
#'@return Output to console
#'@examples for (i in LETTERS[1:26]) counter(iterable = i)
#'@export
counter = function(msg = 'Doing loop ', vec = NULL, iterable = NULL){
  
  if(!exists('counter_idx')) counter_idx <<- 1
  if(!is.null(vec)){
    msg = paste(msg, vec[counter_idx])
  } else if(!is.null(iterable)){
    msg = paste0(msg, ' ', iterable)
  } else {
    msg = paste0(msg, counter_idx)
  }
  
  cat("\r", msg)
  counter_idx <<- counter_idx + 1
  Sys.sleep(0.2)
  
}


#'@title Get the non-matching objects from two vectors
#'
#'@description Returns objects not common among two vectors. The opposite to intersect.
#'@param x (vector)
#'@param y (vector)
#'@return vector
#'@examples outersect(LETTERS[1:10], LETTERS[5:20])
#'@export
#'@notes https://www.r-bloggers.com/outersect-the-opposite-of-rs-intersect-function/
outersect = function(x, y) {
  sort(c(setdiff(x, y), setdiff(y, x)))
}


#'@title Unzip file downloaded from the ALA and read CSV into memory
#
#'@description Unzips file and handles exceptions. 
#'@param fn (string) Filepath to .zip archive.
#'@param memory_only (boolean) If TRUE (default) will delete unzipped data once read into memory
#'@returns data.frame
#'@examples output = unzipper('C:/Users/species_x.zip')
#'@export
unzipper = function(fn, memory_only = TRUE){
  temp_dump = paste(dirname(fn), gsub('.zip', '', basename(fn)), sep = '/')
  dir.create(temp_dump)
  df = tryCatch(unzip(fn, exdir = temp_dump)[1] %>% 
                  read.csv(stringsAsFactors=FALSE),
                error = function(e) e)
  if(memory_only) unlink(temp_dump, recursive = TRUE) 
  return(df)
}


#'@title Ensures filepath has trailing forwardslash
#
#'@description Adds trailing forwardslash to filepath if not present so that the filepath can be added to for writing.
#'@param filepath (string) 
#'@returns string (filepath with trailing forwardslash)
#'@examples fn = checkfilepath('C:/Users')
#'@export
check_filepath = function(filepath){
  check_nchar = substr(filepath, nchar(filepath), nchar(filepath))
  if(check_nchar == '/'){
    return(filepath)
  } else if(check_nchar == '\\\\') {
    return(filepath)
  } else {
    return(paste0(filepath, '/'))
  }
}
