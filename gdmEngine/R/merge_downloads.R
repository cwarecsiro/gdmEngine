#'@title Merge datasets downloaded from the ALA
#'
#'@description Unzip and load CSV files downloaded from the ALA. 
#'
#'@param src (string) Filepath to directory where loadeds are saved. Must not be other .zip files stored in src.
#'@param keep_unzip (boolean) Save the unzipped contents of each file?
#'@param parallel (boolean) If TRUE (default) will use all CPU (-2) available to load data.
#'@param verbose (boolean) Print messages to console. Default TRUE.
#'
#'@returns List. $data is a data.frame with the combined records. $log is a data.frame with the taxa names searched, 
#'records returned, and any errors.
#'
#'@examples output = merge_downloads('C:/Users/raw_files')
#'
#'@export
merge_downloads = function(src, keep_unzip = FALSE, parallel = TRUE, 
                           verbose = TRUE){
  
  filelist = list.files(src, full.names = TRUE, pattern = '.zip')  
  
  unzipper = function(fn, memory_only = TRUE){
    temp_dump = paste(dirname(fn), gsub('.zip', '', basename(fn)), sep = '/')
    dir.create(temp_dump)
    df = tryCatch(unzip(fn, exdir = temp_dump)[1] %>% 
                    read.csv(stringsAsFactors=FALSE),
                  error = function(e) e)
    if(memory_only) unlink(temp_dump, recursive = TRUE) 
    return(df)
  }
  
  
  if(parallel){
    no_cores = parallel::detectCores() - 2
    cl = parallel::makeCluster(no_cores)
    clusterExport(cl, list(), envir=environment())
    clusterExport(cl, list('download_args', 'download_occurrences',
                           'check_filepath', 'counter', 'outersect', 'unzipper', 'filelist'))
    clusterEvalQ(cl, library(dplyr))  
    clusterEvalQ(cl, library(gdmEngine))  
    cat(paste0('Extracting data from ', length(filelist),
               ' files...'))
    
    ptm = proc.time()
    datasets = parLapply(cl, filelist, unzipper)
    proc.time()-ptm
    
    gc()
    
    stopCluster(cl = cl)
    
  } else { # end if parallel
    
    datasets = lapply(filelist, unzipper)
    
  } # end if serial
  
  ## check for anything that is not a data.frame
  check_df = sapply(datasets, function(x) is(x, 'data.frame'))
  not_df = which(!check_df)
  
  ## create log
  log = data.frame(matrix(nrow = length(filelist), ncol = 0))
  log$taxa = sapply(filelist, function(x) gsub('.zip', '', basename(x)))
  log$n_records = sapply(datasets, nrow)
  
  ## deal with not_df
  for (x in not_df){
    if(is(datasets[[x]], 'simpleError')){
      log$n_records[[x]] = gsub(' ', '_', datasets[[x]]$message)
    } else {
      log$n_records[[x]] = 'ERROR'
    }
  }
  
  ## rm NULL datasets and rbind
  datasets = datasets[-c(not_df)]
  data_bound = rbindlist(datasets)
  
  ## output
  bundle = list(data = as.data.frame(data_bound),
                log = log)
  if(verbose){
    msg1 = 'Returned object is a list.'
    msg2 = paste('$data is a data.frame with ', nrow(data_bound), 'records')
    msg3 = '$log is a summary data.frame.'
    cat(paste(msg1, msg2, msg3, sep = '\n'))
  }
  
  return(bundle)
  
}
