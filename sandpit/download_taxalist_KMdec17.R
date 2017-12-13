#'@title Download occurrences for a list of taxa from the ALA 
#'
#'@description Takes a vector of taxonomic names and searches the ALA for these. Can save the data and/or load to memory. 
#'
#'@param specieslist (vector) Vector of taxonomic names to search for.
#'@param output.folder (string) Optional dstination filepath to write data to. If missing, tempdir is used.
#'@param parallel (boolean) Split list across multiple CPUs (all available -2). Default FALSE.
#'@param ... Additional named arguments to pass to gdmEngine::download_ala which is called internally.
#'
#'@return std.output
#'
#'@examples download_specieslist(c('spp1', 'spp2'), 'tempdir()'C:/Users/xyx, read = FALSE)
#'
#'@export
download_taxalist = function(specieslist, 
                             output.folder = NULL, 
                             parallel = FALSE, 
                             ...)
  {
  ## add in start-fresh or overwrite flag
  ## make an option for just downloading zip files (probably with a log file).
  
  ## check args
  if(!is.null(output.folder)){
    if(!dir.exists(output.folder)) dir.create(output.folder)
    } else {
    output.folder = tempdir()
    }# end else !is.null...
  ## in case the list is derived from a data.frame
  if(is(specieslist, 'factor')){
    specieslist = as.character(paste(specieslist))
    }# end if is(specieslist...
  
  ## check filepath ending
 output.folder = check_filepath(output.folder)
  
  ## config output.folder for raw downloads
  raw_files = paste0(output.folder , 'raw_files')
  if(dir.exists(raw_files)){
    cat('Warning: output.folder exists and contents will be written over...')
    } else {
    dir.create(raw_files)
    }#end else dir.exists... 
  
  ## default args
  download_args = list(
    method = "indexed", 
    download_reason_id = 7,
    output.folder = raw_files,
    read = FALSE,
    verbose = FALSE)
  
  ## check for user download args
  user_args = list(...)
  if (length(user_args)) {
    ## all opts of gdmEngine::download_occurrences
    opts = formals(download_occurrences)
    for (i in 1:length(user_args)) {
      this_opt <- names(user_args)[i]
      if (! (this_opt %in% names(opts))) {
        cat(paste0("'", this_opt, "'", ' is not a valid option. Should be one of:'), 
            names(opts), sep = '\n')
        } else {
        download_args[[this_opt]] = paste(user_args[i])
        }
      } # end user_opts 
   } # end if length(user_args)
  
  ## log container
  log = list()
  log$ALA_DOWNLOAD = Sys.Date()
  log_f = file.path(output.folder, paste0('ALA_download_R_log_', Sys.Date(), '.RData'))
  save(log, file = log_f)

  ## run
  try(repeat{
    
    load(log_f)
    not_done <- outersect(specieslist, names(log)[-1])
    
    ## Break check
    if(length(not_done)==0) break 
    
    ## Not finished...? Start looping again.
    
    if(parallel){
      no_cores = parallel::detectCores() - 2
      cl = parallel::makeCluster(no_cores)
      clusterExport(cl, list(), envir=environment())
      clusterExport(cl, list('download_args', 'download_occurrences',
                             'check_filepath', 'counter', 'outersect'))
      #clusterEvalQ(cl, library(ALA4R))
      clusterEvalQ(cl, library(gdmEngine))
     
      cat(paste0('Searching the ALA for records of ', length(not_done),
                 ' taxa...'))
      
      check = parLapply(cl, not_done, function(x){
        tryCatch({
          do.call(download_occurrences, 
          c(list(taxon = x), download_args))
        }, 
        error = function(e) {
          e
        })
      })
      
      ## log
      check = lapply(check, function(x)
        if(is.null(x)) x = 'Successfully downloaded file')
      log = c(log, check)
      names(log)[2:length(log)] = not_done
      
      gc()
      
      stopCluster(cl = cl)
    
      } else { # end if parallel
    
      ## loop
      for(spp in not_done){
        
        ## counter to print to console
        counter(msg = 'Searching the ALA for records of', iterable = spp)
        
        ## Download records from ALA 
        ## Log any errors in downloads by catching them and passing them
        ## through to the log. 
        check <- tryCatch(ala <- do.call(download_occurrences, 
                                         c(list(taxon = spp), download_args)),
                          error = function(e) e)
        
        
        ## log
        if(!is.null(check)){ ## signifies error
          log[[spp]] = check
        } else {
          log[[spp]] = 'Successfully downloaded file'
        } 
        
        ## Give the ALA a break  
        ## Sys.sleep(5)
        
      } # end spp loop
    
    gc()
    
    } # end else serial
    
    ## save log
    save(log, file = log_f)
  
  }) #end try(repeat...
  
  #**************************************************************************
  # write a text log file describing how the data was created *************************************
  fileConn<-file(file.path(output.folder,paste0("TaxaList_Download_",Sys.Date(),"_log_file.txt")),'w')
  writeLines("#######################################################################",con = fileConn)
  writeLines("###",con = fileConn)
  writeLines("### ALA data download log file ",con = fileConn)
  writeLines("###",con = fileConn)
  writeLines(paste0("### Created ",Sys.time()," using the download_taxalist() function."),con = fileConn)
  writeLines("###",con = fileConn)
  writeLines("#######################################################################",con = fileConn)
  writeLines("",con = fileConn)
  writeLines(paste0("Output data folder = ", raw_files),con = fileConn)
  writeLines(paste0("Number of species queried = ", length(specieslist)),con = fileConn)
  writeLines(paste0("TaxonName DownloadSuccessful"),con = fileConn)
  for(i.spp in 2:length(log))
    {
    dld.result <- "no"
    if(log[[i.spp]]=='Successfully downloaded file')
      {
      dld.result <- "yes"
      } # end if
    writeLines(paste0(names(log[2]), " ", dld.result),con = fileConn)
    }# end for i.spp
  writeLines("#######################################################################",con = fileConn)
  close(fileConn) #**********************************************************
  
  
  ## summary
  search_summary = paste('Found records for ', 
                         length(log[log == 'Successfully downloaded file']),
                         'of', length(specieslist), 'listed taxa')
  cat('\n')
  cat('Completed occurrence searches', sep = '\n')
  cat(search_summary)
  
}

