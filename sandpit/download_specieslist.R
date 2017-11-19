#'@title Download occurrences for a list of taxa from the ALA 
#'
#'@description Takes a vector of taxonomic names and searches the ALA for these. Can save the data and/or load to memory. 
#'
#'@param specieslist (vector) Vector of taxonomic names to search for.
#'@param dst (string) Destination filepath to write data to.
#'@param ... Additional named arguments to pass to download_ala
#'
#'@return std.output
#'
#'@examples download_specieslist(c('spp1', 'spp2'), tempdir(), read = FALSE)
#'
#'@export


download_specieslist = function(specieslist, dst, ...){
  
  ## add in start-fresh or overwrite flag
  ## make an option for just downloading zip files (probably with a log file).
  
  ## checks
  if(!dir.exists(dst)) dir.create(dst)
  assert_that(dir.exists(dst))
  ## in case the list is derived from a data.frame
  if(is(specieslist, 'factor')){
    specieslist = as.character(paste(specieslist))
  }
  dst = check_filepath(dst)
  
  ## config
  spp.ala.table <- data.frame("Species" = specieslist,
                              "Other_names" = NA,
                              "Synonym"= NA,
                              "New_name" = NA,
                              "Cultivated" = NA,
                              "Records" = NA)
  
  ## progress table
  progress = paste0(dst, 'APC_Progress_table.RData')
  save(spp.ala.table, file = progress)
  
  ## native records dir
  native = paste0(dst, 'ALA_Records_of_Non_cultivated_APC_Spp')
  if(dir.exists(native)){
    print('Warning: directory exists and contents will be written over...')
  } else {
    dir.create(native)
  }
  
  ## raw downloads dir
  raw_files = paste0(dst, 'raw_files')
  if(dir.exists(raw_files)){
    print('Warning: directory exists and contents will be written over...')
  } else {
    dir.create(raw_files)
  }
  
  ## default args
  download_args = list(
    method = "indexed", 
    download_reason_id = 7,
    dst = raw_files,
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
  } # end 
  
  ## fails container
  fails = NULL
  
  ## other names container
  possible_synonyms = list()
  
  ## fails, and synonyms dst paths
  fails_fp = paste0(dst, 'failed_downloads.RData')
  possible_synonyms_fp = paste0(dst, 'possible_synonyms.RData')
  
  flipper_idx <<- 1
  
  ## run
  try(repeat{
    
    load(progress)
    not.done <- which(is.na(spp.ala.table$Records))
    spp      <- specieslist[not.done]
    
    ## Break check
    if(length(spp)==0) break 
    
    ## Not finished?... which species do we still need to work through?
    
    counter = function(flipper_idx){
      #if(flipper_idx == 6) {
      #  flipper_idx = 1
      #}
      #flipper = c('|', '/', '_', '\\', '|' )
      msg = paste0('downloading taxa ' flipper_idx)
      #cat("\r", flipper[flipper_idx])
      cat("\r", msg)
      flipper_idx <<- flipper_idx + 1
    }
    counter(flipper_idx)
    warnings()
    ## LOOP!
    for(sp.n in spp){
      
      # sp.n = spp[2]
      
      ## Match the species to a row in the table
      R.id = match(sp.n, as.character(spp.ala.table$Species)) 
      #match(sp.n, gsub("/","",as.character(spp.ala.table$Species)))
      
      ## 1. First, download all records from ALA 
      ## Some records in the ALA are blank which causes an error. 
      ## Get around this by passing the error to 'check' and calling 
      ## the error '-2.' If the return from ALA is good, check #3 will 
      ## inheret the usual listed 'occurrences' and 'list' objects 
      ## (e.g class = list; length = 2)
      
      #sp.n = 'Abelmoschus manihot'
      #length(ala)
      #ala$data
      #check <- tryCatch(ala <- occurrences(taxon = sp.n ,download_reason_id=7),
      #                  error=function(e) -2)
      check <- tryCatch(ala <- do.call(download_occurrences, 
                                       c(list(taxon = sp.n), download_args)),
                        error = function(e) -2)
      
      #if(length(check) == 2){
      if(is(check, 'data.frame')){  
        #if(dim(ala$data)[1] != 0){
        if(dim(ala)[1] != 0){
          
          if(!file.exists(paste0(native, '/', sp.n, '_ALA_records.RData'))){
            
            ## Flag records that do not have the species name asked for 
            ## (suggesting ALA has determined this is a synonym and redirected
            ## me),or, multiple species names (not inlcuding blanks and NA's).
            #if(any(sp.n %in% ala$data$species)){
            
            ## 0. check col names for expected species field:
            check_names = grep('scientificName', names(ala))
            ## if fail, throw error for now:
            if(!length(check_names)){
              stop('The field names are not as expected - needs addressing')
            }
            if(any(sp.n %in% ala$scientificName)){
              
              ## 1. The species name matches, so presumably is not a synonym, 
              ## but are there other names as well?
              spp.names = as.character(unique(c(
                sp.n, ala$scientificName)))
              if( length(spp.names) > 1 ){ 
                # Enter data in progress table
                spp.ala.table$Other_names[R.id] = 1
                possible_synonyms$sp.n = spp.names 
              }
              
              ## 2. Dicth all of those that are in Australia but are
              ## noted as being "Cultivated..Escapee"
              cultivated = names(ala)[grep('Cultivated', names(ala))]
              ## coerce to R logical
              ala[[cultivated]] = as.logical(paste(ala[[cultivated]]))
              if(any(ala[[cultivated]])){
                ala = ala[!ala[[cultivated]],]
              }
              
              ## 3. Save.
              save(ala, file=paste0(native, '/', sp.n, '_ALA_records.RData'))
              spp.ala.table$Records[R.id] = nrow(ala)
              
            } else {
              ## Species name does not match, 
              ## so presumably the recorded one is an out of date synonym
              spp.names = as.character(unique(c(ala$scientificName, sp.n)))
              # Enter data in progress table
              spp.ala.table$Synonym[R.id] = 1
              spp.ala.table$Records[R.id] = -1
              if(length(spp.names)){
                spp.ala.table$New_name[R.id] = paste(spp.names[1])
              }
            }
            
          } else {
            # Load existing RData file
            load(paste0(native, '/', sp.n, '_ALA_records.RData'))
            spp.ala.table$Records[R.id] = nrow(ala)
          }
          
        } else {
          spp.ala.table$Records[R.id] = 0
          print(paste('The ALA record was empty for', sp.n))
        } # end if else no records
        
      } else {
        if(check == -2){
          spp.ala.table$Records[R.id] = check
          print(paste('No data was available from the ALA for', sp.n))
          fails = c(fails, sp.n)
        } else {
          stop('Error - no records, and the error was not -2...')
        }
      } # end if else try
      
      # Save progress table
      save(spp.ala.table, file = progress)
      # Tidy up
      if(exists('R.id')) rm(R.id)
      if(exists('ala')) rm(ala)
      if(exists('spp.names')) rm(spp.names)
      # Clean out
      gc()
      # Pause
      Sys.sleep(5)
    } # Loop through spp
    
  })
  
  save(possible_synonyms, file = possible_synonyms_fp)
  save(fails, file = fails_fp)
  
  ## summary
  search_summary = paste('Found records for ', 
                         nrow(spp.ala.table) - length(fails), 
                         ' of ', nrow(spp.ala.table), ' listed taxa')
  cat('Completed occurrence searches', sep = '\n')
  cat(search_summary)
  
}




## here's the APC data
root = 'Z:/users/bus07c/GBACC/Evo_GDM_Gaussian_Niche/Input/Plants'
APC <- read.csv(paste(root, '/APC_Unique_Native_Spp.csv', sep = ''))
APC <- sort(unique(c(paste( APC$Genus," ", APC$Species,sep = ''))))

dst = 'O:/source/biol/vascular_plants'
ptm = proc.time()
download_taxalist(APC, dst = dst, parallel = TRUE)
proc.time()-ptm


src = paste0(dst, '/raw_files')
