## ----------------------------------------------------------------------------
##		   
## ALA download by species list 
##
## Plant code (becuase there are lots of plants, and lots of weedy plants)
##
## Adapted from Alex Bush's code:
## (Z:\users\bus07c\GBACC\Evo_GDM_Gaussian_Niche\Script)
##
## ----------------------------------------------------------------------------


##
library(ALA4R)
ala_config(caching="on")
source('Z:/users/bitbucket/utils/R/utils_cw.R')
library(stringr)
library(RCurl)
library(assertthat)

## get latest APNI species list
getAPC = function(){
  ## Download latest APC species list and filter
  url = 'http://www.cpbr.gov.au/apni/apni.html'
  assert_that(url.exists(url))
  get_apc_webpage = readLines(url)
  ## filter out non-species components
  filter_list = grep('taxon_id', get_apc_webpage)
  apc_spp = get_apc_webpage[filter_list[1]:filter_list[length(filter_list)]]
  ## strip html
  apc_spp = gsub("<.*?>", "", get_apc_spp)
  ## apply some filters:
  ## - rm anything after second white space (i.e. below species)
  apc_spp = word(apc_spp, 1, 2)
  ## - rm any record with a '.' --> might be too heavy handed?
  check_dot = grep('\\.', apc_spp)
  check_com = grep(',', apc_spp)
  checks = c(check_dot, check_com)
  if(length(checks)){
    apc_spp = apc_spp[-c(checks)]
  }
  ## - rm anything with only one character
  check_char = strsplit(apc_spp, ' ')
  check_char = lapply(check_char, function(x) nchar(x[2]))
  apc_spp = apc_spp[-c(which(check_char == 1))]
  ## - rm species with uppercase
  check_upper = word(apc_spp, 2) 
  check_upper = substr(check_upper, 1, 1)
  check_upper = unname(sapply(check_upper, function(x) x %in% LETTERS))
  apc_spp = apc_spp[!check_upper]
  ## - rm duplicates
  apc_spp = unique(apc_spp)
  ## format
  apc_spp = data.frame(spp = apc_spp, 
                       genus = word(apc_spp, 1),
                       species = word(apc_spp, 2), row.names = NULL)
  apc_spp = apply(apc_spp, 2, function(x) as.character(paste(x))) 
  apc_spp = apply(apc_spp, 2, function(x) as.character(paste(x)))
  ## done
  return(apc_spp)
} 


this_url = 'http://biocache.ala.org.au/ws/occurrences/index/download?q=Abelmoschus+manihot&reasonTypeId=7'

APC[1]
df = download_occurrences(APC[1], method = "indexed", 
                          download_reason_id = 7, verbose = TRUE)

library(RCurl)
library(httr)

download_specieslist = function(specieslist, dst, ...){
  stop('new version of this in download_ala.R')
  ## args:
  ## specieslist: vector of species to download records for
  ## dst: a filepath to write to
  ## todo:
  ## add in start-fresh or overwrite flag
  
  ## checks
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
  
  ## run
  try(repeat{
    
    load(progress)
    not.done <- which(is.na(spp.ala.table$Records))
    spp      <- specieslist[not.done]
    
    ## Break check
    if(length(spp)==0) break 
    
    ## Not finished?... which species do we still need to work through?
    
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

ptm = proc.time()
download_specieslist(APC, dst)
proc.time()-ptm



## params:
## ----------------------------------------------------------------------------
dst = 'O:/source/biol/vascular_plants'
APC <- read.csv(paste(ab_root, '/APC_Unique_Native_Spp.csv', sep = ''))
APC <- sort(unique(c(paste( APC$Genus," ", APC$Species,sep = '')))) 
## These species are trouble.... boot them:
## spp.ala.table[9863:9865,1]
## Gompholobium sp. Dave's Creek (P.I.Forster+ PIF5979)
## Gompholobium sp. Pilbara (N.F.Norris 908)           
## Gompholobium sp. Point Archer J.Wrigley+ NQ1301) 
APC = APC[-c(9863:9865)]

## run:
## ----------------------------------------------------------------------------
download_specieslist(APC, dst)

x = c(x, c('|', '/', '_', '\\', '|' ))
for (i in 1:5) cat(format(x[i], justify = 'left'))

for (i in 1:10) {Sys.sleep(0.5); cat("\r",x[i])}
  warnings()



xx = read.table(ll[15])
library(gdmEngine)

ala_config()$user_agent

download_to_file
thisfile <- cached_get(url = this_url, type = "binary_filename", 
                      verbose = TRUE)

library(dplyr)
tmp = tempfile()
response = GET(this_url, write_disk(tmp, overwrite=TRUE))
if(response$status_code != 200) cat(response$status_code)
dump_into = tempdir()
## assumes data.csv is always first...
df <- unzip(tmp, exdir = dump_into)[1] %>% read.csv()


cached_get=function(url,type="text",caching=ala_config()$caching,verbose=ala_config()$verbose,on_redirect=NULL,on_client_error=NULL,on_server_error=NULL,encoding=ala_config()$text_encoding) {
  assert_that(is.notempty.string(url))
  assert_that(is.string(type))
  type=match.arg(tolower(type),c("text","json","filename","binary_filename"))
  assert_that(is.string(caching))
  caching=match.arg(tolower(caching),c("on","off","refresh"))
  assert_that(is.flag(verbose))
  
  ## strip newlines or multiple spaces from url: these seem to cause unexpected behaviour
  url=str_replace_all(url,"[\r\n ]+"," ")
  if (nchar(url)>getOption("ALA4R_server_config")$server_max_url_length) warning("URL length may be longer than is allowed by the server")
  
  if (identical(caching,"off") && !(type %in% c("filename","binary_filename"))) {
    ## if we are not caching, get this directly without saving to file at all
    if (verbose) { cat(sprintf("  GETting URL %s\n",url)) }
    
    ## if use RCurl directly
    h=basicHeaderGatherer()
    x=getURL(url=url,useragent=ala_config()$user_agent,header=FALSE,headerfunction=h$update)
    diag_message=""
    if ((substr(h$value()[["status"]],1,1)=="5") || (substr(h$value()[["status"]],1,1)=="4")) {
      ## do we have any useful diagnostic info in x?
      diag_message=get_diag_message(x)
    }
    check_status_code(h$value()[["status"]],extra_info=diag_message,on_redirect=on_redirect,on_client_error=on_client_error,on_server_error=on_server_error)
    
    ## else use httr
    ##x=GET(url=url,user_agent(ala_config()$user_agent)) ## use httr's GET wrapper around RCurl
    ##check_status_code(x,on_redirect=on_redirect,on_client_error=on_client_error,on_server_error=on_server_error)
    ##x=content(x,as="text")
    
    if (identical(type,"json")) {
      if (nchar(x)<1) {
        ## empty string, fromJSON will throw error, so just return NULL
        x=NULL
      } else {
        x=jsonlite::fromJSON(x) ## do text-then-conversion, rather than content(as="parsed") to avoid httr's default use of RJSONIO::fromJSON
      }
    }
    x
  } else {
    ## use caching
    thisfile=download_to_file(url,binary_file=identical(type,"binary_filename"),verbose=verbose,on_redirect=on_redirect,on_client_error=on_client_error,on_server_error=on_server_error)
    if (!file.exists(thisfile)) {
      ## file does not exist
      NULL
    } else {
      if (type %in% c("json","text")) {
        if (!(file.info(thisfile)$size>0)) {
          NULL
        } else {
          ## for json, previously we read as string first, then passed this to fromJSON. This was compatible with either RJSON's or rjsonlite's version of fromJSON (only RJSON's version can take a filename directly). However, it introduced issues with readLines and handling of encoding. Now we pass the file directly to jsonlite::fromJSON
          if (type=="json") {
            ## convert directly from file - this also allows jsonlite to handle encoding issues
            jsonlite::fromJSON(thisfile)
          } else {
            fid=file(thisfile, "rt")
            out=readLines(fid,warn=FALSE,encoding=encoding)
            close(fid)
            out
          }
        }
      } else if (type %in% c("filename","binary_filename")) {
        thisfile
      } else {
        ## should not be here! did we add an allowed type to the arguments without adding handler code down here?
        stop(sprintf("unrecognized type %s in cached_get",type))
      }
    }
  }
}                      
                      
                      

### probably delete below...
# Alex stored relevant data here:
ab_root <- 'Z:/users/bus07c/GBACC/Evo_GDM_Gaussian_Niche/Input/Plants'
# Put my data here:
cw_root <- 'Z:/users/war42q/ALA/ala_downloads/plants'
if(dir.exists(cw_root)) print(paste('output data goes here:', cw_root))

# Load APC table and create a vector of names
APC <- read.csv(paste(ab_root, '/APC_Unique_Native_Spp.csv', sep = ''))
APC <- sort(unique(c(paste( APC$Genus," ", APC$Species,sep = '')))) # vector of 22135 spp names

head(APC)

# Set coordinate precision threshold
coord_threshold <- 1500

# run number
run_n = 4
new_dir = paste(cw_root, '/ALA_Records_of_Non_cultivated_APC_Spp', run_n, sep = '')
if(dir.exists(new_dir)){
  print('Warning: new_dir exists and contents will be written over...')
} else {
  dir.create(new_dir)
}

# if not done:
# Create table to store progress of loop
spp.ala.table <- data.frame("Species"=APC,"Synonym"=NA,"New.name"="0","Cultivated"=NA,"Records"=NA)
spp.ala.table$New.name <- as.character(spp.ala.table$New.name)
save(spp.ala.table, file= paste(cw_root, '/ALA_Records_of_Non_cultivated_APC_Spp', run_n, 
                                '/APC_Progress_table2.RData', sep = ''))

# Create clock?
#clock <- 0
#while(x < 1000){
#  ptm <- proc.time()
#  timer <- (proc.time()-ptm)[3]
# clock <- clock + timer
#}

# custom assertion args
assertion_args <- read.csv('Z:/users/war42q/ALA/assertion_args.csv', h = TRUE)


# These species are trouble.... boot them:
# spp.ala.table[9863:9865,1]
# Gompholobium sp. Dave's Creek (P.I.Forster+ PIF5979)
# Gompholobium sp. Pilbara (N.F.Norris 908)           
# Gompholobium sp. Point Archer J.Wrigley+ NQ1301) 
spp.ala.table[9863:9865,5] <- -1
head(spp.ala.table)
tail(APC)
tail(spp.ala.table)


#######
# Possible fix for embedded nul stuff
# file <- "file.csv"
# tt <- tempfile()  # or tempfile(tmpdir="/dev/shm")
# system(paste0("tr < ", file, " -d '\\000' >", tt))
# fread(tt)
#######


# DO
try(repeat{
  
  load(paste(cw_root, "/ALA_Records_of_Non_cultivated_APC_Spp", run_n, "/APC_Progress_table2.RData", sep = ''))
  not.done <- which(is.na(spp.ala.table$Records))
  spp      <- APC[not.done]
  
  # Break check
  if(length(spp)==0) stop() 
  
  # Not finished?... which species do we still need to work through?
  
  # LOOP!
  for(sp.n in spp){
    # sp.n = spp[2]
    
    # Match the species to a row in the table
    R.id = match(sp.n, as.character(spp.ala.table$Species)) #match(sp.n, gsub("/","",as.character(spp.ala.table$Species)))
    
    # 1. First, download all records from ALA 
    # Some records in the ALA are blank which causes an error. Get around this by
    # passing the error to 'check' and calling ther error '-2.' If the return from ALA is good, check
    # will inheret the usual listed 'occurrences' and 'list' objects (e.g class = list; length = 2)
    check <- tryCatch(ala <- occurrences(taxon = sp.n ,download_reason_id=7), error=function(e) -2)
    
    if(length(check) == 2){
      
      if(dim(ala$data)[1] != 0){
        
        if(!file.exists(paste(cw_root,"/ALA_Records_of_Non_cultivated_APC_Spp", run_n, "/",sp.n,"_ALA_records2.RData",sep=""))){
          
          # Flag records that do not have the species name asked for (suggesting ALA has determined this is a synonym and redirected me),
          # or, multiple species names (not inlcuding blanks and NA's).
          if(any(sp.n %in% ala$data$species)){
            
            # The species name matches, so presumably is not a synonym, but are there other names as well?
            spp.names = as.character(unique(c(sp.n,ala$species)))
            # Remove any NAs or blanks
            spp.names = spp.names[nchar(spp.names)>2]
            if( length(spp.names) > 1 ){ 
              # Enter data in progress table
              spp.ala.table$New.name[R.id] = paste(spp.names)
            }
            
            # 2. Second set aside all of those that are in Australia but are noted as being "occCultivatedEscapee"
            if(any(ala$data$occCultivatedEscapee)){
              cultivated_records = ala$dat[which(ala$data$occCultivatedEscapee),]
              ala$data = ala$data[!ala$data$occCultivatedEscapee,]
            }
            
            # 3. Third, remove records missing Lat/Long information
            ala$data = ala$data[!is.na(ala$data$longitude) & !is.na(ala$data$latitude),]
            
            # 4. Fourth, remove records checks indicate are "fatal" (suspicious).
            #xa = check_assertions(ala)
            #x_afcols = names(ala$data) %in% xa$occurColnames[xa$fatal]                      # Columns of x corresponding to a fatal assertion
            #x_afrows = if(sum(x_afcols)>1){which(apply(ala$data[,x_afcols],1,any))} else {  # Rows of x that have a fatal assertion
            #  ifelse(sum(x_afcols)==1, which(ala$data[,x_afcols]),NA)}
            #if(!any(is.na(x_afrows))){ ala$data = ala$data[-x_afrows,] }                    # Reduce to the "clean" data (data rows without fatal assertions)
            
            # 4.1 Add in custom QA filtering
            ala$data = filterQA(ala$data, assertion_args, FALSE)
            
            
            # 5. Fifth, convert object to table
            ala = ala$data
            
            if(nrow(ala)==0){
              #     cat(paste(sp.n," records are all dubious.\n",sep=""),file=namescript) 
              spp.ala.table$Records[R.id] = 0
              
            } else {
              
              # 6. Sixth, restrict this dataset to only those from Australia
              ala = ala[ala$longitude >= 112.9 & ala$longitude <= 154 & ala$latitude >= -43.74 & ala$latitude <= -9, ]
              
              # 7. Seventh, retain only records with reasonable (<1000m) spatial uncertainty (still includes NA's)
              ala = ala[ala$coordinateUncertaintyInMetres <= coord_threshold | is.na(ala$coordinateUncertaintyInMetres), ]
              
              if(nrow(ala)==0){ 
                #         cat(paste(sp.n," has no records with good precision in Australia.\n",sep=""),file=namescript) 
                spp.ala.table$Records[R.id] = 0
              } else {
                # Remove forward slashes from species name
                sp.n = gsub("/","",sp.n)
                
                # 9. Ninth, save.
                save(ala, file=paste(cw_root, "/ALA_Records_of_Non_cultivated_APC_Spp", run_n, "/",sp.n,"_ALA_records2.RData",sep=""))
                spp.ala.table$Records[R.id] = nrow(ala)
                
                # 10. Tenth, if there were extra cultivated records, save them too.
                if(exists("cultivated_records")){
                  print(paste(sp.n," has ",sum(!is.na(cultivated_records$longitude))," cultivated records with latitude-longitude.",sep=""))
                  cultivated_records = cultivated_records[cultivated_records$coordinateUncertaintyInMetres <= 1000 | is.na(cultivated_records$coordinateUncertaintyInMetres) & !is.na(cultivated_records$longitude) & !is.na(cultivated_records$latitude) & cultivated_records$longitude >= 112.9 & cultivated_records$longitude <= 154 & cultivated_records$latitude >= -43.74 & cultivated_records$latitude <= -9, ]
                  if(nrow(cultivated_records)>0){ save(cultivated_records, file=paste(cw_root, "/ALA_Records_of_Non_cultivated_APC_Spp", run_n, "/",sp.n,"_ALA_cultivated_records2.RData",sep="")) ; spp.ala.table$Cultivated[R.id] = nrow(cultivated_records) }
                  rm(cultivated_records)
                }
                
              } # End of exit from point 6 & 7
              
            } # End of exit from point 5.
            
          } else {
            # Species name does not match, so presumably the recorded one is an out of date synonym
            spp.names = as.character(unique(c(ala$data$species, sp.n)))
            # Remove any NAs or blanks
            spp.names = spp.names[nchar(spp.names) > 2]
            
            # Enter data in progress table
            spp.ala.table$Synonym[R.id] = 1
            spp.ala.table$Records[R.id] = -1
            spp.ala.table$New.name[R.id] = paste(spp.names)
            
          }
          
        } else {
          # Load existing RData file
          load(paste(cw_root, "/ALA_Records_of_Non_cultivated_APC_Spp", run_n, "/",sp.n,"_ALA_records2.RData",sep=""))
          spp.ala.table$Records[R.id] = nrow(ala)
          rm(ala)
        }
        
      } else {
        spp.ala.table$Records[R.id] = 0
        print(paste('The ALA record was empty for', sp.n))
      } # end if else no records
      
    } else {
      if(check == -2){
        spp.ala.table$Records[R.id] = check
        print(paste('No data was available from the ALA for', sp.n))
      } else {
        stop('Error - ALA did not return records, and the error was not -2...')
      }
    } # end if else try
    
    # Save progress table
    save(spp.ala.table, file=paste(cw_root, "/ALA_Records_of_Non_cultivated_APC_Spp", run_n, "/APC_Progress_table2.RData", sep = ''))
    # Tidy up
    if(exists('R.id')) rm(R.id)
    if(exists('ala')) rm(ala)
    if(exists('xa')) rm(xa)
    if(exists('x_afrows')) rm(x_afrows)
    if(exists('x_afcols')) rm(x_afcols)
    if(exists('spp.names')) rm(spp.names)
    # Clean out
    gc()
    # Pause
    Sys.sleep(5)
  } # Loop through spp
  
}) # End of repeat
# -->> got to Mastigolejeunea humilis_ALA_records2
# ---->> ugh corrupt RData
# -->> start again





# Run: 1.
# First run got to Acacia salicina. Bumped up the pause here to 1sec as there was an SO suggestion that this would 
# dodge the cryptic error 'Error in function (type, msg, asError = TRUE)  : 
# transfer closed with outstanding read data remaining

# Run: 2.
# Got to Actephila longipedicellata

# Run: 3.
# Got to Alyogyne pinoniana



# - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = #
#               COMBINE
# - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = #

rm(list=ls())
library(raster)
# Load APC table and create a vector of names
APC = read.csv("Input/Plants/APC_Unique_Native_Spp.csv")
APC = sort(unique(c(paste( APC$Genus," ",APC$Species,sep=""))))

# Load the latest version of the plant data
load(paste("Input/Plants/Plant_Thermal_Niche_Width_Table_Step8.RData",sep=""))
nursery.plants = unique(as.character(plant.TNW$Species))
# Make a second list of species in the APC for which I do not already have plant records as part of the nurseries list
APC2 = APC[is.na(match(APC, nursery.plants))]

# Load a 1km grid of Australia
map = raster("Input/Aus_1km_basemap.asc")

# Simple table of details needed from ALA
plantmat = data.frame("Species"="S","Longitude"=0,"Latitude"=0,"CellID"=0)
cultimat = data.frame("Species"="S","Longitude"=0,"Latitude"=0,"CellID"=0)

for(sp.n in APC2){
  # Load ala records
  if(file.exists(paste("Input/Plants/ALA_Records_of_Non_cultivated_APC_Spp/",sp.n,"_ALA_records.RData",sep=""))){
    load(paste("Input/Plants/ALA_Records_of_Non_cultivated_APC_Spp/",sp.n,"_ALA_records.RData",sep=""))
    ala   = ala[,c("longitude","latitude")]
    cells = cellFromXY(map, ala)
    ala   = cbind(rep(sp.n,nrow(ala)), ala, cells)
    ala   = as.data.frame(ala)
    names(ala) = c("Species","Longitude","Latitude","CellID")
    plantmat = rbind(plantmat, ala)
  }
  
  #   # If they exist also add on the cultivated records
  #   if(file.exists(paste("Input/Plants/ALA_Records_of_Non_cultivated_APC_Spp/",sp.n,"_cultivated_records_ALA_records.RData",sep=""))){
  #     load(paste("Input/Plants/ALA_Records_of_Non_cultivated_APC_Spp/",sp.n,"_cultivated_records_ALA_records.RData",sep=""))
  #     ala   = ala[,c("longitude","latitude")]
  #     cells = cellFromXY(map, ala)
  #     ala   = cbind(rep(sp.n,nrow(ala)), ala, cells)
  #     ala = as.data.frame(ala)
  #     names(ala) = c("Species","Longitude","Latitude","CellID")
  #     cultimat = rbind(cultimat, ala)
  #   }
  print(sp.n)
}

plantmat = plantmat[-1,]
cultimat = cultimat[-1,]

# Now add records from species that are in the cultivated dataset
for(sp.n in APC){
  # Load ala records
  if(file.exists(paste("Input/Plants/ALA_Records/",sp.n,"_ALA_records.RData",sep=""))){
    load(paste("Input/Plants/ALA_Records/",sp.n,"_ALA_records.RData",sep=""))
    ala   = ala[,c("longitude","latitude")]
    cells = cellFromXY(map, ala)
    ala   = cbind(rep(sp.n,nrow(ala)), ala, cells)
    ala = as.data.frame(ala)
    names(ala) = c("Species","Longitude","Latitude","CellID")
    plantmat = rbind(plantmat, ala)
  }
  print(sp.n)
}

# Some of the cells may not be on land
land = which(!is.na(map[1:length(map)]))
plantmat = plantmat[plantmat$CellID%in%land,]

# Save tables
save(plantmat, file=paste("Input/Plants/Aggregated_ALA_Plant_Records.RData",sep=""))
save(cultimat, file=paste("Input/Plants/Aggregated_ALA_Cultivated_Plant_Records.RData",sep=""))


