###########################################################################################
##
## GDM Workflow functions - to support implementation of the GDM workflow
##
## Karel Mokany et al, (CSIRO) - Draft 30 August 2017
##
###########################################################################################


# Install libraries


# Load libraries
library(ALA4R)
library(raster)
library(data.table)
library(docstring)
library(assertthat)
library(httr)
source('//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/utilities.R')

## push a file to the gdm bit bucket repo
## added more bling docstrings
?gitr.push # e.g. gitr.push('all')

##______________________________________________________________________________________##

download_occurrences <- function (taxon, wkt, fields, qa, method = "offline", 
                                  email, download_reason_id, dry_run = FALSE, 
                                  verbose = TRUE) {
  
  #'@title Download occurrences from the ALA
  #'
  #'@description Simplified and customised ALA occurrence download function. 
  #'Largely a hack of the ALA4R function 'occurrences'
  #'
  #'@param taxon (string) Can be a free text search, or taxonomic rank can be used to focus search.
  #'@param wkt (WKT formatted string) Bounding box to filter geographic locations of occurrences to (optional).
  #'@param fields (CSV) Fields to return. Defaults to subset of fields.
  #'@param qa (CSV) Quality assertion fields. Currently set to 'includeall.'
  #'@param method (string) One of 'offline' or indexed. 'offline' only method available currently.
  #'@param email (string) Email address to receive path to .zip archive.
  #'@param download_reason (int) Valid ALA download reason. See website. Accepts numeric codes only (0:12). 
  #'@param dry_run (bool) Return URL to be requested without sending the request. Default: FALSE
  #'@param verbose (bool) Switch print statements and status checking from the ALA servers on. D
  #'Default: TRUE. 
  #'
  #'@returns if verbose, query status
  #'
  #'@note strips away much of the functionality of the ALA4R function
  #'
  #'@section ToDo: add indexed method
  #'
  #'@section Depends: pkg: assertthat, httr
  #'
  #'@examples 
  #'download_occurrences(taxon="genus:Wandella", wkt="POLYGON((112.0+44.0,154.0+44.0,154.0+9.0,112.0+9.0,112.0+44.0))", 
  #'method = "offline",email="karel.mokany@csiro.au", download_reason_id="7", verbose = TRUE)
  #'
  
  pkg_check = suppressWarnings(lapply(c('httr', 'assertthat'), require, character.only = TRUE))
  if(!all(unlist(pkg_check))) stop('Could not load one of required packages httr or asserthat')
  
  ## possibly leave this in case smaller downloads are required
  method <- match.arg(tolower(method), c("offline", "indexed"))
  if (method == "indexed") {
    valid_fields_type <- "occurrence_stored"
    stop('indexed currently not available: use ALA4R package instead')
  } else {
    valid_fields_type <- "occurrence"
  } 
  
  is.notempty.string = function(x) {
    is.string(x) && !is.na(x) && nchar(x)>0
  }
  
  ## query container
  this_query <- list()
  
  ## must be included
  if (!missing(taxon)) {
    if (is.factor(taxon)) {
      taxon <- as.character(taxon)
    }
    assert_that(is.notempty.string(taxon))
    this_query$q <- taxon
  }
  
  # optional
  if (!missing(wkt)) {
    assert_that(is.notempty.string(wkt))
    this_query$wkt <- wkt
  }
  
  ## ammend if I leave out the above
  if (length(this_query) == 0) {
    stop("invalid request: need at least one of taxon, fq, or wkt to be specified")
  }
  
  ## leave for now
  if (method == "offline") {
    #if (record_count_only) 
    #    stop("record_count_only can only be used with method=\"indexed\"")
    if (missing(email) || !is.string(email) || nchar(email) < 1) 
      stop("email is required for method=offline")
  }
  
  ## here I think I can just assume that a numeric reason is always entered.
  ## could even hard-code one.
  reason_ok <- !is.na(download_reason_id)
  if (reason_ok) {
    ## Should take a code from 0-12
    reason_ok = reason_ok %in% seq(0:12)
  }
  
  ## modify later
  if (!reason_ok) {
    stop("download_reason_id must be a valid reason_id. See...")
  }
  
  ## This guy. I think for our purposes we just need a standard set of 
  ## field names. These could be called from a file, or loaded in the 
  ## function. Assume for now it is loaded into memory and called by 
  ## fields arg.
  if (!missing(fields)) {
    assert_that(is.character(fields))
    this_query$fields <- paste(fields, collapse = ",")
  } else {
    if(verbose) cat('... Requesting default fields', sep = '\n')
    fields = as.character(c('occurrenceID',
                            'kingdom.p',
                            'phylum.p',
                            'classs.p',
                            'order.p',
                            'family.p',
                            'genus.p',
                            'subgenus.p',
                            'specificEpithet',
                            'scientificName.p',
                            'acceptedNameUsage',
                            'taxonConceptID.p',
                            'taxonRank.p',
                            'eventDate.p',
                            'year.p',
                            'decimalLatitude.p',
                            'decimalLongitude.p',
                            'coordinateUncertaintyInMeters.p'))
    this_query$fields <- paste(fields, collapse = ",")
    
  }
  
  ## Think this should error if missing. As above it will probably be a file
  ## that needs to be in memory (or set to a hard-coded path).
  if (!missing(qa)){
    assert_that(is.character(qa))
    this_query$qa <- paste(qa, collapse = ",")
  } else {    
    if(verbose) cat('... Requesting default quality assertion fields', sep = '\n')
    qa = as.character(c('zeroCoordinates',
                        #'nameNotSupplied',
                        'nameNotRecognised',
                        'decimalLatLongCalculationFromEastingNorthingFailed',
                        'zeroLatitude',
                        #'decimalLatLongCalculationFromVerbatimFailed',
                        'coordinatesCentreOfCountry',
                        'processingError',
                        'occCultivatedEscapee',
                        'coordinatesCentreOfStateProvince',
                        'decimalLatLongConversionFailed',
                        'zeroLongitude'))
    qa = 'includeall'
    this_query$qa <- paste(qa, collapse = ",")
  }
  
  if (method == "offline") this_query$email <- email
  this_query$reasonTypeId <- download_reason_id
  #this_query$sourceTypeId <- ala_sourcetypeid()
  this_query$esc <- "\\"
  ## This guy! 
  #this_query$sep <- "\t"
  this_query$sep <- ","
  ## Make the file name something intuitive
  #this_query$file <- "data"
  ## Build file name from taxon search and date
  fn = gsub("[[:punct:]]", "_", taxon)
  fn = paste(fn, Sys.Date(), sep = '_')
  fn = gsub(' ', '_', fn)
  this_query$file <- fn
  
  ## url builder (source utilities_internal.R)
  build_url_from_parts <- function(base_url,path=NULL,query=list()) {
    this_url <- parse_url(base_url)
    this_url$path <- clean_path(this_url$path,path)
    if (length(query)>0) {
      this_url$query <- query
    }
    build_url(this_url)
  }    
  
  ## internal url builder function (source utilities_internal.R)
  clean_path <- function(...,sep="/") {
    path1 <- sapply(list(...),FUN=function(z)paste(z,sep=sep,collapse=sep)) ## collapse individual arguments
    ## workaround to avoid replacing "http://" with "http:/", since this is now used in GUID strings (July 2016)
    path <- paste(path1,sep=sep,collapse=sep) ## paste parts together
    path <- gsub("http://","http:@@",path,fixed=TRUE)
    path <- gsub(paste0("[",sep,"]+"),sep,path) ## remove multiple slashes
    path <- gsub("http:@@","http://",path,fixed=TRUE)
    sub(paste0("^",sep),"",path) ## remove leading slash
  }
  
  ## assign this directly
  base_url = 'http://biocache.ala.org.au/ws/'
  # getOption("ALA4R_server_config")$base_url_biocache    
  
  ## I think this is all good. 
  ## If I want to simply, remove 'indexed' option
  if (method == "indexed"){ 
    this_url <- build_url_from_parts(base_url, c("occurrences", "index", "download"), 
                                     query = this_query)
  } else {
    this_url <- build_url_from_parts(base_url, c("occurrences", "offline", "download"),
                                     query = this_query)
  }
  
  if(dry_run) {
    
    return(this_url)
    
  } else {
    
    ## Think a lot of this can be junked. 
    ## Build in a pinger to the ALA server.
    if (method == "offline") {
      ## send request	
      result = GET(this_url)
      if(verbose){
        if(status_code(result) == 200) {
          msg = 'ALA query is being processed'
          fmat = paste(rep('-', nchar(msg)), collapse = '')
          requestStatus = paste('Request status: ', content(result)$status, sep = '')
          queueSize = paste('Queue size....: ', content(result)$queueSize, sep = '')
          statusCheck = paste('Status check..: ', content(result)$statusUrl, sep = '')
          cat('', fmat, msg, fmat, requestStatus, queueSize, statusCheck, sep = '\n')
          shell.exec(paste(content(result)$statusUrl))
        } else {
          msg = paste(': Something has gone wrong with the request ',
                      '(status ', status_code(result), ')', sep = '')
          errorType = content(result)$errorType
          cat(errorType, msg, sep = '')
        }
      }
    }
  }
}

##______________________________________________________________________________________##
## A function to read in, filter and format data downloaded from the ALA

ALA.download.file <- "J:\\source\\biol\\class_Reptilia\\class_Reptilia_2017-10-04\\class_Reptilia_2017-10-04.csv"

# Clean_ALA_data <- function(ALA.download.file, # path to the downloaded .csv file
#                            output.folder = getwd(), # a folder to save the outputs to
#                            output.name = "filtered_data", # a name to use in saving the outputs 
#                            domain.mask,      # a raster layer specifying the analysis domain
#                            earliest.year = 1970,
#                            spatial.uncertainty.m = 2000)
#)
#  {
## Read in the data
ALA.data <- read.csv(ALA.download.file)  

# catch the original number of records
n.recs.start <- nrow(ALA.data)

## Filter the data by QA assertions
ALA.data<-ALA.data[!ALA.data$Name.not.supplied=='true',]
ALA.data<-ALA.data[!ALA.data$Name.not.recognised=='true',]
ALA.data<-ALA.data[!ALA.data$Cultivated...escapee=='true',]
ALA.data<-ALA.data[!ALA.data$Unparseable.verbatim.coordinates=='true',]
ALA.data<-ALA.data[!ALA.data$Decimal.latitude.longitude.conversion.failed=='true',]
ALA.data<-ALA.data[!ALA.data$Unable.to.convert.UTM.coordinates=='true',]
ALA.data<-ALA.data[!ALA.data$Supplied.coordinates.are.zero=='true',]
ALA.data<-ALA.data[!ALA.data$Zero.latitude=='true',]
ALA.data<-ALA.data[!ALA.data$Zero.longitude=='true',]
ALA.data<-ALA.data[!ALA.data$Coordinates.centre.of.country=='true',]
ALA.data<-ALA.data[!ALA.data$Supplied.coordinates.centre.of.state=='true',]
ALA.data<-ALA.data[!ALA.data$processing.Error=='true',]

## Filter by taxonomic resolution ##ADD AN ARGUMENT TO SPECIFY THIS
ALA.data <- ALA.data[(ALA.data$Taxon.Rank...matched == 'species' | ALA.data$Taxon.Rank...matched == 'subspecies'),]

## Filter by date range 
ALA.data <- ALA.data[(ALA.data$Year...parsed >= earliest.year),]

## Filter by spatial Uncertainty
ALA.data <- ALA.data[(ALA.data$Coordinate.Uncertainty.in.Metres...parsed <= spatial.uncertainty.m),]

## Filter by presence of spatial coordinates (necessary data)
ALA.data <- ALA.data[!is.na(ALA.data$Latitude...processed),]
ALA.data <- ALA.data[!is.na(ALA.data$Longitude...processed),]

## Filter by spatial domain 
ALA.data <- ALA.data[(ALA.data$Longitude...processed > domain.mask@extent@xmin),]
ALA.data <- ALA.data[(ALA.data$Longitude...processed < domain.mask@extent@xmax),]
ALA.data <- ALA.data[(ALA.data$Latitude...processed > domain.mask@extent@ymin),]
ALA.data <- ALA.data[(ALA.data$Latitude...processed < domain.mask@extent@ymax),]

## Filter by location on spatial grid
grd.pts<-extract(domain.mask,cbind(ALA.data$Longitude...processed,ALA.data$Latitude...processed))
ALA.data <- ALA.data[!is.na(grd.pts),]

## Create cleaned data file for modelling + log file of how it was created
# select only the data we will/might use later
ALA.data <- ALA.data[,c(1,10,13,15:18)]

# write the data to file
out.path <- file.path(output.folder,paste0(output.name,"_",Sys.Date(),".csv")) 
write.csv(ALA.data, out.path)

# write a log file describing how the data was created *************************************
fileConn<-file(file.path(output.folder,paste0(output.name,"_",Sys.Date(),"_log_file.txt")))
writeLines("#######################################################################",fileConn)
writeLines("###",fileConn)
writeLines("### ALA data filtration log file ",fileConn)
writeLines("###",fileConn)
writeLines(paste0("### Created ",Sys.time()," using the Clean_ALA_data() function."),fileConn)
writeLines("###",fileConn)
writeLines("#######################################################################",fileConn)
writeLines("",fileConn)
writeLines(paste0("Input data file = ", ALA.download.file),fileConn)
writeLines(paste0("Output data file = ", out.path),fileConn)
writeLines(paste0("Domain mask applied = ", domain.mask),fileConn)
writeLines(paste0("Data restricted to after ",earliest.year),fileConn)
writeLines(paste0("Data restricted to spatial uncertainty < ",spatial.uncertainty.m, "m"),fileConn)
writeLines(paste0("Number of records before filtering = ", n.recs.start),fileConn)
writeLines(paste0("Number of records after filtering = ", nrow(ALA.data)),fileConn)
close(fileConn) #**************************************************************************

#} # end function Clean_ALA_data


    
#----- Below = old code, when we were dealing with multiple files -------#
#  ## Standardise ALA input ######################
#  # Get the name/s of the data file/s in this folder
#  folder.files <- list.files(path = ALA.download.folder)
#  data.files <- folder.files[substr(folder.files,1,4) == "data"]   
#  # If there's only one data file, read it in
#  if(length(data.files) == 1)
#    {
#    ALA.data <- read.csv(paste0(ALA.download.folder,"\\",data.files[1]), na.strings=c("na","NA",""), header=TRUE, sep="\t")
#    } # end if length(data.files) == 1
#  # if there are multiple data files, read them in & join them together
#  if(length(data.files) > 1)
#    {
#    #ALA.data <- read.table(paste0(ALA.download.folder,"\\",data.files[1]), header = TRUE, sep = "\t", na.strings = c("","NaN","NA","nan"), quote = "", fill = TRUE)
#    
#    ALA.data <- read.csv(paste0(ALA.download.folder,"\\",data.files[1]), na.strings=c("na","NA",""), header=TRUE)#, sep="\t")
#    
#    #ALA.data <- fread(paste0(ALA.download.folder,"\\",data.files[1]), sep="\t", na.strings=c("na","NA",""),header=TRUE, data.table=FALSE)
#   
#    #ALA.data <- fread(paste0(ALA.download.folder,"\\",data.files[1]), sep=",", na.strings=c("na","NA",""),header=TRUE, data.table=FALSE)#,sep="\t")
#    
#    for(i in 2:length(data.files))
#      {
## csv      
#next.read <- read.csv(paste0(ALA.download.folder,"\\",data.files[i]), na.strings=c("na","NA",""," "), header=FALSE)#, sep="\t")
#     
## fread
##next.read <- fread(paste0(ALA.download.folder,"\\",data.files[i]),  sep = "\t", header = TRUE, fill = TRUE, blank.lines.skip = TRUE)
##next.read <- fread(paste0(ALA.download.folder,"\\",data.files[i]),  sep = "\t", fill = TRUE, blank.lines.skip = TRUE)
#n.file.rows<-length(readLines(paste0(ALA.download.folder,"\\",data.files[i])))
#???next.read <- fread(paste0(ALA.download.folder,"\\",data.files[i]), sep = "\t", nrows=(n.file.rows-1), header = FALSE, fill = TRUE, blank.lines.skip = TRUE)
#colnames(next.read) <- colnames(ALA.data)
#next.read <- as.data.frame(next.read)
## top.check <- readLines(paste0(ALA.download.folder,"\\",data.files[i]), n=5)

# read.table
#next.read<-read.table(paste0(ALA.download.folder,"\\",data.files[i]), header = FALSE, sep = "\t", na.strings = c("","NaN","NA","nan"))
#NAMES <- read.table(paste0(ALA.download.folder,"\\",data.files[i]), nrow = 1,  sep = "\t")

      ALA.data <- rbind(ALA.data,next.read)
      } # end for i 
    } # end if length(data.files) > 1

  
  


