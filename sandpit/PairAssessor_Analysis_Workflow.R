###########################################################################################
##
## Sitepair sample attribute assessor workflow
##
## Karel Mokany et al, (CSIRO) - Draft 28 April 2018
##
###########################################################################################

# Install libraries


# Load libraries
library(ALA4R)
library(raster)
library(gdmEngine)
library(data.table)
library(foreach)
library(doParallel)
#library(dplyr)
#library(magrittr)
#library(plyr)
#library(assertthat)
#library(spatstat)

## ESTABLISH KEY INPUTS ------------------------------------------------------------------##
# Read in a spatial raster specifying the domain and resolution to be modelled
Aus.domain.mask <- raster("//ces-10-cdc/OSM_CDC_GISDATA_work/AUS0025/CLIM/MASK/MASK0.flt")

# SPECIFY ALA DATA FILTERING THRESHOLDS
data.start.year = 1970
location.uncertainty.limit = 2000

# Specify Environmental layers
climate.files <- list.files(path = "//lw-osm-02-cdc/OSM_CBR_LW_R51141_GPAA_work/ENV/A/OUT/1990", full.names=TRUE, pattern = ".flt")
terrain.files <- list.files(path = "//lw-osm-02-cdc/OSM_CBR_LW_R51141_GPAA_work/ENV/A/OUT/LAND", full.names=TRUE, pattern = ".flt")
env.files <- c(climate.files, terrain.files)
env.files <- env.files[(substr(env.files, nchar(env.files)-3, nchar(env.files)) == ".flt")] # to remove some arcmap filenames
env.files <- env.files[-c(20,21,32,35,36,38,39,40,43,44,45,46)] # remove grids we don't want to assess in the modelling
env.stk <- stack(env.files)

# PLANTS INPUTS
species.names.file <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/vascular_plants/APC_and_Orchid_SpeciesNames.csv"
species.names <- read.csv(species.names.file)
species.names <- as.character(species.names[,1])
species.records.folder <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/vascular_plants"
species.records.folder.raw <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/vascular_plants/raw_files"
data.processing.folder <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants"
agg.cell.rad <- 2.25
min.rich.limit <- 10
max.rich.limit <- 400
min.rich.rad <- 200
min.rich.proportion <- 0.25
n.pairs.model <- 144000 # equates to each site used 10 times
train.proportion <- 0.8
n.pairs.test <- 36000   # equates to each site used 10 times

# AMPHIBIANS INPUTS
species.names.file <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/amphibians/AFD-20171211T130458.csv"
species.names <- read.csv(species.names.file)
species.names <- paste(species.names$GENUS, species.names$SPECIES)
species.names <- unique(species.names)
species.records.folder <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/amphibians"
species.records.folder.raw <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/amphibians/raw_files"
data.processing.folder <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/amphibians"
agg.cell.rad <- 2.25
min.rich.limit <- 2
max.rich.limit <- 50
min.rich.rad <- 200
min.rich.proportion <- 0.25
n.pairs.model <- 100000
train.proportion <- 0.8
n.pairs.test <- 20000

# LAND SNAIL INPUTS
species.names.file <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/land_snails/AusLandSnails_ALASpeciesList_9Mar18.csv"
species.names <- read.csv(species.names.file)
species.names <- species.names$Species.Name
species.names <- unique(species.names)
species.records.folder <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/land_snails"
species.records.folder.raw <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/land_snails/raw_files"
data.processing.folder <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/land_snails"
agg.cell.rad <- 2.25
min.rich.limit <- 2
max.rich.limit <- 50
min.rich.rad <- 50
min.rich.proportion <- 0.25
n.pairs.model <- 50000
train.proportion <- 0.8
n.pairs.test <- 10000


# REPTILE INPUTS
species.names.file <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/reptiles/AFD-20171211T113438.csv"
species.names <- read.csv(species.names.file)
species.names <- paste(species.names$GENUS, species.names$SPECIES)
species.names <- unique(species.names)
species.records.folder <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/reptiles"
species.records.folder.raw <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/source/biol/reptiles/raw_files"
data.processing.folder <- "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/reptiles"
agg.cell.rad <- 2.25
min.rich.limit <- 3
max.rich.limit <- 50
min.rich.rad <- 200
min.rich.proportion <- 0.25
n.pairs.model <- 100000
train.proportion <- 0.8

## DOWNLOAD BIOLOGICAL DATA FROM ALA -----------------------------------------------------##
# Download the species records from ALA
download_taxalist(specieslist = species.names,
                  dst = species.records.folder)

## MERGE THE BIOLOGICAL DATA FOR EACH SPECIES --------------------------------------------##
All.records <- merge_downloads(src=species.records.folder.raw,
                               output.folder = data.processing.folder,
                               parallel = FALSE)

## FILTER THE BIOLOGICAL DATA ------------------------------------------------------------##
Filtered.records <- filter_ALA_data(ALA.download.data = All.records$data,             
                                    output.folder = data.processing.folder,       
                                    domain.mask = Aus.domain.mask,                   
                                    earliest.year = data.start.year,
                                    spatial.uncertainty.m = location.uncertainty.limit)

## AGGREGATE THE BIOLOGICAL DATA TO GRID CELLS -------------------------------------------##
Aggregated.records <- aggregate_ALA_data(ALA.filtered.data = Filtered.records,
                                          domain.mask = Aus.domain.mask,
                                          agg.radius.ncells = agg.cell.rad,
                                          output.folder = data.processing.folder)

## REFINE BIOLOGICAL DATA TO GRID CELLS CONTAINING SUITABLE NUMBERS OF SPECIES RECORDS----##
Selected.records <- select_gridcells_composition(ALA.aggregated.data = Aggregated.records ,
                                                domain.mask = Aus.domain.mask,
                                                min.richness.threshold = min.rich.limit,
                                                max.richness.threshold = max.rich.limit,
                                                reference.radius.ncells = min.rich.rad,
                                                min.proportion.max.richness = min.rich.proportion,
                                                output.folder = data.processing.folder)

## EXTRACT ENVIRONMENTAL DATA FOR SELECTED GRID CELLS ------------------------------------##
Site.Env.Data <- extract_env_data(ALA.composition.data = Selected.records,             
                                  environment.stk = env.stk,
                                  output.folder = data.processing.folder)

##TEMP##
#AMPHIBIANS -------
#Selected.records <- read.csv("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/amphibians/selected_gridcell_composition_2018-03-05.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/amphibians/site_env_data_2018-03-05.csv")
#VASCULAR PLANTS -------
Selected.records <- read.csv("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/selected_gridcell_composition_2018-03-07.csv")
Site.Env.Data <- read.csv("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/site_env_data_2018-03-07.csv")
#LAND SNAILS -------
#Selected.records <- read.csv("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/land_snails/selected_gridcell_composition_2018-03-09.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/land_snails/site_env_data_2018-03-09.csv")
#Reptiles ---------
#Selected.records <- read.csv("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/reptiles/selected_gridcell_composition_2018-03-15.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/reptiles/site_env_data_2018-03-15.csv")
##ENDTEMP##

# Make a single focussed sample ----------------------------------------------------------##
geowt.sample <- sitepair_sample_assessor(site.env.data=Site.Env.Data, 
                                        composition.data=Selected.records,
                                        n.pairs.per.site = 10,
                                        n.crossvalid.tests = 10,
                                        sample.method = 'geowt',
                                        #b.used.factor=2,
                                        #b.dpair.factor=0.5,
                                        #b.epair.factor=1,
                                        #sigma.spair=NULL,
                                        #spair.factor=1,
                                        domain.mask=Aus.domain.mask,
                                        pcs.projargs="+init=epsg:3577", 
                                        bandwidth.geowt=150000,
                                        bandwidth.skip=2,
                                        bandwidth.DistFact=1,
                                        geowt.RndProp=0.05,
                                        output.folder = "//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/amphibians/SitePairSampleAssessment",       
                                        output.name = "Sitepair_sample_assessor_output_geowt1") 


##----------------------------------------------------------------------------------------##
## NOW RUN THE SITEPAIR SAMPLING ATTRIBUTE ASSESSMENT ------------------------------------##
taxa <- "vascular_plants"

#### PARAMETERS FOR RANDOM ####
parameter.tbl <- expand.grid(p.sample.method='random',
                             p.n.pairs.per.site=c(5,10,20), 
                             p.b.used.factor=c(0.1,0.5,1,2,3,4,5),
                             stringsAsFactors=FALSE)
parameter.tbl<-unique(parameter.tbl)
parameter.tbl$run.name<-paste0(taxa,"_",parameter.tbl$p.sample.method,"_",c(1:nrow(parameter.tbl))) ## **CHANGE FOR EACH TAXA**
analysis.out.folder<-paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/",parameter.tbl$p.sample.method[1])  ## **CHANGE FOR EACH TAXA**
write.csv(parameter.tbl,paste0(analysis.out.folder,"/",taxa,"_parameters_",parameter.tbl$p.sample.method[1],".csv"),row.names = FALSE) ## **CHANGE FOR EACH TAXA**
# run the parallel loop over parameter combinations
cl<-makeCluster(25) #setup parallel backend to use 12 processors
registerDoParallel(cl) 
foreach(i.run=1:nrow(parameter.tbl), .packages='gdmEngine') %dopar% {
  this.sample <- sitepair_sample_assessor(site.env.data = Site.Env.Data,
                                          composition.data = Selected.records ,
                                          n.pairs.per.site = parameter.tbl$p.n.pairs.per.site[i.run],
                                          n.crossvalid.tests = 3,
                                          sample.method = parameter.tbl$p.sample.method[i.run],
                                          b.used.factor=parameter.tbl$p.b.used.factor[i.run],
                                          domain.mask=Aus.domain.mask,
                                          pcs.projargs="+init=epsg:3577",
                                          output.folder = analysis.out.folder,
                                          output.name = parameter.tbl$run.name[i.run])
} # end for i.run
stopCluster(cl)

#### PARAMETERS FOR GEODIST ####
parameter.tbl <- expand.grid(p.sample.method='geodist',
                             p.n.pairs.per.site=c(5,10,20), 
                             p.b.used.factor=c(0.1,0.5,1,2,3,4,5),
                             p.b.dpair.factor=c(0.1,0.5,1,2,3,4,5),
                             stringsAsFactors=FALSE)
parameter.tbl<-unique(parameter.tbl)
parameter.tbl$run.name<-paste0(taxa,"_",parameter.tbl$p.sample.method,"_",c(1:nrow(parameter.tbl))) ## **CHANGE FOR EACH TAXA**
analysis.out.folder<-paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/",parameter.tbl$p.sample.method[1])  ## **CHANGE FOR EACH TAXA**
write.csv(parameter.tbl,paste0(analysis.out.folder,"/",taxa,"_parameters_",parameter.tbl$p.sample.method[1],".csv"),row.names = FALSE) ## **CHANGE FOR EACH TAXA**
# run the parallel loop over parameter combinations
cl<-makeCluster(25) #setup parallel backend to use 12 processors
registerDoParallel(cl) 
foreach(i.run=1:nrow(parameter.tbl), .packages='gdmEngine') %dopar% {
  this.sample <- sitepair_sample_assessor(site.env.data = Site.Env.Data,
                                          composition.data = Selected.records ,
                                          n.pairs.per.site = parameter.tbl$p.n.pairs.per.site[i.run],
                                          n.crossvalid.tests = 3,
                                          sample.method = parameter.tbl$p.sample.method[i.run],
                                          b.used.factor=parameter.tbl$p.b.used.factor[i.run],
                                          b.dpair.factor=parameter.tbl$p.b.dpair.factor[i.run],
                                          domain.mask=Aus.domain.mask,
                                          pcs.projargs="+init=epsg:3577",
                                          output.folder = analysis.out.folder,
                                          output.name = parameter.tbl$run.name[i.run])
} # end for i.run
stopCluster(cl)
 
#### PARAMETERS FOR ENVDIST ####
parameter.tbl <- expand.grid(p.sample.method='envdist',
                             p.n.pairs.per.site=c(5,10,20), 
                             p.b.used.factor=c(0.1,0.5,1,2,3,4,5),
                             p.b.epair.factor=c(0.1,0.5,1,2,3,4,5),
                             stringsAsFactors=FALSE)
parameter.tbl<-unique(parameter.tbl)
parameter.tbl$run.name<-paste0(taxa,"_",parameter.tbl$p.sample.method,"_",c(1:nrow(parameter.tbl))) ## **CHANGE FOR EACH TAXA**
analysis.out.folder<-paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/",parameter.tbl$p.sample.method[1])  ## **CHANGE FOR EACH TAXA**
write.csv(parameter.tbl,paste0(analysis.out.folder,"/",taxa,"_parameters_",parameter.tbl$p.sample.method[1],".csv"),row.names = FALSE) ## **CHANGE FOR EACH TAXA**
# run the parallel loop over parameter combinations
cl<-makeCluster(25) #setup parallel backend to use 12 processors
registerDoParallel(cl) 
foreach(i.run=1:nrow(parameter.tbl), .packages='gdmEngine') %dopar% {
  this.sample <- sitepair_sample_assessor(site.env.data = Site.Env.Data,
                                          composition.data = Selected.records ,
                                          n.pairs.per.site = parameter.tbl$p.n.pairs.per.site[i.run],
                                          n.crossvalid.tests = 3,
                                          sample.method = parameter.tbl$p.sample.method[i.run],
                                          b.used.factor=parameter.tbl$p.b.used.factor[i.run],
                                          b.epair.factor=parameter.tbl$p.b.epair.factor[i.run],
                                          domain.mask=Aus.domain.mask,
                                          pcs.projargs="+init=epsg:3577",
                                          output.folder = analysis.out.folder,
                                          output.name = parameter.tbl$run.name[i.run])
} # end for i.run
stopCluster(cl)

#### PARAMETERS FOR GEODENS ####
parameter.tbl <- expand.grid(p.sample.method='geodens',
                             p.n.pairs.per.site=c(5,10,20), 
                             p.b.used.factor=c(0.1,0.5,1,2,3,4,5),
                             p.sigma.spair=c(0.25,0.5,1.0,1.5,2),	
                             p.b.spair.factor=c(0.25,0.5,1,2,4),
                             stringsAsFactors=FALSE)
parameter.tbl<-unique(parameter.tbl)
parameter.tbl$run.name<-paste0(taxa,"_",parameter.tbl$p.sample.method,"_",c(1:nrow(parameter.tbl))) ## **CHANGE FOR EACH TAXA**
analysis.out.folder<-paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/",parameter.tbl$p.sample.method[1])  ## **CHANGE FOR EACH TAXA**
write.csv(parameter.tbl,paste0(analysis.out.folder,"/",taxa,"_parameters_",parameter.tbl$p.sample.method[1],".csv"),row.names = FALSE) ## **CHANGE FOR EACH TAXA**
# run the parallel loop over parameter combinations
cl<-makeCluster(25) #setup parallel backend to use 12 processors
registerDoParallel(cl) 
foreach(i.run=1:nrow(parameter.tbl), .packages='gdmEngine') %dopar% {
  this.sample <- sitepair_sample_assessor(site.env.data = Site.Env.Data,
                                          composition.data = Selected.records ,
                                          n.pairs.per.site = parameter.tbl$p.n.pairs.per.site[i.run],
                                          n.crossvalid.tests = 3,
                                          sample.method = parameter.tbl$p.sample.method[i.run],
                                          b.used.factor=parameter.tbl$p.b.used.factor[i.run],
                                          sigma.spair=parameter.tbl$p.sigma.spair[i.run],
                                          spair.factor=parameter.tbl$p.b.spair.factor[i.run],
                                          domain.mask=Aus.domain.mask,
                                          pcs.projargs="+init=epsg:3577",
                                          output.folder = analysis.out.folder,
                                          output.name = parameter.tbl$run.name[i.run])
} # end for i.run
stopCluster(cl)

#### PARAMETERS FOR GEOWT ####
parameter.tbl <- expand.grid(p.sample.method='geowt',
                             p.n.pairs.per.site=c(5,10,20), 
                             p.bandwidth.geowt=c(100000,150000,200000,300000,400000),
                             p.bandwidth.skip=c(1,2,3,4),
                             p.bandwidth.DistFact=c(0.5,1,1.5,2),
                             p.geowt.RndProp=c(0.01,0.05,0.1,0.2),
                             stringsAsFactors=FALSE)

parameter.tbl<-unique(parameter.tbl)
parameter.tbl$run.name<-paste0(taxa,"_",parameter.tbl$p.sample.method,"_",c(1:nrow(parameter.tbl))) ## **CHANGE FOR EACH TAXA**
analysis.out.folder<-paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/",parameter.tbl$p.sample.method[1])  ## **CHANGE FOR EACH TAXA**
write.csv(parameter.tbl,paste0(analysis.out.folder,"/",taxa,"_parameters_",parameter.tbl$p.sample.method[1],".csv"),row.names = FALSE) ## **CHANGE FOR EACH TAXA**
#  for parallel implementation...
cl<-makeCluster(25) #setup parallel backend to use 14 processors
registerDoParallel(cl) 
# run the parallel loop over parameter combinations
foreach(i.run=1:nrow(parameter.tbl), .packages='gdmEngine') %dopar% {
  this.sample <- sitepair_sample_assessor(site.env.data = Site.Env.Data,
                                          composition.data = Selected.records ,
                                          n.pairs.per.site = parameter.tbl$p.n.pairs.per.site[i.run],
                                          n.crossvalid.tests = 3,
                                          sample.method = parameter.tbl$p.sample.method[i.run],
                                          domain.mask=Aus.domain.mask,
                                          pcs.projargs="+init=epsg:3577",
                                          bandwidth.geowt=parameter.tbl$p.bandwidth.geowt[i.run],
                                          bandwidth.skip=parameter.tbl$p.bandwidth.skip[i.run],
                                          bandwidth.DistFact=parameter.tbl$p.bandwidth.DistFact[i.run],
                                          geowt.RndProp=parameter.tbl$p.geowt.RndProp[i.run],
                                          output.folder = analysis.out.folder,
                                          output.name = parameter.tbl$run.name[i.run])
  } # end for i.run
stopCluster(cl)
##----------------------------------------------------------------------------------------##

##----------------------------------------------------------------------------------------##
## SUMMARISE THE RESULTS FROM THE SITEPAIR SAMPLING ASSESSMENT ---------------------------##
library(gtools)
# Specify the taxa
taxa <- "land_snails"
# Specify the sample types
sample_method <- c("envdist","geodens","geodist","geowt","random")
# Loop through the sample methods, and join the data together
taxa.results.tbl <- NULL
for(i.smp in 1:length(sample_method))
  {
  smp.mth <- sample_method[i.smp]
  # Specify the folder where the output files are held
  analysis.out.folder<-paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/",smp.mth)  ## **CHANGE FOR EACH TAXA**
  analysis.out.files<-as.character(list.files(analysis.out.folder))
  # Read in the parameter table
  parameter.tbl<-read.csv(paste0(analysis.out.folder,"/",taxa,"_parameters_",smp.mth,".csv")) ## **CHANGE FOR EACH TAXA**
  # Create a table to catch the summary data
  results.tbl <- data.frame(ProcTime = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Dissimilarity.Evenness = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Min.Dissimilarity = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Q1.Dissimilarity = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Mdn.Dissimilarity = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Q3.Dissimilarity = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Max.Dissimilarity = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Sites.Geo.Evenness = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Sitepairs.Geo.Evenness = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Min.Sitepairs.GeoDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Q1.Sitepairs.GeoDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Mdn.Sitepairs.GeoDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Q3.Sitepairs.GeoDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Max.Sitepairs.GeoDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Sites.Env.Evenness = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Sitepairs.Env.Evenness = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Min.Sitepairs.EnvDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Q1.Sitepairs.EnvDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Mdn.Sitepairs.EnvDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Q3.Sitepairs.EnvDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Max.Sitepairs.EnvDistance = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Min.nTimes.SitesUsedInPairs = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Q1.nTimes.SitesUsedInPairs = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Mdn.nTimes.SitesUsedInPairs = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Q3.nTimes.SitesUsedInPairs = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Max.nTimes.SitesUsedInPairs = as.numeric(rep(NA,times=nrow(parameter.tbl)))  )
  # Now loop through parameter combos, see if there's an output file for each, read it in and
  # summarise the results, catching in the output table.
  for(i.run in 1:nrow(parameter.tbl))
    {
    #find if there is an output file to match the run.name
    #paste(as.character(parameter.tbl$run.name[i.run]),"_",sep='')
    f.name <- analysis.out.files[pmatch(paste(as.character(parameter.tbl$run.name[i.run]),"_",sep=''), analysis.out.files)]
    if(!is.na(f.name))
      {
      # open the file
      load(file.path(analysis.out.folder,f.name))
      # catch the results 
      results.tbl$ProcTime[i.run] = Sitepair_assessor_results$ProcessingTime[3]
      results.tbl$Dissimilarity.Evenness[i.run] = mean(Sitepair_assessor_results$DissimilarityEvenness)
      results.tbl$Min.Dissimilarity[i.run] = mean(Sitepair_assessor_results$DissimilaritySummary[,1])
      results.tbl$Q1.Dissimilarity[i.run] =  mean(Sitepair_assessor_results$DissimilaritySummary[,2])
      results.tbl$Mdn.Dissimilarity[i.run] =  mean(Sitepair_assessor_results$DissimilaritySummary[,3])
      results.tbl$Q3.Dissimilarity[i.run] =  mean(Sitepair_assessor_results$DissimilaritySummary[,5])
      results.tbl$Max.Dissimilarity[i.run] =  mean(Sitepair_assessor_results$DissimilaritySummary[,6])
      results.tbl$Sites.Geo.Evenness[i.run] = mean(Sitepair_assessor_results$SitesGeoEvenness)
      results.tbl$Sitepairs.Geo.Evenness[i.run] = mean(Sitepair_assessor_results$SitepairsGeoEvenness)
      results.tbl$Min.Sitepairs.GeoDistance[i.run] = mean(Sitepair_assessor_results$SitepairsGeoDistanceSummary[,1])
      results.tbl$Q1.Sitepairs.GeoDistance[i.run] = mean(Sitepair_assessor_results$SitepairsGeoDistanceSummary[,2])
      results.tbl$Mdn.Sitepairs.GeoDistance[i.run] = mean(Sitepair_assessor_results$SitepairsGeoDistanceSummary[,3])
      results.tbl$Q3.Sitepairs.GeoDistance[i.run] = mean(Sitepair_assessor_results$SitepairsGeoDistanceSummary[,5])
      results.tbl$Max.Sitepairs.GeoDistance[i.run] = mean(Sitepair_assessor_results$SitepairsGeoDistanceSummary[,6])
      results.tbl$Sites.Env.Evenness[i.run] = mean(Sitepair_assessor_results$SitesEnvEvenness)
      results.tbl$Sitepairs.Env.Evenness[i.run] = mean(Sitepair_assessor_results$SitepairsEnvEvenness)
      results.tbl$Min.Sitepairs.EnvDistance[i.run] = mean(Sitepair_assessor_results$SitepairsEnvDistanceSummary[,1])
      results.tbl$Q1.Sitepairs.EnvDistance[i.run] = mean(Sitepair_assessor_results$SitepairsEnvDistanceSummary[,2])
      results.tbl$Mdn.Sitepairs.EnvDistance[i.run] = mean(Sitepair_assessor_results$SitepairsEnvDistanceSummary[,3])
      results.tbl$Q3.Sitepairs.EnvDistance[i.run] = mean(Sitepair_assessor_results$SitepairsEnvDistanceSummary[,5])
      results.tbl$Max.Sitepairs.EnvDistance[i.run] = mean(Sitepair_assessor_results$SitepairsEnvDistanceSummary[,6])
      results.tbl$Min.nTimes.SitesUsedInPairs[i.run] = mean(Sitepair_assessor_results$nTimesSitesUsedInPairs[,1])
      results.tbl$Q1.nTimes.SitesUsedInPairs[i.run] = mean(Sitepair_assessor_results$nTimesSitesUsedInPairs[,2])
      results.tbl$Mdn.nTimes.SitesUsedInPairs[i.run] = mean(Sitepair_assessor_results$nTimesSitesUsedInPairs[,3])
      results.tbl$Q3.nTimes.SitesUsedInPairs[i.run] = mean(Sitepair_assessor_results$nTimesSitesUsedInPairs[,5])
      results.tbl$Max.nTimes.SitesUsedInPairs[i.run] = mean(Sitepair_assessor_results$nTimesSitesUsedInPairs[,6])
      # remove the results
      rm(Sitepair_assessor_results)
      }# end if !is.na(f.name)
    } # end for i.run
  parameter.results.tbl <- cbind(parameter.tbl, results.tbl) 
  taxa.results.tbl <- smartbind(taxa.results.tbl,parameter.results.tbl)
  }# end for i.smp

parameter.results.tbl<-taxa.results.tbl[!is.na(taxa.results.tbl$run.name),]
write.csv(parameter.results.tbl, paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/",taxa,"_SampleAssess_Results.csv"), row.names = FALSE)

## Plot results ###############################################
library(ggplot2)
#taxa <- "amphibians" 
#parameter.results.tbl<-read.csv(paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/",taxa,"_SampleAssess_Results.csv"))
figures.out.folder<-paste0("//osm-23-cdc.it.csiro.au/OSM_CBR_LW_DEE_work/processing/biol/",taxa,"/SitePairSampleAssessment/Figures")
# convert factors to factors
parameter.results.tbl$p.sample.method<-as.factor(parameter.results.tbl$p.sample.method)
parameter.results.tbl$p.n.pairs.per.site<-as.factor(parameter.results.tbl$p.n.pairs.per.site)
parameter.results.tbl$p.b.used.factor<-as.factor(parameter.results.tbl$p.b.used.factor)
parameter.results.tbl$p.b.epair.factor<-as.factor(parameter.results.tbl$p.b.epair.factor)
parameter.results.tbl$p.sigma.spair<-as.factor(parameter.results.tbl$p.sigma.spair)
parameter.results.tbl$p.b.spair.factor<-as.factor(parameter.results.tbl$p.b.spair.factor)
parameter.results.tbl$p.b.dpair.factor<-as.factor(parameter.results.tbl$p.b.dpair.factor)
parameter.results.tbl$p.bandwidth.geowt<-as.factor(parameter.results.tbl$p.bandwidth.geowt)
parameter.results.tbl$p.bandwidth.skip<-as.factor(parameter.results.tbl$p.bandwidth.skip)
parameter.results.tbl$p.bandwidth.DistFact<-as.factor(parameter.results.tbl$p.bandwidth.DistFact)
parameter.results.tbl$p.geowt.RndProp<-as.factor(parameter.results.tbl$p.geowt.RndProp)

## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## First compare overall sampling methods, to see which performs best over all
p1<-ggplot(parameter.results.tbl, aes(x=p.sample.method, y=Dissimilarity.Evenness)) +
  geom_boxplot() + labs(y = "Dissim.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p2<-ggplot(parameter.results.tbl, aes(x=p.sample.method, y=Mdn.Dissimilarity)) +
  geom_boxplot() + labs(y = "Mdn.Dissim") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p3<-ggplot(parameter.results.tbl, aes(x=p.sample.method, y=Sitepairs.Geo.Evenness)) +
  geom_boxplot() + labs(y = "Geo.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p4<-ggplot(parameter.results.tbl, aes(x=p.sample.method, y=Sitepairs.Env.Evenness)) +
  geom_boxplot() + labs(y = "Env.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p5<-ggplot(parameter.results.tbl, aes(x=p.sample.method, y=Mdn.Sitepairs.EnvDistance)) +
  geom_boxplot() + labs(y = "Mdn.Env.Dist") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p6<-ggplot(parameter.results.tbl, aes(x=p.sample.method, y=Mdn.Sitepairs.GeoDistance)) +
  geom_boxplot() + labs(y = "Mdn.Geo.Dist") 
png(paste0(figures.out.folder,"/SampleMethod_PerformaceSummary.png"),height=1000,width=500)
multiplot(p1, p2, p3, p4, p5, p6, cols=1)
dev.off()#___

## GEOWT ANALYSIS ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## for p.bandwidth.geowt
p1<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.geowt, y=Dissimilarity.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Dissim.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p2<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.geowt, y=Mdn.Dissimilarity, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Dissim") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p3<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.geowt, y=Sitepairs.Geo.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Geo.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p4<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.geowt, y=Sitepairs.Env.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Env.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p5<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.geowt, y=Mdn.Sitepairs.EnvDistance, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Env.Dist") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p6<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.geowt, y=Mdn.Sitepairs.GeoDistance, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Geo.Dist") 
png(paste0(figures.out.folder,"/p_bandwidth_geowt.png"),height=1000,width=500)
multiplot(p1, p2, p3, p4, p5,p6,cols=1)
dev.off()#___

## for p.bandwidth.skip
p1<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.skip, y=Dissimilarity.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Dissim.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p2<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.skip, y=Mdn.Dissimilarity, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Dissim") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p3<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.skip, y=Sitepairs.Geo.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Geo.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p4<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.skip, y=Sitepairs.Env.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Env.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p5<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.skip, y=Mdn.Sitepairs.EnvDistance, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Env.Dist") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p6<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.skip, y=Mdn.Sitepairs.GeoDistance, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Geo.Dist") 
png(paste0(figures.out.folder,"/p_bandwidth_skip.png"),height=1000,width=500)
multiplot(p1, p2, p3, p4, p5,p6,cols=1)
dev.off()#___

## for p.bandwidth.DistFact
p1<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.DistFact, y=Dissimilarity.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Dissim.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p2<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.DistFact, y=Mdn.Dissimilarity, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Dissim") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p3<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.DistFact, y=Sitepairs.Geo.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Geo.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p4<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.DistFact, y=Sitepairs.Env.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Env.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p5<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.DistFact, y=Mdn.Sitepairs.EnvDistance, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Env.Dist") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p6<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.DistFact, y=Mdn.Sitepairs.GeoDistance, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Geo.Dist") 
png(paste0(figures.out.folder,"/p_bandwidth_distFact.png"),height=1000,width=500)
multiplot(p1, p2, p3, p4, p5,p6,cols=1)
dev.off()#___

## for p.geowt.RndProp
p1<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.geowt.RndProp, y=Dissimilarity.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Dissim.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p2<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.geowt.RndProp, y=Mdn.Dissimilarity, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Dissim") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p3<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.geowt.RndProp, y=Sitepairs.Geo.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Geo.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p4<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.geowt.RndProp, y=Sitepairs.Env.Evenness, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Env.Even") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p5<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.geowt.RndProp, y=Mdn.Sitepairs.EnvDistance, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Env.Dist") + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p6<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.geowt.RndProp, y=Mdn.Sitepairs.GeoDistance, fill=p.n.pairs.per.site)) +
  geom_boxplot() + labs(y = "Mdn.Geo.Dist") 
png(paste0(figures.out.folder,"/p_bandwidth_randProp.png"),height=1000,width=500)
multiplot(p1, p2, p3, p4, p5,p6,cols=1)
dev.off()#___

# RUNTIME
p1<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.geowt, y=ProcTime, fill=p.n.pairs.per.site)) +
  geom_boxplot() + ylim(0, 10000) +  theme(legend.position="none")
p2<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.skip, y=ProcTime, fill=p.n.pairs.per.site)) +
  geom_boxplot() + ylim(0, 10000) +  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), legend.position="none")
p3<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.bandwidth.DistFact, y=ProcTime, fill=p.n.pairs.per.site)) +
  geom_boxplot() + ylim(0, 10000) +  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), legend.position="none")
p4<-ggplot(parameter.results.tbl[parameter.results.tbl$p.sample.method == "geowt",], aes(x=p.geowt.RndProp, y=ProcTime, fill=p.n.pairs.per.site)) +
  geom_boxplot() + ylim(0, 10000) +  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), legend.position="none")
png(paste0(figures.out.folder,"/p_bandwidth_ProcTime.png"),height=500,width=1000)
multiplot(p1, p2, p3, p4, cols=4)
dev.off()#___

## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##








## for p.geo
png(paste0(figures.out.folder,"/p_geo_rndMnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.geo, y=rnd.Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_geo_rndMneRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.geo, y=rnd.Mn.eRMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_geo_MnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.geo, y=Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_geo_rndMnMAE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.geo, y=rnd.Mn.MAE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_geo_MnMAE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.geo, y=Mn.MAE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
## for p.n.pairs.train
png(paste0(figures.out.folder,"/p_nPairsTrain_d2FinMod.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.n.pairs.train, y=D2.FinalMod, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_nPairsTrain_rndMnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.n.pairs.train, y=rnd.Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_nPairsTrain_rndMneRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.n.pairs.train, y=rnd.Mn.eRMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_nPairsTrain_MnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.n.pairs.train, y=Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_nPairsTrain_rndMnMAE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.n.pairs.train, y=rnd.Mn.MAE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_nPairsTrain_MnMAE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.n.pairs.train, y=Mn.MAE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
## for p.b.used
png(paste0(figures.out.folder,"/p_bUsedFactor_rndMnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.b.used.factor, y=rnd.Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_bUsedFactor_rndMneRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.b.used.factor, y=rnd.Mn.eRMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_bUsedFactor_MnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl, aes(x=p.b.used.factor, y=Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
## for p.b.dpair.factor with geo.dist method
png(paste0(figures.out.folder,"/p_bdPairFactor_GeoDist_rndMnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.b.dpair.factor),], aes(x=p.b.dpair.factor, y=rnd.Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_bdPairFactor_GeoDist_rndMneRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.b.dpair.factor),], aes(x=p.b.dpair.factor, y=rnd.Mn.eRMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_bdPairFactor_GeoDist_MnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.b.dpair.factor),], aes(x=p.b.dpair.factor, y=Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
## for p.b.epair.factor with env dist method
png(paste0(figures.out.folder,"/p_bePairFactor_EnvDist_rndMnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.b.epair.factor),], aes(x=p.b.epair.factor, y=rnd.Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_bePairFactor_EnvDist_rndMneRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.b.epair.factor),], aes(x=p.b.epair.factor, y=rnd.Mn.eRMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
png(paste0(figures.out.folder,"/p_bePairFactor_EnvDist_MnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.b.epair.factor),], aes(x=p.b.epair.factor, y=Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.sample.method)
dev.off()#___
## for p.sigma.spair & p.b.spair.factor with geo.dens method
png(paste0(figures.out.folder,"/p_sPairSigma_GeoDens_rndMnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.sigma.spair),], aes(x=p.sigma.spair, y=rnd.Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.b.spair.factor)
dev.off()#___
png(paste0(figures.out.folder,"/p_sPairSigma_GeoDens_rndMneRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.sigma.spair),], aes(x=p.sigma.spair, y=rnd.Mn.eRMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.b.spair.factor)
dev.off()#___
png(paste0(figures.out.folder,"/p_sPairSigma_GeoDens_MnRMSE.png"),height=500,width=500)
ggplot(parameter.results.tbl[!is.na(parameter.results.tbl$p.sigma.spair),], aes(x=p.sigma.spair, y=Mn.RMSE, fill=p.n.predictors.min)) +
  geom_boxplot() +
  facet_grid(. ~ p.b.spair.factor)
dev.off()#___

## End of plotting



#############################################################
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
