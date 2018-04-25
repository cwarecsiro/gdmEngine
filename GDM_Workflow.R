###########################################################################################
##
## Generailsed Dissimilarity Modelling workflow
##
## Karel Mokany et al, (CSIRO) - Draft 30 August 2017
##
###########################################################################################

# Install libraries


# Load libraries
library(gdmEngine)
library(ALA4R)
library(raster)
library(data.table)
library(dplyr)
library(magrittr)
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
#soil.files <- list.files(path = "//lw-osm-02-cdc/OSM_CBR_LW_HCAS_work/HCAS2.0/HCAS2.0a/ENV/SOIL/TOP", full.names=TRUE, pattern = ".flt")
soil.files <- list.files(path = "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/env/SOIL/TOP", full.names=TRUE, pattern = ".flt")
env.files <- c(climate.files, terrain.files, soil.files)
env.files <- env.files[(substr(env.files, nchar(env.files)-3, nchar(env.files)) == ".flt")] # to remove some arcmap filenames
env.files <- env.files[-c(26,29,30,32,33,34,37,38,39,40)] # remove grids we don't want to assess in the modelling
env.stk <- stack(env.files, quick=TRUE) #env.stk <- stack(env.files)

# PLANTS INPUTS
species.names.file <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/vascular_plants/APC_and_Orchid_SpeciesNames.csv"
species.names <- read.csv(species.names.file)
species.names <- as.character(species.names[,1])
species.records.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/vascular_plants"
species.records.folder.raw <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/vascular_plants/raw_files"
data.processing.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants"
agg.cell.rad <- 2.25
min.rich.limit <- 10
max.rich.limit <- 400
min.rich.rad <- 200
min.rich.proportion <- 0.25
n.pairs.model <- 100000
train.proportion <- 0.8
n.pairs.test <- 20000


# AMPHIBIANS INPUTS
species.names.file <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/amphibians/AFD-20171211T130458.csv"
species.names <- read.csv(species.names.file)
species.names <- paste(species.names$GENUS, species.names$SPECIES)
species.names <- unique(species.names)
species.records.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/amphibians"
species.records.folder.raw <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/amphibians/raw_files"
data.processing.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians"
agg.cell.rad <- 2.25
min.rich.limit <- 2
max.rich.limit <- 50
min.rich.rad <- 200
min.rich.proportion <- 0.25
n.pairs.model <- 100000
train.proportion <- 0.8
n.pairs.test <- 20000

# LAND SNAIL INPUTS
species.names.file <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/land_snails/AusLandSnails_ALASpeciesList_9Mar18.csv"
species.names <- read.csv(species.names.file)
species.names <- species.names$Species.Name
species.names <- unique(species.names)
species.records.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/land_snails"
species.records.folder.raw <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/land_snails/raw_files"
data.processing.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/land_snails"
agg.cell.rad <- 2.25
min.rich.limit <- 2
max.rich.limit <- 50
min.rich.rad <- 50
min.rich.proportion <- 0.25
n.pairs.model <- 50000
train.proportion <- 0.8
n.pairs.test <- 10000


# REPTILE INPUTS
species.names.file <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/reptiles/AFD-20171211T113438.csv"
species.names <- read.csv(species.names.file)
species.names <- paste(species.names$GENUS, species.names$SPECIES)
species.names <- unique(species.names)
species.records.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/reptiles"
species.records.folder.raw <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/reptiles/raw_files"
data.processing.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles"
agg.cell.rad <- 2.25
min.rich.limit <- 3
max.rich.limit <- 50
min.rich.rad <- 200
min.rich.proportion <- 0.25
n.pairs.model <- 100000
train.proportion <- 0.8
n.pairs.test <- 20000


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
Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/selected_gridcell_composition_2018-03-05.csv")
Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/site_env_data_2018-03-05.csv")
#VASCULAR PLANTS -------
#Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/selected_gridcell_composition_2018-03-07.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/site_env_data_2018-03-07.csv")
#LAND SNAILS -------
#Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/land_snails/selected_gridcell_composition_2018-03-09.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/land_snails/site_env_data_2018-03-09.csv")
#Reptiles ---------
Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles/selected_gridcell_composition_2018-03-15.csv")
Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles/site_env_data_2018-03-15.csv")


## DERIVE A GDM --------------------------------------------------------------------------##
ptm <- proc.time()
Final.GDM <- gdm_builder(site.env.data = Site.Env.Data, 
                        composition.data = Selected.records ,
                        geo=TRUE,
                        n.pairs.train = n.pairs.model,
                        n.pairs.test = n.pairs.test,
                        selection.metric = 'D2',
                        sample.method = 'random',
                        Indiv.Dev.Explained.Min = 1.0,
                        n.predictors.min = 9,
                        domain.mask=Aus.domain.mask,
                        pcs.projargs="+init=epsg:3577",
                        output.folder = data.processing.folder,       
                        output.name = "gdm_builder_FinMod") 
proc.time() - ptm

## Assess sitepair samples --------------------------------------------------------------------------##
SitePairOut <- sitepair_sample_assessor(site.env.data = Site.Env.Data, 
                                        composition.data = Selected.records ,
                                        n.pairs.train = n.pairs.model,
                                        sample.method = 'random',
                                        domain.mask=Aus.domain.mask,
                                        pcs.projargs="+init=epsg:3577",
                                        output.folder = data.processing.folder,       
                                        output.name = "sitepair_assess_amph") 


## ADITIONAL STUFF ##-----------------------???----------------------------------------------##
#Random sample
Pairs.Table.Rnd <- sitepair_sample_random(site.env.data = Site.Env.Data,
                                          n.pairs.target = n.pairs.model)
#Geographic distance based sample
Pairs.Table.Geo <- sitepair_sample_geographic(site.env.data = Site.Env.Data,
                           n.pairs.target = n.pairs.model,
                           a.used=0.05, 
                           b.used.factor=2, 
                           c.used=3, 
                           a.dpair=0.05, 
                           b.dpair.factor=2, 
                           c.dpair=3)

# Environmental distance based sample
Pairs.Table.Env <- sitepair_sample_environment(site.env.data = Site.Env.Data,
                                               n.pairs.target = n.pairs.model,
                                               env.colnames = c('PTA','TXX'),
                                               b.used.factor=2, 
                                               c.used=3, 
                                               b.epair.factor=1, 
                                               c.epair=3)

# Neighbourhood site-density based sample (weighted towards less sampled areas)
Pairs.Table.Dens <- sitepair_sample_density(site.env.data = Site.Env.Data,
                                            n.pairs.target = n.pairs.model,
                                            domain.mask = Aus.domain.mask,
                                            a.used=0.05, 
                                            b.used.factor=2, 
                                            c.used=3, 
                                            sigma.spair=0.5,
                                            a.spair=0.05, 
                                            b.spair.factor=1.0, 
                                            c.spair=1)

