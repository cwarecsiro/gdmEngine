###########################################################################################
##
## Generailsed Dissimilarity Modelling workflow
##
## Karel Mokany et al, (CSIRO) - Draft 30 August 2017
##
###########################################################################################

# Install libraries


# Load libraries
library(ALA4R)
library(raster)
library(gdmEngine)
library(data.table)
#library(dplyr)
library(magrittr)
library(plyr)
library(assertthat)
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
species.names.file <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/vascular_plants/APC_and_Orchid_SpeciesNames.csv"
species.names <- read.csv(species.names.file)
species.names <- as.character(species.names[,1])
species.records.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/vascular_plants"
species.records.folder.raw <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/biol/vascular_plants/raw_files"
data.processing.folder <- "//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants"

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
Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/selected_gridcell_composition_2017-12-14.csv")
Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/site_env_data_2018-01-05.csv")

## DERIVE A GDM --------------------------------------------------------------------------##
ptm <- proc.time()
Final.GDM <- gdm_builder(site.env.data = Site.Env.Data, 
                        composition.data = Selected.records ,
                        geo=TRUE,
                        train.proportion = 0.8,
                        n.pairs.train = n.pairs.model,
                        n.pairs.test = n.pairs.test,
                        n.crossvalid.tests = 10,
                        correlation.threshold = 0.7,
                        selection.metric = 'RMSE',
                        Indiv.Dev.Explained.Min = 1.0,
                        n.predictors.min = 8,
                        output.folder = data.processing.folder,       
                        output.name = "gdm_builder_FinMod") 
proc.time() - ptm

## ADITIONAL STUFF ##---------------------------------------------------------------------##
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