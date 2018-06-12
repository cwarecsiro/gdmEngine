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
#library(data.table)
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
#soil.files <- list.files(path = "//lw-osm-02-cdc/OSM_CBR_LW_HCAS_work/HCAS2.0/HCAS2.0a/ENV/SOIL/TOP", full.names=TRUE, pattern = ".flt")
soil.files <- list.files(path = "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/env/SOIL/TOP", full.names=TRUE, pattern = ".flt")
env.files <- c(climate.files, terrain.files, soil.files)
env.files <- env.files[(substr(env.files, nchar(env.files)-3, nchar(env.files)) == ".flt")] # to remove some arcmap filenames
env.files <- env.files[-c(3,11,12,26,29,30,31,32,33,34,37,38,39,40)] # remove grids we don't want to assess in the modelling
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
n.pairs.model <- 144000 # equates to each site used 10 times
train.proportion <- 0.8
n.pairs.test <- 36000   # equates to each site used 10 times


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
n.pairs.model <- 67000 # equates to each site used 10 times
train.proportion <- 0.8
n.pairs.test <- 17000  # equates to each site used 10 times

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
n.pairs.model <- 14000 # equates to each site used 10 times
train.proportion <- 0.8
n.pairs.test <- 3500   # equates to each site used 10 times


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
n.pairs.model <- 60500  # equates to each site used 10 times
train.proportion <- 0.8
n.pairs.test <- 15000   # equates to each site used 10 times


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
Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/site_env_data_2018-06-04.csv")
#VASCULAR PLANTS -------
#Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/selected_gridcell_composition_2018-03-07.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/site_env_data_2018-05-31.csv")
#LAND SNAILS -------
Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/land_snails/selected_gridcell_composition_2018-03-09.csv")
Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/land_snails/site_env_data_2018-06-04.csv")
#Reptiles ---------
Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles/selected_gridcell_composition_2018-03-15.csv")
Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles/site_env_data_2018-06-04.csv")
##ENDTEMP##

## DERIVE A GDM ------------------------------------------------------------------------------------##
ptm <- proc.time()
GDM.Selection <- gdm_builder(site.env.data = Site.Env.Data, 
                             composition.data = Selected.records,
                             geo=FALSE,
                             n.pairs.train = n.pairs.model,
                             n.pairs.test = n.pairs.test,
                             correlation.threshold = 0.8,
                             selection.metric = 'D2',
                             sample.method = 'geowt',
                             Indiv.Dev.Explained.Min = 1.0,
                             n.predictors.min = 5,
                             domain.mask=Aus.domain.mask,
                             pcs.projargs="+init=epsg:3577",
                             bandwidth.geowt=150000,
                             bandwidth.skip=2,
                             bandwidth.DistFact=1,
                             geowt.RndProp=0.05,
                             output.folder = data.processing.folder,       
                             output.name = "gdm_mod_builder_results_GeowtSamp_NoGeo_V2") 
proc.time() - ptm

## SELECT A SET OF PREDICTORS FOR A GDM & APPLY SIGNIFICANCE TEST -----------------------------------##
final.mod.preds <- GDM.Selection$Mean.Final.GDM$predictors
geo.in = TRUE
if(final.mod.preds[1] == 'Geographic')
  {
  final.mod.preds <- final.mod.preds[-1]
  geo.in = TRUE
  }# end if final.mod.preds[1] == 'Geographic'
# or specify directly, for example: 
# final.mod.preds <- c('WDA','TXM','PTX','ELVR1000','SNDT','ECET','TNI','PTOT') #PLANTS
 final.mod.preds <- c('TXM','EAAS','TRI','ECET','ELVR1000','SNDT') #AMPHIBIANS
 final.mod.preds <- c('WDA','EAAS','TRA','EPI','BDWT','ELVR1000','CLYT') #geo.in=TRUE #LANDSNAILS
 
## ASSUMING YOU'RE HAPPY WITH A SET OF PREDICTORS, FIT A FINAL MODEL, INCLUDING CROSS-VALIDATION 
## ASSESSMENT AND SIGNIFICANCE TEST -----------------------------------------------------------------##
final.model.test <- gdm_build_single_model(site.env.data = Site.Env.Data, 
                                            composition.data = Selected.records,
                                            predictor.names = final.mod.preds,
                                            geo=geo.in,
                                            n.pairs.train = n.pairs.model,
                                            n.pairs.test = n.pairs.test,
                                            sample.method = 'geowt',
                                            b.used.factor=2,
                                            domain.mask=Aus.domain.mask,
                                            pcs.projargs="+init=epsg:3577",
                                            bandwidth.geowt=150000,
                                            bandwidth.skip=2,
                                            bandwidth.DistFact=1,
                                            geowt.RndProp=0.05,
                                            output.folder = file.path(data.processing.folder,"Final_GDM"),
                                            output.name = "gdm_build_FinMod_land_snails")
 
# ## SIGNIFICANCE TEST
# ptm <- proc.time()
# gdm_ext.sigtest(dllpath="//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/GDM_EXT_For_Karel/GDM4Rext.dll",
#                 wdpath = data.processing.folder,
#                 datatable = "//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/gdm_builder_FinMod_RandSamp_GDMtable_2018-05-10.csv", # GDM input table saved to .csv
#                 outname = "FinMod_RandSamp_sig_test",
#                 iterations = 100,
#                 do_geo = TRUE)
# proc.time() - ptm

## NOW TRANSFORM THE GRIDS BASED ON THE SELECTED MODEL ----------------------------------------------##
final.gdm <- final.model.test$Mean.Final.GDM
TransformGrids(gdm.model=final.gdm,
               env.grids.stk=env.stk,
               extrap.method = 'Conservative',
               output.folder = file.path(data.processing.folder,"Final_GDM")) 


## And plot the transformed grids to see if there are any spatial issues with the model projections
trans.grids <- list.files(path = file.path(data.processing.folder,"Final_GDM"), full.names=TRUE, pattern = ".flt")
for(i in 1:length(trans.grids))
  {
  next.ras <- raster(trans.grids[i])
  png(paste0(file.path(data.processing.folder,"Final_GDM"),"/plot_",names(next.ras),".png"),height=1000,width=1000)
  plot(next.ras)
  dev.off()#___
  } # end for i


## Create a plot of the main axes of compositional variation from the GDM transformed layers
trans.grids <- list.files(path = file.path(data.processing.folder,"Final_GDM"), full.names=TRUE, pattern = ".flt")
trans.grid.stk <- stack(trans.grids)
rastDat <- sampleRandom(trans.grid.stk, 200000)
pcaSamp <- prcomp(rastDat)
pcaRast <- predict(trans.grid.stk, pcaSamp, index=1:3)
# scale rasters
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255
#plotRGB(pcaRast, r=1, g=2, b=3)
# KM - rescale pca axes to represent their contribution
PC1.rng <- max(pcaSamp$x[,1]) - min(pcaSamp$x[,1])
PC2.rng <- max(pcaSamp$x[,2]) - min(pcaSamp$x[,2])
PC3.rng <- max(pcaSamp$x[,3]) - min(pcaSamp$x[,3])
PC2.scl <- PC2.rng/PC1.rng
PC3.scl <- PC3.rng/PC1.rng
pcaRast[[2]] <- 255 - (pcaRast[[2]] * PC2.scl)
pcaRast[[3]] <- 255 - (pcaRast[[3]] * PC3.scl)
plotRGB(pcaRast, r=1, g=2, b=3)



## ADITIONAL STUFF ##-----------------------???----------------------------------------------##
## Assess sitepair samples --------------------------------------------------------------------------##
SitePairOut <- sitepair_sample_assessor(site.env.data = Site.Env.Data, 
                                        composition.data = Selected.records ,
                                        n.pairs.train = n.pairs.model,
                                        sample.method = 'random',
                                        domain.mask=Aus.domain.mask,
                                        pcs.projargs="+init=epsg:3577",
                                        output.folder = data.processing.folder,       
                                        output.name = "sitepair_assess_amph") 

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