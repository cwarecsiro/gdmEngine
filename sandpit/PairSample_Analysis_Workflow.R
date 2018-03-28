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
#Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/selected_gridcell_composition_2018-03-05.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/site_env_data_2018-03-05.csv")
#VASCULAR PLANTS -------
#Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/selected_gridcell_composition_2018-03-07.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/vascular_plants/site_env_data_2018-03-07.csv")
#LAND SNAILS -------
#Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/land_snails/selected_gridcell_composition_2018-03-09.csv")
#Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/land_snails/site_env_data_2018-03-09.csv")
#Reptiles ---------
Selected.records <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles/selected_gridcell_composition_2018-03-15.csv")
Site.Env.Data <- read.csv("//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles/site_env_data_2018-03-15.csv")
##ENDTEMP##

## NOW RUN THE SITEPAIR SAMPLING PARAMETER ASSESSMENT ------------------------------------##
# Now create a parameter table listing the possible combinations 
parameter.tbl <- expand.grid(p.sample.method=c('random','geodist','envdist','geodens'),
                             p.selection.metric=c('RMSE'), #add eRMSE?
                             p.geo=c(TRUE,FALSE),
                             p.n.predictors.min=c(5,10),
                             p.Indiv.Dev.Explained.Min=c(0.5,1.0),
                             p.n.pairs.train=c(100000,200000),
                             p.n.pairs.test=c(20000,40000),
                             p.b.used.factor=c(1,2,3),
                             p.b.dpair.factor=c(1,2,3),	
                             p.b.epair.factor=c(1,2,3),	
                             p.sigma.spair=c(0.25,0.5,1.0),	
                             p.b.spair.factor=c(1,2,3),
                             stringsAsFactors=FALSE)
# Remove irrelevant parameter combinations
parameter.tbl$p.b.dpair.factor[(parameter.tbl$p.sample.method != 'geodist')]<-NA
parameter.tbl$p.b.epair.factor[(parameter.tbl$p.sample.method != 'envdist')]<-NA
parameter.tbl$p.b.spair.factor[(parameter.tbl$p.sample.method != 'geodens')]<-NA
parameter.tbl$p.sigma.spair[(parameter.tbl$p.sample.method != 'geodens')]<-NA
parameter.tbl$p.n.predictors.min[(parameter.tbl$p.Indiv.Dev.Explained.Min==1)]<-5
parameter.tbl$p.n.predictors.min[(parameter.tbl$p.Indiv.Dev.Explained.Min==0.5)]<-10
parameter.tbl$p.n.pairs.train[(parameter.tbl$p.n.pairs.test==40000)]<-200000
parameter.tbl$p.n.pairs.train[(parameter.tbl$p.n.pairs.test==20000)]<-100000
parameter.tbl<-unique(parameter.tbl)
parameter.tbl<-parameter.tbl[order(parameter.tbl$p.sample.method),]
parameter.tbl$run.name<-paste0("Amph_",parameter.tbl$p.sample.method,"_",c(1:nrow(parameter.tbl)))
# Now run all the parameter combinations
analysis.out.folder<-"//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles/SitePairSampleTesting"
write.csv(parameter.tbl,paste0(analysis.out.folder,"/reptile_parameters.csv"),row.names = FALSE)
# for(i.run in 1:nrow(parameter.tbl))
#   {
#   ptm <- proc.time()
#   This.GDM <- gdm_builder(site.env.data = Site.Env.Data,
#                            composition.data = Selected.records ,
#                            geo=parameter.tbl$p.geo[i.run],
#                            n.pairs.train = parameter.tbl$p.n.pairs.train[i.run],
#                            n.pairs.test = parameter.tbl$p.n.pairs.test[i.run],
#                            correlation.threshold = 0.75,
#                            selection.metric = parameter.tbl$p.selection.metric[i.run],
#                            sample.method=parameter.tbl$p.sample.method[i.run],
#                            Indiv.Dev.Explained.Min = parameter.tbl$p.Indiv.Dev.Explained.Min[i.run],
#                            n.predictors.min = parameter.tbl$p.n.predictors.min[i.run],
#                            b.used.factor=parameter.tbl$p.b.used.factor[i.run],
#                            b.dpair.factor=parameter.tbl$p.b.dpair.factor[i.run],
#                            b.epair.factor=parameter.tbl$p.b.epair.factor[i.run],
#                            sigma.spair=parameter.tbl$p.sigma.spair[i.run],
#                            spair.factor=parameter.tbl$p.b.spair.factor[i.run],
#                           domain.mask=Aus.domain.mask,
#                           output.folder = analysis.out.folder,
#                           output.name = parameter.tbl$run.name[i.run])
#   proc.time() - ptm
#   }# end for i.run
# Or for parallel implementation...
library(foreach)
library(doParallel)
cl<-makeCluster(12) #setup parallel backend to use 10 processors
registerDoParallel(cl) 
# run the parallel loop over parameter combinations
foreach(i.run=1:nrow(parameter.tbl), .packages='gdmEngine') %dopar% {
  This.GDM <- gdm_builder(site.env.data = Site.Env.Data, 
                          composition.data = Selected.records,
                          geo=parameter.tbl$p.geo[i.run],
                          n.pairs.train = parameter.tbl$p.n.pairs.train[i.run],
                          n.pairs.test = parameter.tbl$p.n.pairs.test[i.run],
                          correlation.threshold = 0.7,
                          selection.metric = parameter.tbl$p.selection.metric[i.run],
                          sample.method=parameter.tbl$p.sample.method[i.run],
                          Indiv.Dev.Explained.Min = parameter.tbl$p.Indiv.Dev.Explained.Min[i.run],
                          n.predictors.min = parameter.tbl$p.n.predictors.min[i.run],
                          b.used.factor=parameter.tbl$p.b.used.factor[i.run],
                          b.dpair.factor=parameter.tbl$p.b.dpair.factor[i.run],
                          b.epair.factor=parameter.tbl$p.b.epair.factor[i.run],
                          sigma.spair=parameter.tbl$p.sigma.spair[i.run],
                          spair.factor=parameter.tbl$p.b.spair.factor[i.run],
                          domain.mask=Aus.domain.mask,
                          output.folder = analysis.out.folder,       
                          output.name = parameter.tbl$run.name[i.run])  
} # end for i.run
stopCluster(cl)

## SUMMARISE THE RESULTS FROM THE SITEPAIR SAMPLING ASSESSMENT ---------------------------##
# Specify the folder where the output files are held
analysis.out.folder<-"//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/reptiles/SitePairSampleTesting"
analysis.out.files<-as.character(list.files(analysis.out.folder))
# Read in the parameter table
parameter.tbl <- read.csv(paste0(analysis.out.folder,"/reptile_parameters.csv"))
# Create a table to catch the summary data
results.tbl <- data.frame(ProcTime = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          nvars.FinalMod = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          n.pairs.used = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          D2.FinalMod = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Mn.RMSE = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Mn.eRMSE = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Mn.MAE = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          rnd.Mn.RMSE = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          rnd.Mn.eRMSE = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          rnd.Mn.MAE = as.numeric(rep(NA,times=nrow(parameter.tbl))),
                          Mn.Dissimilarity = as.numeric(rep(NA,times=nrow(parameter.tbl))) )
# Now loop through parameter combos, see if there's an output file for each, read it in and
# summarise the results, catching in the output table.
for(i.run in 1:nrow(parameter.tbl))
  {
  #find if there is an output file to match the run.name
  paste(as.character(parameter.tbl$run.name[i.run]),"_",sep='')
  f.name <- analysis.out.files[pmatch(paste(as.character(parameter.tbl$run.name[i.run]),"_",sep=''), analysis.out.files)]
  if(!is.na(f.name))
    {
    # open the file
    load(file.path(analysis.out.folder,f.name))
    # catch the results 
    results.tbl$ProcTime[i.run] = GDM_Builder_results$ProcessingTime[3]
    results.tbl$nvars.FinalMod[i.run] = length(GDM_Builder_results$Predictors)
    results.tbl$n.pairs.used[i.run] = length(GDM_Builder_results$Mean.Final.GDM$observed)
    results.tbl$D2.FinalMod[i.run] = mean(GDM_Builder_results$Deviance.Explained)
    results.tbl$Mn.RMSE[i.run] = mean(GDM_Builder_results$Root.Mean.Squre.Error)
    results.tbl$Mn.eRMSE[i.run] = mean(GDM_Builder_results$Equalised.RMSE)
    results.tbl$Mn.MAE[i.run] = mean(GDM_Builder_results$Mean.Absolute.Error)
    results.tbl$rnd.Mn.RMSE[i.run] = mean(GDM_Builder_results$rnd.Root.Mean.Squre.Error)
    results.tbl$rnd.Mn.eRMSE[i.run] = mean(GDM_Builder_results$rnd.Equalised.RMSE)
    results.tbl$rnd.Mn.MAE[i.run] = mean(GDM_Builder_results$rnd.Mean.Absolute.Error)
    results.tbl$Mn.Dissimilarity[i.run] = mean(GDM_Builder_results$dissimilarity.summary[,4])
    # remove the results
    rm(GDM_Builder_results)
    }# end if !is.na(f.name)
  } # end for i.run
parameter.results.tbl <- cbind(parameter.tbl, results.tbl) 
write.csv(parameter.results.tbl, paste0(analysis.out.folder,"/reptile_parameters_results.csv"), row.names = FALSE)


## Plot results ###############################################
library(ggplot2)
parameter.results.tbl<-read.csv('//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/SitePairSampleTesting/amph_parameters_results.csv')
figures.out.folder<-"//osm-23-cdc/OSM_CBR_LW_DEE_work/processing/biol/amphibians/SitePairSampleTesting/Figures"
# convert factors to factors
parameter.results.tbl$p.geo<-as.factor(parameter.results.tbl$p.geo)
parameter.results.tbl$p.n.predictors.min<-as.factor(parameter.results.tbl$p.n.predictors.min)
parameter.results.tbl$p.n.pairs.train<-as.factor(parameter.results.tbl$p.n.pairs.train)
parameter.results.tbl$p.b.used.factor<-as.factor(parameter.results.tbl$p.b.used.factor)
parameter.results.tbl$p.b.dpair.factor<-as.factor(parameter.results.tbl$p.b.dpair.factor)
parameter.results.tbl$p.b.epair.factor<-as.factor(parameter.results.tbl$p.b.epair.factor)
parameter.results.tbl$p.sigma.spair<-as.factor(parameter.results.tbl$p.sigma.spair)
parameter.results.tbl$p.b.spair.factor<-as.factor(parameter.results.tbl$p.b.spair.factor)

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

