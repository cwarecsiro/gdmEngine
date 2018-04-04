#####################################################################################################
## 
## Convert grids to Australian Albers
##
#####################################################################################################

library(raster)
library(foreach)
library(doParallel)
library(sp)

# specify the files
climate.files <- list.files(path = "//lw-osm-02-cdc/OSM_CBR_LW_R51141_GPAA_work/ENV/A/OUT/1990", full.names=TRUE, pattern = ".flt")
terrain.files <- list.files(path = "//lw-osm-02-cdc/OSM_CBR_LW_R51141_GPAA_work/ENV/A/OUT/LAND", full.names=TRUE, pattern = ".flt")
soil.files <- list.files(path = "//osm-23-cdc/OSM_CBR_LW_DEE_work/source/env/SOIL/TOP", full.names=TRUE, pattern = ".flt")
env.files <- c(climate.files, terrain.files, soil.files)
env.files <- env.files[(substr(env.files, nchar(env.files)-3, nchar(env.files)) == ".flt")] # to remove some arcmap filenames
env.files <- env.files[-c(26,29,30,32,33,34,37,38,39,40)] # remove grids we don't want to assess in the modelling

# specify the output folder
out.path<-"//osm-23-cdc/OSM_CBR_LW_DEE_work/source/env/ENV_Albers"

# loop through the files, reproject into albers, save the reprojected grids
#for(i.ras in 1:length(env.files)){ #for serial processing
cl<-makeCluster(10) #setup parallel backend to use 12 processors
registerDoParallel(cl) 
# run the parallel loop over parameter combinations
foreach(i.ras=32:length(env.files), .packages='raster') %dopar% {
  this.ras<-raster(env.files[i.ras])
  projectRaster(this.ras,
                res=250,
                crs=CRS("+init=epsg:3577"),
                filename=paste0(out.path,"/",this.ras@data@names,".flt"))
  } # end for i.ras
stopCluster(cl)

# Just project the extent of a GCS raster to Albers
Aus.Albers.ext <- projectExtent(domain.mask,
                                crs=CRS("+init=epsg:3577"))

# Project points in GCS to Albers
site.pts.gcs <- SpatialPoints(coords=site.env.data[,c(3,4)],
                              proj4string=domain.mask@crs)

site.pts.pcs <- spTransform(site.pts.gcs,
                            CRSobj=CRS("+init=epsg:3577"))

# & test back-conversion
site.pts.pcs.gcs <- spTransform(site.pts.pcs,
                                CRSobj=domain.mask@crs)

