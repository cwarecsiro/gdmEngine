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


# Read in a spatial raster specifying the domain and resolution to be modelled
domain.mask <- raster("S:\\AUS0025\\CLIM\\MASK\\MASK0.flt")


## DOWNLOAD BIOLOGICAL DATA --------------------------------------------------------------##
# Download the biological data for a specified taxonomic group from the ALA, over the 
# specified spatial domain. Records will be available for download through a link emailed
# to the specified address.
ala_config(caching = "off")
occurrences(taxon = "class:Equisetopsida", # specify the taxon here
            wkt = paste0("POLYGON((",domain.mask@extent@xmin," ",domain.mask@extent@ymin,",",
                                    domain.mask@extent@xmin," ",domain.mask@extent@ymax,",",
                                    domain.mask@extent@xmax," ",domain.mask@extent@ymax,",",
                                    domain.mask@extent@xmax," ",domain.mask@extent@ymin,",",
                                    domain.mask@extent@xmin," ",domain.mask@extent@ymin,"))"),
            fq = "class:Equisetopsida", # specify the taxon here
            qa = "all", 
            method = "offline", 
            email = "mok010@csiro.au", # specify your email address here
            download_reason_id=7) 

## FILTER & FORMAT THE BIOLOGICAL DATA ---------------------------------------------------##
# For the specified data downloaded from ALA, prepare the data into a suitable format, 
# filtering out records that don't meet specified criteria.



## AGGREGATE THE BIOLOGICAL DATA TO GRID CELLS -------------------------------------------##
# Using specified criteria, aggregate the biological data to cells on the spacial grid, 
# removing data for cells where there are insuffucient records to be deemed a 'community sample'.









occurrences(taxon = focal.taxon,
            wkt ="POLYGON((112.9 -43.7425,112.9 -8,154 -8,154 -43.7425,112.9 -43.7425))",
            fq = focal.taxon,
            qa = "all", 
            method = "offline", 
            email = email.address, 
            download_reason_id=7) 
