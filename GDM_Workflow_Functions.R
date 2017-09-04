###########################################################################################
##
## GDM Workflow functions - to support implementation of the GDM workflow
##
## Karel Mokany et al, (CSIRO) - Draft 30 August 2017
##
###########################################################################################


# Install libraries


# Load libraries
#library(ALA4R)
library(raster)

##______________________________________________________________________________________##
# A function to read in, filter and format data downloaded from the ALA
#Clean_ALA_data <- function(
ALA.download.folder <- "M:\\users\\mok010\\DoEE_GDM_Capacity_Building\\source\\Biological\\genus_Wandella\\data"
ALA.download.folder <- "M:\\users\\mok010\\DoEE_GDM_Capacity_Building\\source\\Biological\\class_Equisetopsida\\data"

#)
  {
  
  ## Standardise ALA input ######################
  # Get the name/s of the data file/s in this folder
  folder.files <- list.files(path = ALA.download.folder)
  data.files <- folder.files[substr(folder.files,1,4) == "data"]   
  # If there's only one data file, read it in
  if(length(data.files) == 1)
    {
    ALA.data <- read.csv(paste0(ALA.download.folder,"\\",data.files[1]), na.strings=c("na","NA","", " "),header=TRUE,sep="\t")
    } # end if length(data.files) == 1
  # if there are multiple data files, read them in & join them together
  if(length(data.files) > 1)
    {
    ALA.data <- read.csv(paste0(ALA.download.folder,"\\",data.files[1]), na.strings=c("na","NA","", " "),header=TRUE,sep="\t")
    for(i in 2:length(data.files))
      {
      next.read <- read.csv(paste0(ALA.download.folder,"\\",data.files[i]), na.strings=c("na","NA","", " "),header=TRUE,sep="\t")
      ALA.data <- rbind(ALA.data,next.read)
      } # end for i 
    } # end if length(data.files) > 1
  
  
  
  ## QA assertions
  
  ## Filter by spatial domain
  
  ## Filter by taxonomic resolution
  
  ## Filter by date range 
  
  ## Filter by spatial Uncertainty
  
  ## Filter by spatial coordinates (necessary data)
  
  ## Create cleaned data file for modelling + log file of how it was created
  
  
  
  
  
  
  
  
  
  } # end function Clean_ALA_data  


