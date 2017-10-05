
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
## Make R package for gdm workflow
##  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## package libs
library(devtools)
library(roxygen2)

## root
pkg_root = '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/gdmEngine'

## DESCRIPTION file
DESCRIPTION = c('Package: gdmEngine',
                'Version: 0.01',
                paste('Date:', Sys.Date()),
                'Title: Workflow for GDM',
                'Description: Functions used to develop GDMs',
                'Author: Macroecological Modelling Team (CSIRO)',
                'Maintainer: Chris Ware <chris.ware@csiro.au>',
                'SystemRequirements: git with shell distribution'
                )

sink(paste(pkg_root, 'DESCRIPTION', sep = '/'))
cat(DESCRIPTION, sep = '\n')
sink()

## Build with devtools
setwd(pkg_root)
document()
build()
install()
## check() update examples before this is run

use_vignette('gdmEngine') # for later really...
## convert to html
library(knitr)
setwd(paste(getwd(), 'vignettes', sep = '/'))
rmarkdown::render(list.files())

## test if it's working
gitr.push(files = 'all')
?download_occurrences
?gitr.push
#rm(list = ls())
#detach(package:gdmEngine, unload = TRUE)

