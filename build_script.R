
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
## Make R package for gdm workflow
##  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## This first part of the code can just be hightlighted and run as a block
## starting here -->
## package libs


library(devtools)
library(roxygen2)
library(Rcpp)


update_build = function(){
  pkg_root = '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/gdmEngine'
  
  ## write DESCRIPTION file
  DESCRIPTION = c('Package: gdmEngine',
                  'Version: 0.01',
                  paste('Date:', Sys.Date()),
                  'Title: Workflow for GDM',
                  'Description: Functions used to develop GDMs',
                  paste('Author:', unname(Sys.info()['user'])),
                  'Maintainer: Chris Ware <chris.ware@csiro.au>',
                  'SystemRequirements: git with shell distribution'
                  #paste('Authors@R:', unname(Sys.info()['user']))
  )
  sink(paste(pkg_root, 'DESCRIPTION', sep = '/'))
  cat(DESCRIPTION, sep = '\n')
  sink()
  
  ## Build with devtools
  setwd(pkg_root)
  document()
  build()
  install()
}

get_latest = function(){
  pkg_root = '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow'
  setwd(pkg_root)
  file.copy(paste0('gdmEngine'), .libPaths(), recursive = TRUE)
  cat(paste0('gdmEngine copied to ', .libPaths()))
  require(gdmEngine)
}


## root
#pkg_root = '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/gdmEngine'

## write DESCRIPTION file
#DESCRIPTION = c('Package: gdmEngine',
#                'Version: 0.01',
#                paste('Date:', Sys.Date()),
#                'Title: Workflow for GDM',
#                'Description: Functions used to develop GDMs',
#                paste('Author:', unname(Sys.info()['user'])),
#                'Maintainer: Chris Ware <chris.ware@csiro.au>',
#                'SystemRequirements: git with shell distribution'
#                #paste('Authors@R:', unname(Sys.info()['user']))
#                )
#sink(paste(pkg_root, 'DESCRIPTION', sep = '/'))
#cat(DESCRIPTION, sep = '\n')
#sink()

## Build with devtools
#setwd(pkg_root)
#document()
#build()
#install()
## check() update examples before this is run

## <-- ending here


