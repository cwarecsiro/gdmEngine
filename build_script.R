
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
                  'SystemRequirements: git with shell distribution',
                  'Licence: errr',
                  #paste('Authors@R:', unname(Sys.info()['user']))
                  'Imports: Rcpp (>= 0.11.4)',
                  'LinkingTo: Rcpp'
  )
  sink(paste(pkg_root, 'DESCRIPTION', sep = '/'))
  cat(DESCRIPTION, sep = '\n')
  sink()
  
  ## Build with devtools
  setwd(pkg_root)
  document()
  build()
  install(quick = TRUE)
}

#file.exists('\\\\ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/gdmEngine_0.01.tar.gz')
#install(pkg = '.', args = '-no-multiarch')
#system('R CMD INSTALL --library=Z:/users/bitbucket/gdm_workflow/gdmEngine gdmEngine_0.01.tar.gz ')

    
#get_latest = function(){
#  pkg_root = '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow'
#  setwd(pkg_root)
#  file.copy(paste0('gdmEngine'), .libPaths(), recursive = TRUE)
#  cat(paste0('gdmEngine copied to ', .libPaths()))
#  require(gdmEngine)
#}


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


