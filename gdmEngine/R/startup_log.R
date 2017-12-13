.onAttach <- function(libname, pkgname) {
  pkg_root = '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/gdmEngine'
  check_edits = readLines(paste(pkg_root, 'DESCRIPTION', sep = '/'))
  last_mod = gsub('Date: ', '', check_edits[grep('Date', check_edits)])
  last_ed = gsub('Author: ', '', check_edits[grep('Author', check_edits)])
  
  msg1 = paste0('Last build on ', last_mod)
  msg2 = paste0('Last modified by ', last_ed)
  msg3 = rep('-', nchar(msg2) + 1)
  msg4 = '\n'
  
  packageStartupMessage(msg4)
  packageStartupMessage(msg3)
  packageStartupMessage(msg1)
  packageStartupMessage(msg2)
  packageStartupMessage(msg3)
  packageStartupMessage(msg4)
  
}

