
#' @title Push files to git
#'
#' @description Wrappers around git functions. Defaults set up for gdm project
#'
#' @param repo File path to any repository. Default is set to the gdm_workflow repo. Any other path will currently error.
#' @param files Files to be tracked. Default: this file (GDM_Workflow_functions.R). 
#' Paths to a file should be given, and should be relative to the repository root. Also accepts an argument 'all'.
#' @param msg Commit message (string) -m flag is set for this git function, so an argument is required. 
#' Defaults are set provided if left blank.
#' @param verbose Print standard output from git. Default TRUE
#'
#' @return std.output from git
#' 
#' @note make general...?
#' 
#' @examples gitr.push()
#' @examples gitr.push(files = 'all')
#' 
#' @export
gitr.push = function(repo = 'gdm_workflow', files = 'GDM_Workflow_Functions.R',
                     msg = 'updating workflow functions', verbose = TRUE){
  
  if(!file.exists('C:\\Program Files\\Git\\bin\\sh.exe')) 
    stop('Cannot find git - install it!')
  
  ## exec files should be on path of package
  sh_path = paste(.libPaths(), 'gdmEngine/exec/gitr.sh', sep = '/')
  bat_path = paste(.libPaths(), 'gdmEngine/exec/gitr.bat', sep = '/')
  assert_that(file.exists(sh_path))
  assert_that(file.exists(bat_path))
  ## account for any whitespace
  bat_path = paste('"', bat_path, '"', sep = '')
  sh_path = paste('"', sh_path, '"', sep = '')
  
  ## check repo
  if(repo == 'gdm_workflow') {
    rep_root = 
      '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow'
  } else {
    stop('pushing to repos other than gdm_workflow is on the to do list...')
  }
  
  ## check files to be tracked
  assert_that(length(files) == 1)
  if(files == 'all') {
    git_add = '--all'
  } else if(files != 'GDM_Workflow_Functions.R') {
    test_path = paste(rep_root, files, sep = '/')
    if(!file.exists(test_path)) 
      stop(paste('files arg is not a path relative to', rep_root))
    git_add = files
  } else {
    git_add = 'GDM_Workflow_Functions.R'
  }
  
  ## msg arg must be present
  assert_that(length(msg) == 1)
  if(git_add != 'GDM_Workflow_Functions.R' &  
     msg == 'updating workflow functions' & 
     git_add == '--all') {
    msg = 'updating all files'
  }
  ## current workaround for messages and white space which doesn't seem to 
  ## flow from batch to bash...
  msg = gsub("[[:space:]]", '_', msg)
  
  ## args
  send_args = paste(
    bat_path,
    sh_path,
    rep_root, 
    git_add, 
    msg, 
    sep = ' ')
  
  ## run
  if (verbose) {
    system(send_args, intern = TRUE)
  } else {
    system(send_args, intern = FALSE)
  } 
  
}






