
#' @title Configure git call 
#'
#' @description Configure parameters to push file to git, or save parameters for future. Defaults set up for gdm project
#'
#' @param ... List of named args to pass. 
#' @section Valid options:
#' \itemize{
#'  \item{"shell_exe"}
#'  \item{"repo_root"}
#'  \item{"repo"}
#'  \item{"files"}
#'  \item{"msg"}
#'  \item{"verbose"}
#' }
#'
#' @return list of git config options
#' 
#' @note make general...?
#' 
#' @examples gitr.config()
#' 
#' @export			 
gitr.config = function(...){
  
  ## get user options
  user_opts <- list(...)
  
  ## get info about user and check for a permanent file dst
  if(unname(Sys.info()['sysname']) == 'Windows'){
    dst = '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/utils/R'
    assert_that(dir.exists(dst))
  } else if (unname(Sys.info()['sysname']) == 'Linux'){
    ## what should this be?
    dst = paste0('...') 
    assert_that(dir.exists(dst))
  } else {
    stop('Cannot find a location to write config file to')
  }
  
  ## default options
  default_opts = list(
    shell_exe = 'C:\\Program Files\\Git\\bin\\sh.exe',
    repo_root =  '//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow',
    repo = 'gdm_workflow',
    files = 'GDM_Workflow_Functions.R',
    msg = 'updating workflow functions',
    verbose = TRUE
  )
  
  ## check if config file exists
  fn = paste0('config_', unname(Sys.info()['user']), '_gitr.yaml')
  src = list.files(dst, pattern = fn, full.names = TRUE)
  if(!length(src)) { 
    ## obj to be modified by user args if there are any
    opts = default_opts
    
    ## if file exists, read it  
  } else {
    src = paste0(dst, '/', fn)
    opts = yaml::yaml.load_file(src) 
  }
  
  if (length(user_opts)) {
    for (i in 1:length(user_opts)) {
      this_opt <- names(user_opts)[i]
      if (! (this_opt %in% names(default_opts))) {
        cat(paste0("'", this_opt, "'", ' is not a valid option. Should be one of:'), 
            names(opts), sep = '\n')
      } else {
        opts[[this_opt]] = paste(user_opts[i])
      }
    } # end user_opts 
  } # end check for user_opts
  
  ## if not default, write
  if(!identical(opts, default_opts)){
    sink(src)
    cat(as.yaml(opts))
    sink()
  }
  
  return(opts)
  
}

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
gitr.push = function(verbose = TRUE, ...){
  
  if("package:gdmEngine" %in% search()){
    ## exec files should be on path of package
    sh_path = paste(.libPaths(), 'gdmEngine/exec/gitr.sh', sep = '/')
    bat_path = paste(.libPaths(), 'gdmEngine/exec/gitr.bat', sep = '/')
    assert_that(file.exists(sh_path))
    assert_that(file.exists(bat_path))
  } else {
    ## find these paths
    sh_path = "//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/lib/gdmEngine/exec/gitr.sh"
    bat_path = "//ces-10-cdc/OSM_CDC_MMRG_work/users/bitbucket/gdm_workflow/lib/gdmEngine/exec/gitr.bat"
    assert_that(file.exists(sh_path))
    assert_that(file.exists(bat_path))
  }  
  ## account for any whitespace
  bat_path = paste('"', bat_path, '"', sep = '')
  sh_path = paste('"', sh_path, '"', sep = '')
  
  user_opts <- list(...)
  ## check for user config
  opts = gitr.config()
  ## check for sh.exe
  if(!file.exists(opts[['shell_exe']])) stop('Cannot find git - install it!')
  
  ## update opts with user options  
  if (length(user_opts)){
    opts = gitr.config()
    for (i in 1:length(user_opts)) {
      this_opt <- names(user_opts)[i]
      if (! (this_opt %in% names(opts))) {
        cat(paste0("'", this_opt, "'", ' is not a valid option. Should be one of:'), 
            names(opts), sep = '\n')
      } else {
        opts[[this_opt]] = paste(user_opts[i])
      }
    } # end user_opts 
  }  
  
  ## check files to be tracked
  assert_that(length(opts$files) == 1)
  if(opts$files == 'all') {
    opts$files = '--all'
  } else if(opts$files != 'GDM_Workflow_Functions.R') {
    test_path = paste(opts$repo_root, opts$files, sep = '/')
    if(!file.exists(test_path)) 
      stop(paste("'", opts$files, "'", ' is not a path relative to ', opts$repo_root))
  } 
  
  ## msg arg must not be empty
  assert_that(length(opts$msg) == 1)
  
  ## current workaround for messages and white space which doesn't seem to 
  ## flow from batch to bash...
  opts$msg = gsub("[[:space:]]", '_', paste(opts$msg))
  
  ## args
  send_args = paste(
    bat_path,
    sh_path,
    opts$repo_root, 
    opts$files, 
    opts$msg, 
    sep = ' ')
  
  ## run
  if (verbose) {
    system(send_args, intern = TRUE)
  } else {
    system(send_args, intern = FALSE)
  } 
  
}




