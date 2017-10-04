
###############################################################################
##
## Utility functions to support gdm workflow functions
##
###############################################################################


doc <- function(fn){
  'Summary'
  '-------'
  'Returns doc string' 
  ''
  'Args'
  '----'
  'fn: function'
  ''
  'Returns'
  '-------'
  ''
  'Notes'
  'library(docstring) has a more sophisticated approach'
  'This is a minimalist approach'
  'end doc string'
  # Example:
  #   square <- function(x){
  #     # whatever
  #     'Summary'
  #     '-------'
  #     'Squares a number'
  #     ''
  #     'Args'
  #     '----'
  #     'x: number'
  #     ''
  #     'Returns'
  #     '-------'
  #     'number'
  #     'end doc string'
  #     return(x ^ 2)
  #   }
  #   doc(square)
  dstring <- body(fn)
  endoc <- grep('end doc string', dstring)
  if(length(endoc) == 0) stop('Need an <end doc sting> statement')
  elements <- dstring[2:(endoc-1)]
  elements <- paste(elements, rep('\n', length(elements)))
  indents <- c(grep('--', elements), grep('--', elements)-1) 
  #indents <- outersect(seq(1:length(elements)), indents)
  indents <- sort(c(setdiff(seq(1:length(elements)), indents), 
                    setdiff(indents, seq(1:length(elements)))))
  elements[indents] <- paste('\  ', elements[indents], sep = ' ')
  cat(elements, sep = '')
}


gitr.push = function(repo = 'gdm_workflow', files = 'GDM_Workflow_Functions.R',
                     msg = 'updating workflow functions', verbose = TRUE){
  
  'Summary'
  '-------'
  'Wrappers for git calls'
  ''
  'Args'
  '----'
  'repo: any repo. Default gdm_workflow. Any other will currently error'
  'files: files to add to be tracked. Default: this file. Also accepts all.'
  'msg: string. Message for commit.'
  'verbose = print output from git. Default TRUE'
  ''
  'Returns'
  '-------'
  'output from git'
  ''
  'Depends'
  '-------'
  'pkgs: assertthat'
  'apps: git'
  ''
  'Examples'
  '--------'
  'TODO: make general...?'
  'end doc string'
  
  if(!file.exists('C:\\Program Files\\Git\\bin\\sh.exe')) 
    stop('Cannot find git - install it!')
  
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
    "//ces-10-cdc//OSM_CDC_MMRG_work//users//bitbucket//gdm_workflow//gitr.bat", 
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
