#'@title Negative exponential link function
#'
#'@description A negative exponential link function for GLM.
#'
#'@return GLM link.
#'
#'@examples vv <- negexp()
#'
#'@export
negexp <- function()
  {
	## link
    linkfun <- function(mu) -log(1-mu)
	## inverse link
    linkinv <- function(eta) 1-exp(-eta)
	## derivative of inverse link wrt eta
    mu.eta <- function(eta) exp(-eta)
    valideta <- function(eta) all(is.finite(eta))
    link <- paste0("negexp")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
  } # end negexp
