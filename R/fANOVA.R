#' Bayesian functional ANOVA
#'
#' This function carries out Bayesian functional ANOVA using the Normal Inverse Gamma Markov Grove method of 
#' Ma and Soriano (2016).
#'
#' @param W An object of class \code{DWT}.
#' @param alpha Hyperparameter.
#' @param beta Hyperparameter.
#' @param eta Hyperparameter.
#' @param gamma Hyperparameter.
#' @param is.alpha.fixed If \code{TRUE}, alpha is fixed. If \code{FALSE} alpha is defined using Empirical Bayes 
#' starting from value at \code{alpha}.
#' @param is.eta.fixed If \code{TRUE}, eta is fixed. If \code{FALSE} eta is defined using Empirical Bayes.
#' @param is.gamma.fixed If \code{TRUE}, gamma is fixed. If \code{FALSE} gamma is defined using Empirical Bayes.
#' @param nu Hyperparameter.
#' @param n_samp Number of posterior draws.
#' @param transition_mode Type of transition. The two options are \code{Markov} or \code{Independent}.
#' @param method Method used for find maxmimum of marginal likelihood. 
#' 
#' @return An object of class \code{grove}.
#' @export
#' @examples
#' a <- 'TO DO'

fANOVA <- function(W, 
                   formula, 
                   X, 
                   alpha = 0.5, 
                   beta = 1, 
                   nu = 5, 
                   n_samp = 500, 
                   transition_mode = "Markov", 
                   method ="Nelder-Mead") {
  
  if (is.fixed.eta & !is.fixed.gamma & !is.fixed.alpha) {
    output <- .groveEB.fixed.eta(W, 
                                 formula, 
                                 X, 
                                 alpha,
                                 beta,
                                 nu, 
                                 eta, 
                                 n_samp, 
                                 verbose = FALSE, 
                                 transition_mode , 
                                 method)
    
  } else if (is.fixed.eta & is.fixed.gamma & !is.fixed.alpha){
    output <- .groveEB.fixed.eta.and.gamma(W, 
                                           formula, 
                                           X, 
                                           alpha,
                                           beta, 
                                           nu, 
                                           eta, 
                                           gamma, 
                                           n_samp, 
                                           verbose = FALSE, 
                                           transition_mode, 
                                           method)
    
  } else if (!is.fixed.eta & !is.fixed.gamma & is.fixed.alpha) {
    output <- .groveEB.fixed.alpha(W, 
                                    formula, 
                                    X, 
                                    alpha, 
                                    beta, 
                                    nu, 
                                    n_samp, 
                                    verbose = FALSE, 
                                    transition_mode, 
                                    method)
    
  } else if (!is.fixed.eta & !is.fixed.gamma & !is.fixed.alpha) {
    output <- .groveEB(W, 
                        formula, 
                        X, 
                        alpha, 
                        beta, 
                        nu, 
                        n_samp, 
                        verbose = FALSE, 
                        transition_mode, 
                        method)
  } else {
    str.1 <- "ERROR: The selected combination of is.fixed.eta,"
    str.2 <- "is.fixed.gamma and is.fixed.alpha is not allowed."
    print(paste(str.1, str.2))
    return(0)
  }
  
  return(output)
  
}