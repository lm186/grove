#' Bayesian wavelet denoising
#'
#' This function carries out Bayesian wavelet denoising using the Normal Inverse Gamma Markov Tree method of 
#' Ma and Soriano (2016).
#'
#' @param W An object of class \code{DWT}.
#' @param alpha Hyperparameter.
#' @param beta Hyperparameter.
#' @param is.alpha.fixed If \code{TRUE}, alpha is fixed. If \code{FALSE} alpha is defined using Empirical Bayes 
#' starting from value at \code{alpha}.
#' @param is.beta.fixed If \code{TRUE}, beta is fixed. If \code{FALSE} beta is defined using Empirical Bayes 
#' starting from value at \code{beta}. If \code{is.alpha.fixed = TRUE}, 
#' then also \code{is.beta.fixed} must be \code{TRUE}.
#' @param nu Hyperparameter.
#' @param n_samp Number of posterior draws.
#' @param transition_mode Type of transition. The two options are \code{Markov} or \code{Independent}.
#' @param method Method used for find maxmimum of marginal likelihood. 
#' 
#' @return An object of class \code{grove}.
#' @export
#' @examples
#' a <- 'TO DO'

Denoise <- function(W, 
                    alpha = 0.5, 
                    beta = 1, 
                    is.alpha.fixed = TRUE,
                    is.beta.fixed = TRUE,
                    nu = 5, 
                    n_samp = 500, 
                    transition_mode = "Markov",
                    method = "Nelder-Mead") {
  
  formula <- formula(~1)
  X <- model.matrix(frml, matrix(1, ncol = 1, nrow = nrow(W)))
  
  if (is.alpha.fixed & is.beta.fixed) {
    output <- .groveEB.denoising.fixed.alpha.beta(W, 
                                                  formula, 
                                                  X, 
                                                  alpha, 
                                                  beta, 
                                                  nu, 
                                                  n_samp, 
                                                  FALSE,
                                                  transition_mode,
                                                  method) 
    
  } else if (!is.alpha.fixed & is.beta.fixed) {
    output <- .groveEB.denoising.fixed.beta(W, 
                                            formula, 
                                            X, 
                                            alpha, 
                                            beta, 
                                            nu, 
                                            n_samp, 
                                            FALSE,
                                            transition_mode,
                                            method) 
    
  } else if (!is.alpha.fixed & !is.beta.fixed) {
    output <- .groveEB.denoising(W, 
                                 formula, 
                                 X, 
                                 alpha, 
                                 beta, 
                                 nu, 
                                 n_samp, 
                                 FALSE,
                                 transition_mode,
                                 method) 
    
  } else {
    print("ERROR: It is not possible to fix beta but not alpha.")
    return(0)
  }
  
  return(output)
  
}