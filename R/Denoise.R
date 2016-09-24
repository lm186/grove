#' Bayesian wavelet denoising
#'
#' This function carries out Bayesian wavelet denoising using the Normal 
#' Inverse Gamma Markov Tree method of Ma and Soriano (2016).
#'
#' @param W An object of class \code{DWT}.
#' @param alpha Hyperparameter.
#' @param nu Hyperparameter controlling variance heterogeneity. If \code(Inf),
#' then the variance is identical for all nodes.
#' @param n_samp Number of posterior draws.
#' @param transition_mode Type of transition. 
#' The two options are \code{Markov} or \code{Independent}.
#' @param method Method used for find maxmimum of marginal likelihood. 
#' 
#' @return An object of class \code{grove}.
#' @export
#' @examples
#' a <- 'TO DO'

Denoise <- function(W, 
                    alpha = 0.5, 
                    nu = 5, 
                    n_samp = 500, 
                    transition_mode = "Markov",
                    method = "Nelder-Mead") {
  
  formula <- formula(~ 1)
  X <- model.matrix(frml, matrix(1, ncol = 1, nrow = nrow(W)))
   
  output <- .groveEB.all.random(W = W, 
                                formula = formula, 
                                X = X, 
                                alpha = alpha, 
                                nu = nu, 
                                eta.rho = .InitEtaRhoPar(),
                                eta.kappa = .InitEtaKappaPar(),
                                gamma.rho = .InitGammaRhoPar(),
                                gamma.kappa = .InitGammaKappaPar(),
                                n_samp = n_samp, 
                                verbose = FALSE,
                                transition_mode = transition_mode,
                                method = method)
  return(output)
}  