#' Bayesian functional ANOVA
#'
#' This function carries out Bayesian functional ANOVA using the 
#' Normal Inverse Gamma Markov Grove method of Ma and Soriano (2016).
#'
#' @param W An object of class \code{DWT}.
#' @param TODO
#' @param is.kappa.fixed If \code{TRUE}, gamma_kappa and eta_kappa are fixed. 
#' If \code{FALSE} gamma_kappa and eta_kappa are defined using Empirical Bayes.
#' @param n_samp Number of posterior draws.
#' @param transition_mode Type of transition. The two options are \code{Markov} 
#' or \code{Independent}.
#' @param method Method used for find maxmimum of marginal likelihood. 
#' 
#' @return An object of class \code{grove}.
#' @export
#' @examples
#' a <- 'TO DO'

fANOVA <- function(W, 
                   formula, 
                   X, 
                   nu = 5, 
                   is.kappa.fixed = FALSE,
                   gamma.kappa = 0.3,
                   eta.kappa = 0.1,
                   n_samp = 500, 
                   transition_mode = "Markov", 
                   method = "Nelder-Mead") {
  
  if (is.kappa.fixed) {
    output <- .groveEB.fixed.kappa(W = W, 
                                   formula = formula, 
                                   X = X, 
                                   alpha = .InitAlphaPar(),
                                   nu = nu, 
                                   eta.rho = .InitEtaRhoPar(),
                                   eta.kappa = eta.kappa,
                                   gamma.rho = .InitGammaRhoPar(),
                                   gamma.kappa = gamma.kappa,
                                   n_samp = n_samp, 
                                   verbose = FALSE,
                                   transition_mode = transition_mode,
                                   method = method)
  } else {
    output <- .groveEB.all.random(W = W, 
                                  formula = formula, 
                                  X = X, 
                                  alpha = .InitAlphaPar(), 
                                  nu = nu, 
                                  eta.rho = .InitEtaRhoPar(),
                                  eta.kappa = eta.kappa,
                                  gamma.rho = .InitGammaRhoPar(),
                                  gamma.kappa = gamma.kappa,
                                  n_samp = n_samp, 
                                  verbose = FALSE,
                                  transition_mode = transition_mode,
                                  method = method)
  }
  return(output)
}