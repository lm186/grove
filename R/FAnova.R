#' @title Bayesian functional ANOVA
#'
#' @description This function carries out Bayesian functional ANOVA using the 
#' Normal Inverse Gamma Markov Grove method of Ma and Soriano (2016).
#'
#' @param W An object of class \code{DWT}.
#' @param X Design matrix.
#' @param formula An object of class formula.
#' @param nu Hyperparameter.
#' @param is.kappa.fixed If \code{TRUE}, gamma.kappa and eta.kappa are fixed. 
#' If \code{FALSE} gamma_kappa and eta_kappa are defined using Empirical Bayes.
#' @param gamma.kappa Hyperparameter.
#' @param eta.kappa Hyperparameter.
#' @param n.samples Number of posterior draws.
#' @param transition.mode Type of transition. The two options are \code{Markov} 
#' or \code{Independent}.
#' @param method Method used for find maxmimum of marginal likelihood. 
#' 
#' @return An object of class \code{grove}.
#' @export
#' @examples
#' data <- GenerateSyntheticAnova(st.dev = 5, n.replicates = 5)
#' W <- DWT(data$noisy.Y)
#' X <- data$X
#' ans <- FAnova(W, X, ~ 1 + factorA + factorB)
#' denoised.data <- InvDWT(ans, x = c(0, 0, 1, 0))
#' PlotFun(denoised.data)

FAnova <- function(W, 
                   X, 
                   formula,                    
                   nu = 5, 
                   is.kappa.fixed = FALSE,
                   gamma.kappa = 0.3,
                   eta.kappa = 0.1,
                   n.samples = 500, 
                   transition.mode = "Markov", 
                   method = "Nelder-Mead") {
  
  if (is.kappa.fixed) {
    output <- .groveEB.fixed.kappa(W = W, 
                                   formula = as.formula(formula), 
                                   X = X, 
                                   alpha = .InitAlphaPar(),
                                   nu = nu, 
                                   eta.rho = .InitEtaRhoPar(),
                                   eta.kappa = eta.kappa,
                                   gamma.rho = .InitGammaRhoPar(),
                                   gamma.kappa = gamma.kappa,
                                   n.samples = n.samples, 
                                   verbose = FALSE,
                                   transition.mode = transition.mode,
                                   method = method)
  } else {
    output <- .groveEB.all.random(W = W, 
                                  formula = as.formula(formula), 
                                  X = X, 
                                  alpha = .InitAlphaPar(), 
                                  nu = nu, 
                                  eta.rho = .InitEtaRhoPar(),
                                  eta.kappa = eta.kappa,
                                  gamma.rho = .InitGammaRhoPar(),
                                  gamma.kappa = gamma.kappa,
                                  n.samples = n.samples, 
                                  verbose = FALSE,
                                  transition.mode = transition.mode,
                                  method = method)
  }
  return(output)
}