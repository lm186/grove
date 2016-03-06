.groveEB <- function(W, 
                    formula, 
                    X, 
                    alpha = 0.5, 
                    beta = 1, 
                    nu = 5, 
                    n_samp = 500, 
                    verbose = FALSE, 
                    transition_mode = "Markov", 
                    method ="Nelder-Mead")
{
  if(class(W)!="DWT")
  {
    print("ERROR: W should be a DWT object")
    return(0)
  }

  if(strsplit(as.character(formula)[2],split="")[[1]][1] != 1)
  {
    print("ERROR: formula should include an intercept.")
    return(0)
  }
  
  if (transition_mode == "Markov") transition_mode = 1
  else if (transition_mode == "Independent") transition_mode = 0
  else {
    transition_mode = 1
    print("WARNING: Unrecognized transition mode. Default to Markov.")
  }
  
  XX = model.matrix(formula, X)
  p = table(attr(XX,"assign"))  
  p_len = length(p)
  
  frml = paste("C", paste(formula, sep="", collapse=" "), sep = " ")
  C_hat = lm(frml, data = cbind(C=W$C,X))$coefficients
  
  init_state = rep(0,p_len)
  sigma_par = mad(tail(t(W$D),n=1))
  ## eta_par = c(2, rep(0.1,p_len-1))
  eta_par = c(2, rep(0.1,max(p_len-1,1)))
  gamma_par = 0.3
  tau_par = rep( (sd( head(t(W$D),n=1))/sigma_par)^2,p_len)
  nu_par = nu
  alpha_par = alpha
  beta_par = beta
 
  fr <- function(input) 
  { 
    tau = input[1:p_len]
    eta.rho = input[p_len+1]
    eta.kap = input[p_len+2]
    gamma.rho = input[p_len+3]
    gamma.kap = input[p_len+4]
    sigma = input[p_len+5]
    alpha = input[p_len+6]

    if(nu_par != Inf) {
      nu = input[p_len+7]
    } else nu = log(999999999)
    
    output = fitGroveML(W$D, 
                        XX, 
                        p = p, 
                        tau_par = exp(tau), 
                        eta_par = c(exp(eta.rho), rep(exp(eta.kap),p_len-1)), 
                        gamma_par = c(0.5*exp(gamma.rho)/(1+exp(gamma.rho)),
                                    rep(0.5*exp(gamma.kap)/(1+exp(gamma.kap)),p_len-1)),
                        init_state = init_state,
                        nu = exp(nu), 
                        sigma = exp(sigma), 
                        alpha = alpha, 
                        beta = beta, 
                        transition_mode = transition_mode)
    if(verbose == TRUE)
      print(-output$marginal_likelihood)
    return(-output$marginal_likelihood)
  }  
  
  if (nu != Inf) {
    
    empirical_bayes = optim(par=c(log(tau_par), 
                                  log(eta_par[1]), 
                                  log(eta_par[2]),
                                  log(gamma_par/(0.5-gamma_par)), 
                                  log(gamma_par/(0.5-gamma_par)),
                                  log(sigma_par),
                                  alpha_par,
                                  log(nu_par)), 
                            fn=fr,  
                            method = method)
  }
  
  else {
    empirical_bayes = optim(par=c(log(tau_par), 
                                  log(eta_par[1]), 
                                  log(eta_par[2]),
                                  log(gamma_par/(0.5-gamma_par)), 
                                  log(gamma_par/(0.5-gamma_par)),
                                  log(sigma_par),
                                  alpha_par), 
                            fn=fr,  
                            method = method)
  }
  
  
  
  tau_par = exp(empirical_bayes$par[1:p_len])
  eta_par = c( exp(empirical_bayes$par[p_len+1]), 
               rep(exp(empirical_bayes$par[p_len+2]),p_len-1))
  gamma_par = c( 0.5*exp(empirical_bayes$par[p_len+3])/(1+exp(empirical_bayes$par[p_len+3])),
                 rep(0.5*exp(empirical_bayes$par[p_len+4])/(1+exp(empirical_bayes$par[p_len+4])),p_len-1))
  sigma_par = exp(empirical_bayes$par[p_len+5])
  alpha_par = empirical_bayes$par[p_len+6]
  beta_par = beta
  
  if (nu_par != Inf) nu_par = exp(empirical_bayes$par[p_len+7])
  else nu_par = 999999999
  
  ans = fitGrove(W$D, 
                 XX, 
                 p = p, 
                 tau_par = tau_par, 
                 eta_par = eta_par, 
                 gamma_par = gamma_par, 
                 init_state = init_state,
                 nu = nu_par, 
                 sigma = sigma_par, 
                 alpha = alpha_par, 
                 beta = beta_par, 
                 n_samp = n_samp, 
                 transition_mode = transition_mode)
  ans$data$formula = formula
  ans$data$X = X
  ans$data$W = W
  colnames( ans$data$design_matrix ) = colnames(XX)
  ans$C_hat = C_hat
  ans$hyperparameters = list(  tau = tau_par, 
                               eta = eta_par, 
                               gamma = gamma_par, 
                               sigma = sigma_par, 
                               nu = nu_par, 
                               alpha = alpha_par, 
                               beta = beta_par)
  

  class(ans) = "grove"
  
  return( ans )
}

.groveEB.fixed.eta <- function(W, 
                              formula, 
                              X, 
                              alpha = 0.5,
                              beta = 1,
                              nu = 5, 
                              eta, 
                              n_samp = 500, 
                              verbose = FALSE, 
                              transition_mode = "Markov", 
                              method = "Nelder-Mead")
{
  if(class(W)!="DWT")
  {
    print("ERROR: W should be a DWT object")
    return(0)
  }
  
  if(strsplit(as.character(formula)[2],split="")[[1]][1] != 1)
  {
    print("ERROR: formula should include an intercept.")
    return(0)
  }
  
  if (transition_mode == "Markov") transition_mode = 1
  else if (transition_mode == "Independent") transition_mode = 0
  else {
    transition_mode = 1
    print("WARNING: Unrecognized transition mode. Default to Markov.")
  }
  
  XX = model.matrix(formula, X)
  p = table(attr(XX,"assign"))  
  p_len = length(p)
  
  frml = paste("C", paste(formula, sep="", collapse=" "), sep = " ")
  C_hat = lm(frml, data = cbind(C=W$C,X))$coefficients
  
  init_state = rep(0,p_len)
  sigma_par = mad(tail(t(W$D),n=1))
  tau_par = rep( (sd( head(t(W$D),n=1))/sigma_par)^2,p_len)
  nu_par = nu
  alpha_par = alpha
  beta_par = beta
  gamma_par = 0.3
  
  eta.rho = eta[1]

  if (p_len > 1) {
    eta.kap = eta[2]
 
  }
  else {
    eta.kap = NA
  }
  
  fr <- function(input) 
  { 
    tau = input[1:p_len]
    sigma = input[p_len+1]
    alpha = input[p_len+2]
    gamma.rho = input[p_len+3]
    gamma.kap = input[p_len+4]
    
    if(nu_par != Inf) {
      nu = input[p_len+5]
    } else nu = log(999999999)
    
    output = fitGroveML(W$D, 
                        XX, 
                        p = p, 
                        tau_par = exp(tau), 
                        eta_par = c(eta.rho, rep(eta.kap,p_len-1)), 
                        gamma_par = c(0.5*exp(gamma.rho)/(1+exp(gamma.rho)),
                                      rep(0.5*exp(gamma.kap)/(1+exp(gamma.kap)),p_len-1)),
                        init_state = init_state,
                        nu = exp(nu), 
                        sigma = exp(sigma), 
                        alpha = alpha, 
                        beta = beta, 
                        transition_mode = transition_mode)
    if(verbose == TRUE)
      print(-output$marginal_likelihood)
    return(-output$marginal_likelihood)
  }  
  
  if (nu != Inf) {
    
    empirical_bayes = optim(par=c(log(tau_par), 
                                  log(gamma_par/(0.5-gamma_par)), 
                                  log(gamma_par/(0.5-gamma_par)),
                                  log(sigma_par),
                                  alpha_par,
                                  log(nu_par)), 
                            fn=fr,  
                            method = method)
  }
  
  else {
    empirical_bayes = optim(par=c(log(tau_par), 
                                  log(gamma_par/(0.5-gamma_par)), 
                                  log(gamma_par/(0.5-gamma_par)),
                                  log(sigma_par),
                                  alpha_par), 
                            fn=fr,  
                            method = "Nelder-Mead")
  }
  
   
  tau_par = exp(empirical_bayes$par[1:p_len])
  eta_par = c( eta.rho, rep(eta.kap,p_len-1))
  gamma_par = c( 0.5*exp(empirical_bayes$par[p_len+1])/(1+exp(empirical_bayes$par[p_len+1])),
                 rep(0.5*exp(empirical_bayes$par[p_len+2])/(1+exp(empirical_bayes$par[p_len+2])),p_len-1))
  sigma_par = exp(empirical_bayes$par[p_len+3])
  alpha_par = empirical_bayes$par[p_len+4]
  beta_par = beta
  
  if (nu_par != Inf) nu_par = exp(empirical_bayes$par[p_len+5])
  else nu_par = 999999999
  
  ans = fitGrove(W$D, 
                 XX, 
                 p = p, 
                 tau_par = tau_par, 
                 eta_par = eta_par, 
                 gamma_par = gamma_par, 
                 init_state = init_state,
                 nu = nu_par, 
                 sigma = sigma_par, 
                 alpha = alpha_par, 
                 beta = beta_par, 
                 n_samp = n_samp, 
                 transition_mode = transition_mode)
  ans$data$formula = formula
  ans$data$X = X
  ans$data$W = W
  colnames( ans$data$design_matrix ) = colnames(XX)
  ans$C_hat = C_hat
  ans$hyperparameters = list(  tau = tau_par, 
                               eta = eta_par, 
                               gamma = gamma_par, 
                               sigma = sigma_par, 
                               nu = nu_par, 
                               alpha = alpha_par, 
                               beta = beta_par)
  
  
  class(ans) = "grove"
  
  return( ans )
}


.groveEB.fixed.eta.and.gamma <- function(W, 
                                        formula, 
                                        X, 
                                        alpha = 0.5,
                                        beta = 1, 
                                        nu = 5, 
                                        eta, 
                                        gamma, 
                                        n_samp = 500, 
                                        verbose = FALSE, 
                                        transition_mode = "Markov", 
                                        method = "Nelder-Mead")
{
  if(class(W)!="DWT")
  {
    print("ERROR: W should be a DWT object")
    return(0)
  }
  
  if(strsplit(as.character(formula)[2],split="")[[1]][1] != 1)
  {
    print("ERROR: formula should include an intercept.")
    return(0)
  }
  
  if (transition_mode == "Markov") transition_mode = 1
  else if (transition_mode == "Independent") transition_mode = 0
  else {
    transition_mode = 1
    print("WARNING: Unrecognized transition mode. Default to Markov.")
  }
  
  XX = model.matrix(formula, X)
  p = table(attr(XX,"assign"))  
  p_len = length(p)
  
  frml = paste("C", paste(formula, sep="", collapse=" "), sep = " ")
  C_hat = lm(frml, data = cbind(C=W$C,X))$coefficients
  
  init_state = rep(0,p_len)
  sigma_par = mad(tail(t(W$D),n=1))
  tau_par = rep( (sd( head(t(W$D),n=1))/sigma_par)^2,p_len)
  nu_par = nu
  alpha_par = alpha
  beta_par = beta
  
  eta.rho = eta[1]
  gamma.rho = gamma[1]
  
  if (p_len > 1) {
    eta.kap = eta[2]
    gamma.kap = gamma[2]
  }
  else {
    eta.kap = NA
    gamma.kap = NA
  }
  
  fr <- function(input) 
  { 
    tau = input[1:p_len]
    sigma = input[p_len+1]
    alpha = input[p_len+2]
    
    if(nu_par != Inf) {
      nu = input[p_len+3]
    } else nu = log(999999999)
    
    output = fitGroveML(W$D, 
                        XX, 
                        p = p, 
                        tau_par = exp(tau), 
                        eta_par = c(eta.rho, rep(eta.kap,p_len-1)), 
                        gamma_par = c(gamma.rho,rep(gamma.kap,p_len-1)),
                        init_state = init_state,
                        nu = exp(nu), 
                        sigma = exp(sigma), 
                        alpha = alpha, 
                        beta = beta, 
                        transition_mode = transition_mode)
    if(verbose == TRUE)
      print(-output$marginal_likelihood)
    return(-output$marginal_likelihood)
  }  
  
  if (nu != Inf) {
    
    empirical_bayes = optim(par=c(log(tau_par), 
                                  log(sigma_par),
                                  alpha_par,
                                  log(nu_par)), 
                            fn=fr,  
                            method = method)
  }
  
  else {
    empirical_bayes = optim(par=c(log(tau_par), 
                                  log(sigma_par),
                                  alpha_par), 
                            fn=fr,  
                            method = method)
  }
  
  
  tau_par = exp(empirical_bayes$par[1:p_len])
  eta_par = c( eta.rho, rep(eta.kap,p_len-1))
  gamma_par = c( gamma.rho,rep(gamma.kap,p_len-1))
  sigma_par = exp(empirical_bayes$par[p_len+1])
  alpha_par = empirical_bayes$par[p_len+2]
  beta_par = beta
  
  if (nu_par != Inf) nu_par = exp(empirical_bayes$par[p_len+3])
  else nu_par = 999999999
  
  ans = fitGrove(W$D, 
                 XX, 
                 p = p, 
                 tau_par = tau_par, 
                 eta_par = eta_par, 
                 gamma_par = gamma_par, 
                 init_state = init_state,
                 nu = nu_par, 
                 sigma = sigma_par, 
                 alpha = alpha_par, 
                 beta = beta_par, 
                 n_samp = n_samp, 
                 transition_mode = transition_mode)
  ans$data$formula = formula
  ans$data$X = X
  ans$data$W = W
  colnames( ans$data$design_matrix ) = colnames(XX)
  ans$C_hat = C_hat
  ans$hyperparameters = list(  tau = tau_par, 
                               eta = eta_par, 
                               gamma = gamma_par, 
                               sigma = sigma_par, 
                               nu = nu_par, 
                               alpha = alpha_par, 
                               beta = beta_par)
  
  
  class(ans) = "grove"
  
  return( ans )
}
