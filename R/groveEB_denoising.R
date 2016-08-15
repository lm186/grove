.groveEB.denoising <- function(W, 
                               formula, 
                               X, 
                               alpha = 0.5, 
                               beta = 1, 
                               nu = 5, 
                               n_samp = 500, 
                               verbose = FALSE,
                               transition_mode = "Markov",
                               method = "Nelder-Mead") {
  if(class(W)!="DWT") {
    print("ERROR: W should be a DWT object")
    return(0)
  }

  if(strsplit(as.character(formula)[2],split="")[[1]][1] != 1)  {
    print("ERROR: formula should include an intercept.")
    return(0)
  }

  if (transition_mode == "Markov") {
    transition_mode = 1
  } else if (transition_mode == "Independent") {
    transition_mode = 0
  } else {
    transition_mode = 1
    print("WARNING: Unrecognized transition mode. Default to Markov.")
  }
  
  XX = model.matrix(formula, X)
  p = table(attr(XX,"assign"))  
  p_len = length(p)
  
  frml = paste("C", paste(formula, sep="", collapse=" "), sep = " ")
  C_hat = lm(frml, data = cbind(C=W$C,X))$coefficients
  
  init_state = rep(0,p_len)
  sigma_par = mad(tail(t(W$D),n=10))
  eta_par = c(2, rep(0.1,max(p_len-1,1)))
  gamma_par = 0.3
  tau_par = rep( (sd( head(t(W$D),n=10))/sigma_par)^2,p_len)
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
    beta = input[p_len+7]
    ## If the nu_par is set to infinity, fit the model with stationary noise;
    ## otherwise take in the initial nu value for EB
    if(nu_par != Inf) {
      nu = input[p_len+8]
    }
    else {
      nu = log(999999999)
    }
    
   
    output = fitGroveML(W$D, 
                        XX, 
                        p = p, 
                        tau_par = exp(tau), 
                        eta_par = c(exp(eta.rho),
                                    rep(exp(eta.kap),p_len-1)), 
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
  
    empirical_bayes = optim(par = c(log(tau_par), 
                                    log(eta_par[1]), 
                                    log(eta_par[2]),
                                    log(gamma_par/(0.5-gamma_par)), 
                                    log(gamma_par/(0.5-gamma_par)),
                                    log(sigma_par),
                                    alpha_par,
                                    beta_par,
                                    log(nu_par)), 
                            fn=fr,  
                            method = method)
  }
  
  else {
    empirical_bayes = optim(par = c(log(tau_par),
                                    log(eta_par[1]), 
                                    log(eta_par[2]),
                                    log(gamma_par/(0.5-gamma_par)), 
                                    log(gamma_par/(0.5-gamma_par)),
                                    log(sigma_par),
                                    alpha_par,
                                    beta_par), 
                            fn=fr,  
                            method = "Nelder-Mead")
  }
  
  tau_par = exp(empirical_bayes$par[1:p_len])
  eta_par = c(exp(empirical_bayes$par[p_len+1]), 
              rep(exp(empirical_bayes$par[p_len+2]), p_len-1))
  gamma_par = c(0.5*exp(empirical_bayes$par[p_len+3])/(1+exp(empirical_bayes$par[p_len+3])),
                rep(0.5*exp(empirical_bayes$par[p_len+4])/(1+exp(empirical_bayes$par[p_len+4])), p_len-1))
  sigma_par = exp(empirical_bayes$par[p_len+5])
  alpha_par = empirical_bayes$par[p_len+6]
  beta_par = empirical_bayes$par[p_len+7]
  if (nu_par != Inf) nu_par = exp(empirical_bayes$par[p_len+8])
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
  
  ## Get the exact posterior mean of the z_{j,k}'s
  pmaps = lapply(ans$samples$S,function (x) {x[2,]})

  J = W$J
  tau.j.vec = rep(2^(-alpha_par*(0:(J-1))),2^(0:(J-1)))*ans$hyperparameters$tau[1]
  D.post.mean.z = 1/(1+1/(t(tau.j.vec)*nrow(W$D))) * t(unlist(pmaps)) * ans$data$W$D
  
  ans$D.post.mean.z = D.post.mean.z
  class(ans) = "grove"
  
  return( ans )
}


.groveEB.denoising.fixed.beta <- function(W, 
                                          formula, 
                                          X, 
                                          alpha = 0.5, 
                                          beta = 1, 
                                          nu = 5, 
                                          n_samp = 500, 
                                          verbose = FALSE,
                                          transition_mode = "Markov",
                                          method = "Nelder-Mead") {
  
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
  sigma_par = mad(tail(t(W$D),n=10))
  eta_par = c(2, rep(0.1,max(p_len-1,1)))
  gamma_par = 0.3
  tau_par = rep( (sd( head(t(W$D),n=10))/sigma_par)^2,p_len)
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


    ## If the nu_par is set to infinity, fit the model with stationary noise;
    ## otherwise take in the initial nu value for EB
    
    if(nu_par != Inf) {
      nu = input[p_len+7]
    }
    else {
      nu = log(999999999)
    }
    
   
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
  
    empirical_bayes = optim(par = c(log(tau_par),
                                    log(eta_par[1]), 
                                    log(eta_par[2]),
                                    log(gamma_par/(0.5-gamma_par)), 
                                    log(gamma_par/(0.5-gamma_par)),
                                    log(sigma_par),
                                    alpha_par,
                                    log(nu_par)), 
                            fn = fr,  
                            method = method)
  } else {
    empirical_bayes = optim(par = c(log(tau_par), 
                                    log(eta_par[1]), 
                                    log(eta_par[2]),
                                    log(gamma_par/(0.5-gamma_par)), 
                                    log(gamma_par/(0.5-gamma_par)),
                                    log(sigma_par),
                                    alpha_par), 
                            fn = fr,  
                            method = method)
  }
  
  tau_par = exp(empirical_bayes$par[1:p_len])
  eta_par = c(exp(empirical_bayes$par[p_len+1]), 
              rep(exp(empirical_bayes$par[p_len+2]), p_len-1))
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
  
  ## Get the exact posterior mean of the z_{j,k}'s
  pmaps = lapply(ans$samples$S,function (x) {x[2,]})

  J = W$J
  tau.j.vec = rep(2^(-alpha_par*(0:(J-1))),2^(0:(J-1)))*ans$hyperparameters$tau[1]
  D.post.mean.z = 1/(1+1/(t(tau.j.vec)*nrow(W$D))) * t(unlist(pmaps)) * ans$data$W$D
  
  ans$D.post.mean.z = D.post.mean.z
  class(ans) = "grove"
  
  return( ans )
}

.groveEB.denoising.fixed.alpha.beta <- function(W, 
                                                formula, 
                                                X, 
                                                alpha = 0.5, 
                                                beta = 1, 
                                                nu = 5, 
                                                n_samp = 500, 
                                                verbose = FALSE,
                                                transition_mode = "Markov",
                                                method = "Nelder-Mead") {
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
  sigma_par = mad(tail(t(W$D),n=10))
  eta_par = c(2, rep(0.1,max(p_len-1,1)))
  gamma_par = 0.3
  tau_par = rep( (sd( head(t(W$D),n=10))/sigma_par)^2,p_len)
  nu_par = nu
  
  fr <- function(input) 
  { 
    tau = input[1:p_len]
    eta.rho = input[p_len+1]
    eta.kap = input[p_len+2]
    gamma.rho = input[p_len+3]
    gamma.kap = input[p_len+4]
    sigma = input[p_len+5]
    
    ## If the nu_par is set to infinity, fit the model with stationary noise;
    ## otherwise take in the initial nu value for EB
    if(nu_par != Inf) {nu = input[p_len+6]} else nu = log(999999999)
    
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
                        beta = beta
                        transition_mode = transition_mode)
    if(verbose == TRUE)
      print(-output$marginal_likelihood)
    return(-output$marginal_likelihood)
  }  
  
  if (nu != Inf) {
  
    empirical_bayes = optim(par=c(log(tau_par), 
                                log(eta_par[1]), log(eta_par[2]),
                                log(gamma_par/(0.5-gamma_par)), 
                                log(gamma_par/(0.5-gamma_par)),
                                log(sigma_par),log(nu_par)), 
                            fn=fr,  
                            method = method)
  }
  
  else {
    empirical_bayes = optim(par=c(log(tau_par), 
                                  log(eta_par[1]), log(eta_par[2]),
                                  log(gamma_par/(0.5-gamma_par)), 
                                  log(gamma_par/(0.5-gamma_par)),
                                  log(sigma_par)), 
                            fn=fr,  
                            method = method)
  }
  
  tau_par = exp(empirical_bayes$par[1:p_len])
  eta_par = c(exp(empirical_bayes$par[p_len+1]), 
              rep(exp(empirical_bayes$par[p_len+2]),p_len-1))
  gamma_par = c(0.5*exp(empirical_bayes$par[p_len+3])/(1+exp(empirical_bayes$par[p_len+3])),
                rep(0.5*exp(empirical_bayes$par[p_len+4])/(1+exp(empirical_bayes$par[p_len+4])),p_len-1))
  sigma_par = exp(empirical_bayes$par[p_len+5])
  if (nu_par != Inf) nu_par = exp(empirical_bayes$par[p_len+6])
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
                 alpha = alpha, 
                 beta = beta, 
                 n_samp = n_samp,
                 transition_mode = transition_mode))
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
                               alpha = alpha, 
                               beta = beta )
  
  ## Get the exact posterior mean of the z_{j,k}'s
  pmaps = lapply(ans$samples$S,function (x) {x[2,]})

  J = W$J
  tau.j.vec = rep(2^(-alpha*(0:(J-1))),2^(0:(J-1)))*ans$hyperparameters$tau[1]
  D.post.mean.z = 1/(1+1/(t(tau.j.vec)*nrow(W$D))) * t(unlist(pmaps)) * ans$data$W$D
  
  ans$D.post.mean.z = D.post.mean.z
  class(ans) = "grove"
  
  return( ans )
}
