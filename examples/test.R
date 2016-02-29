library(grove)

n = 100
J = 5
w = matrix( rnorm( n * ( 2^J - 1) ), nrow = n )
dims = c(1,2,3)

factor_a = factor(sample(1:dims[2], n, replace=TRUE))
factor_b = factor(sample(1:dims[3], n, replace=TRUE))

X = model.matrix(~ factor_a + factor_b)

tau_par = c(3,3,3)
eta_par = c(2,.1,.1)
gamma_par = c(.1,.1,.1)
init_state = rep(0,3) 
output = fitGrove(w, X, p = dims - c(0,1,1), tau_par = tau_par, 
                  eta_par = eta_par, gamma_par = gamma_par, 
                  init_state = init_state,
                  nu = 5, sigma = 10, alpha = 0.5, beta = 1)
class(output) = "grove"

output$prior_null

output$data$formula = formula( ~ 1 + factorA + factorB)
  
plotStates(output, block= "Intercept", legend = TRUE)
plotStates(output, block= "factorA", legend = TRUE)
plotStates(output, block= "factorB", legend = TRUE)


post_sample_means = invDWT(output)
post_sample_means

# plotStates(output$S, states = 3, main = "Intercept", legend = TRUE)
  

