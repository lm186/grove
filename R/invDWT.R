#' Inverse discrete wavelet transform  
#'
#' This function performs the inverse discrete wavelet transform.
#'
#' @param grove_obj An object of class \code{grove}.
#' @param x A vector with the value of the predictor. 
#' @param include_C If \code{TRUE}, C is used to compute to reconstruct the function.
#' @param sample_C If \code{TRUE}, draws from C are used to recontruct the function.

#' @return A matrix each raw represents a draw from the reconstructed signal.
#' @export
#' @examples
#' a <- 'TO DO'


invDWT <- function(grove_obj, 
                   x, 
                   include_C = TRUE, 
                   sample_C = FALSE) {
  
  if (class(grove_obj) != "grove") {
    print("ERROR: input should be a grove class object")
    return(0)
  }
  
  D <- grove_obj$samples$mean
  n_samp <- dim(D)[3]
  
  if (include_C) {
    C <- grove_obj$C_hat
  } else {
    C <- rep(0,length(grove_obj$C_hat))
  }
  
  m <- dim(D)[2]+1 
  
  output <- matrix(NA, ncol = m, nrow = n_samp)
  temp <- wd(rep(1, m))
  
  y <- grove_obj$data$W$C
  X <- model.matrix(grove_obj$data$formula,grove_obj$data$X)
  
  s20 <- summary(lm(y ~ X))$sigma
  nu0 <- 10
  
  for(i in 1:n_samp) {
    temp$D <- rev( t(D[,,i])%*%x )  
    
    if (include_C && sample_C) { ## draw C from its posterior
      
      ##### This code segment is from Peter Hoff's textbook "A first course in Bayesian Statistics"
      g <- length(y)
      S <- 1
      n <- dim(X)[1] 
      p <- dim(X)[2]
      Hg <- ( g /( g+1)) * X%*% solve( t (X)%*%X)%*%t (X)
      SSRg <- t(y)%*%(diag(1,nrow=n)-Hg)%*%y
      s2 <- 1/rgamma(S,(nu0+n )/2 , ( nu0* s20+SSRg)/2 )
      Vb <- g* solve(t(X)%*%X)/(g+1)
      Eb <- Vb%*%t (X)%*%y
      E <- matrix ( rnorm(S*p , 0 , sqrt(s2)),S,p)
      
      C <- t(t(E%*%chol (Vb) ) +c(Eb) )
      #####
    }
    
    temp$C[length(temp$C)] <- sum(C*x)
    output[i,] <- wr(temp)    
  }
  return(output)
}


# 
# invDWT2 <- function(C, D, m)
# {
#     output = matrix(NA, ncol=m, nrow=ncol(D))
#     temp = wd(rep(1,m))
#     for(i in 1:ncol(D))
#     {
#         temp$D = rev(D[,i])
#         temp$C[length(temp$C)] = C[i]
#         output[i,] = wr(temp)
#     }
#     return(output)
# }
# 
