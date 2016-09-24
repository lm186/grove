#' Inverse discrete wavelet transform  
#'
#' This function performs the inverse discrete wavelet transform.
#'
#' @param grove.obj An object of class \code{grove}.
#' @param x A vector with the value of the predictor. 
#' @param include.C If \code{TRUE}, C is used to compute to reconstruct 
#' the function.
#' @param sample.C If \code{TRUE}, draws from C are used to recontruct 
#' the function.

#' @return A matrix each raw represents a draw from the reconstructed signal.
#' @export
#' @examples
#' data <- DJ.EX(n = 512, noisy = TRUE, rsnr = 5)$doppler
#' W <- DWT(data)
#' ans <- Denoise(W)
#' denoised.data <- invDWT(ans)
#' plot(data, type = "l")
#' lines(denoised.data[1, ], col = "red")

invDWT <- function(grove.obj, 
                   x = 1, # TODO: change for anova 
                   include.C = TRUE, 
                   sample.C = FALSE) {
  
  if (class(grove.obj) != "grove") {
    stop("Input should be a grove class object")
  }
  
  D <- grove.obj$samples$mean
  n_samp <- dim(D)[3]
  
  if (include.C) {
    C <- grove.obj$C_hat
  } else {
    C <- rep(0, length(grove.obj$C_hat))
  }
  
  m <- dim(D)[2] + 1 
  
  output <- matrix(NA, ncol = m, nrow = n_samp)
  temp <- wd(rep(1, m))
  
  y <- grove.obj$data$W$C
  X <- model.matrix(grove.obj$data$formula, grove.obj$data$X)
  
  s20 <- summary(lm(y ~ X))$sigma
  nu0 <- 10
  
  for (i in 1:n_samp) {
    if (grove.obj$data$X == 1) {
      temp$D <- rev(D[, , i])
    } else {
      temp$D <- rev(t(D[, , i]) %*% x)
    }
    
    if (include.C && sample.C) { 
      # draw C from its posterior
      
      # This code segment is from Peter Hoff's textbook 
      # "A first course in Bayesian Statistics"
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
      #
    }
    
    temp$C[length(temp$C)] <- sum(C * x)
    output[i, ] <- wr(temp)    
  }
  return(output)
}
