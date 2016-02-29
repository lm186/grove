

library(wavethresh)
library(grove)

m=512

factor.A.1 <- DJ.EX(n=m, noisy=FALSE)$doppler
factor.A.2 <- factor.A.1
factor.A.2[45:60]  <- 1.5*factor.A.1[45:60] 

factor.B.1 <- DJ.EX(n=m, noisy=FALSE)$bumps
factor.B.2  <- factor.B.1 
factor.B.3  <- factor.B.1 
factor.B.2[310:360] = 0.5*factor.B.1[310:360]
## factor.B.2[45:60]  <- 3*factor.B.1[45:60] 
## factor.B.2[310:320] = factor.B.1[310:320] + 10
## factor.B.3[45:60]  <- 1.5*factor.B.1[45:60] 

# factor.C.1 <- DJ.EX(n=m, noisy=FALSE)$blocks
# factor.C.2 <- factor.C.1
# factor.C.2[45:60]  <- 1.5*factor.C.1[45:60] 

plot(1:512,factor.B.2,type="l")


###############################

st_dev = 2
# n = c(1,1,1,1,1)
n.replicate=5
sex = as.factor(c("M","F"))
snp = as.factor(c("aa","ab","bb"))
## ttt = as.factor(c("T","U"))
XX = expand.grid(factorA=sex,factorB=snp)## ,factorC=ttt)
XX = XX[rep(1:nrow(XX),each=n.replicate),]

n_tot = nrow(XX)

small.signal = rbind(factor.A.1, factor.A.2, factor.B.1, factor.B.2, factor.B.3)
                     ## factor.C.1,factor.C.2)
signal = matrix(NA, nrow = n_tot, ncol= m)
noisy.signal = matrix(NA, nrow = n_tot, ncol= m)
D = matrix(NA, nrow = n_tot, ncol= (m-1) )
C = rep(NA, n_tot)



## sex = as.factor( c(rep("M",3), rep("F",3)) )
## snp = as.factor( rep(c("aa","ab","bb"),2 ) )
## XX = data.frame(factorA = sex, factorB = snp)

frml = formula( ~ 1 + factorA + factorB) ## + factorC)
X = model.matrix(frml, XX)

for(i in 1:nrow(X)) { ## for each factor combination
##  for (r in 1:n.replicate) {
  signal[i,] = small.signal[1,] + ( small.signal[2,] - small.signal[1,] )*X[i,2] + 
    small.signal[3,] + ( small.signal[4,] - small.signal[3,] )*X[i,3] + 
    ( small.signal[5,] - small.signal[3,] )*X[i,4] ## + 
  ##  small.signal[6,] + ( small.signal[7,] - small.signal[6,] )*X[i,5]
  ## colSums(small.signal[which(X[i,2:4]==1),])
  noisy.signal[i,] = signal[i,] + rnorm(m,0,st_dev)
##   }

}

for (i in 1:nrow(X)) {
  plot(1:512,noisy.signal[i,],type="l",ylim=range(noisy.signal)) 
}

## range(noisy.signal)
## range(signal)



Y = noisy.signal


W = DWT(noisy.signal)

D = W$D
C = W$C

## X = X[,-c(2,4)]
## p=c(1,1,2)

#####

frml = formula( ~ 1 + factorA + factorB) ##  + factorC)

ans = groveEB(W, frml, XX, alpha = 0.5, beta = 1, nu = 5, n_samp = 500, verbose = TRUE)

  
plotStates(ans, block = "Intercept", legend = FALSE)
plotStates(ans, block = "factorA", legend = FALSE)
plotStates(ans, block = "factorB", legend = TRUE)
## plotStates(ans, block = "factorC", legend = TRUE)


aa = invDWT( ans, x = c(1,0,0,0) )
plotFun(aa)

aa = invDWT( ans, x = c(0,1,0,0) )
plotFun(aa)

aa = invDWT( ans, x = c(0,0,1,0) )
plotFun(aa)

aa = invDWT( ans, x = c(0,0,0,1) )
plotFun(aa)




