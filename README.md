Wavelet functional ANOVA through Markov groves
==============================================

Functional denoising and functional ANOVA through wavelet-domain Markov groves.

### Install
The package can be installed on Linux and Mac using `devtools`:

```S
install.packages('devtools')
library('devtools')
devtools::install_github('grove', 'jacsor')
```

### Use
There are 7 functions in this package, and their descriptions are provided in the help files.

```S
W <- DWT(data)
ans <- Denoise(W)
ans <- FAnova(W, X, ~ 1 + factorA + factorB)
data <- GenerateSyntheticAnova()
denoised.data <- InvDWT(ans)
PlotFun(denoised.data)
PlotStates(ans, block = "factorA")
```

### Reference
Ma L. and Soriano J. (2017). Efficient functional ANOVA through wavelet-domain 
Markov groves. JASA (To appear).
