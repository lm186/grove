#' Function to plot the denoised signal
#'
#' This function plot the denoised signal
#'
#' @param data
#' @param p
#' @param main
#' @param col
#' @param type
#' @param xlab
#' @param ylab
#' @param ylim
#' 
#' @return A plot.
#' @export
#' @examples
#' a <- 'TO DO'

plotFun <- function(data, 
                    p = c(0.025, .5, 0.975), 
                    main = "", 
                    col = "blue", 
                    type = 'l', 
                    ylab = "", 
                    xlab = "", 
                    ylim = NULL,
                    ...)
{
  temp = apply(data, 2, function(x) quantile(x, probs=p))
  x = seq(1, ncol(temp))
  if (is.null(ylim)) ylim = range(temp)
  plot(x, temp[2,], col=col, type=type, ylim=ylim, main=main, ylab = ylab, xlab = xlab,...)
  polygon(c(x,rev(x)),c(temp[3,],rev(temp[1,])), col="grey", border = "darkgrey", lwd=2)
  lines(x, temp[2,], col="blue", lwd=2)
  return(range(temp))
}