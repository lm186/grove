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
                    band.type = "pointwise",
                    main = "", 
                    col = "blue", 
                    type = 'l', 
                    ylab = "", 
                    xlab = "", 
                    ylim = NULL,
                    return.range = FALSE) {
  
  if (length(p) != 3) {
    stop("p should have 3 elements.")
  }
  
  if (!(band.type %in% c("global", "pointwise", "both"))) {
    warning("band.type not recognized. Used default value.")
    band.type <- "pointwise"
  }
  # compute point-wise ci
  point.wise.cis <- apply(data, 2, function(x) quantile(x, probs = p))
  # compute global ci
  if (band.type %in% c("global", "both")) {
    alpha <- p[1] + 1 - p[3]
    # For each time t, compute the distance of each curve normalized 
    # wrt to variance at t.
    dist.matrix <- apply(data, 2, function(x) {
      ((x - mean(x)) / sd(x)) ^ 2
    })
    # Identify curves that are further from the mean based on the notion
    # of distance defined above.
    n <- round(nrow(data) * alpha)
    removed.curves <- rep(NA, n)
    for (i in 1:n) {
      worst.elem <- which(dist.matrix == max(dist.matrix), arr.ind = TRUE) 
      removed.curves[i] <- worst.elem[1, 1]
      # Set distance to zero for the entire curve, otherwise the same
      # curve could be selected multiple times.
      dist.matrix[worst.elem[1, 1], ] <- 0
    }  
    # Compute range for the remaining curves.
    global.cis <- apply(data[-removed.curves, ], 2, function(x) { range(x) })
  } 
  
  x <- seq(1, ncol(point.wise.cis))
  if (is.null(ylim)) {
    if (band.type == "pointwise") {
      ylim <- range(point.wise.cis)
    } else {
      ylim <- range(global.cis)
    }
  }
  plot(x, point.wise.cis[2, ], 
       col = col, 
       type = type, 
       ylim = ylim, 
       main = main, 
       ylab = ylab, 
       xlab = xlab)
  if (band.type %in% c("global", "both")) {    
    polygon(c(x, rev(x)), 
            c(global.cis[2, ], rev(global.cis[1, ])), 
            col = "lightgrey", 
            border = "grey", 
            lwd = 2)
  } 
  if (band.type %in% c("both", "pointwise")) {    
    polygon(c(x, rev(x)), 
            c(point.wise.cis[3, ], rev(point.wise.cis[1, ])), 
            col = "grey", 
            border = "darkgrey", 
            lwd = 2)
  } 
  lines(x, point.wise.cis[2, ], col = "blue", lwd = 2)
  if (return.range) {
    if (band.type %in% c("global", "both")) {    
      return(range(global.cis))
    } else {
      return(range(point.wise.cis))
    }    
  }
}