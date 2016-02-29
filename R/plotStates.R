#' Function to plot the hidden states
#'
#' This function plot on a tree the state of the latent variables.
#'
#' @param grove_obj
#' @param block
#' @param legend
#' @param main
#' @param prior
#' 
#' @return A plot.
#' @export
#' @examples
#' a <- 'TO DO'

plotStates <- function(grove_obj, 
                       block = "Intercept", 
                       legend = FALSE, 
                       main = NULL, 
                       prior = FALSE, ...)
{
  
  if(class(grove_obj)!="grove")
  {
    print("ERROR: input should be a grove class object")
    return(0)
  }
  
  temp = gsub(" + ", "+", as.character(grove_obj$data$formula)[2], fixed = TRUE)
  temp = strsplit(temp,split = "+", fixed = TRUE)[[1]]
  if(block == "Intercept")
    index = 1
  else
  {
    index = which(temp == block) 
    if(length(index)==0)
    {
      print("ERROR: no block with such a name")
      return(0)
    }    

  }


  p = grove_obj$data$p
  pos = rep(0,2^length(p))
  for(i in 0:(2^(length(p))-1) )
  {
    pos[i+1]=as.integer(intToBits(i))[index]
  }
  states = which(pos==1)

  
  if (prior) S = grove_obj$samples$Sprior ## plot prior state probabilities
  else S = grove_obj$samples$S ## plot posterior state probabilities
  
  J = length(S)
  colori = colorRampPalette(c("gray85", "darkblue"))( 100 )  
  idx = 2^(0:(J-1))
  if(legend)
    plot(1, type="n", axes=F, xlab="", ylab="scale", xlim=c(0,1), ylim=c(-2,J+2))
  else
    plot(1, type="n", axes=F, xlab="", ylab="scale", xlim=c(0,1), ylim=c(0,J+2))
  if (is.null(main)) title(main = block)
  else title(main = main)
  i = 0
  for(j in 0:(J-1))
  {
    if(j==0)
      points( seq(0,1,by=1/(2^j+1))[2:(2^j+1)], rep(J-j,2^j), col=colori[ rev(ceiling(sum(S[[j+1]][states])*99))], 
              pch=15, cex=1.5 )
    else
      points( seq(0,1,by=1/(2^j+1))[2:(2^j+1)], rep(J-j,2^j), col=colori[ rev(ceiling(colSums(S[[j+1]][states,])*99))], 
              pch=15, cex=1.5 )
    i = i + 2^j 
  }
  axis(2,  at=1:J, labels=rev(0:(J-1)), las=2)
  
  if(legend)
  {
    points( seq(0.2,0.95,length=11), rep(-.5,11), col=colori[round(seq(1,100,length=11))], 
            pch=15, cex=1.5 )
    text(seq(0.2,0.95,length=11), -0.7, labels=seq(0,1,length=11), pos=1 )
    rect(0.0, -2.2, 1, 0.1)
    if (prior) {
      if(block == "Intercept")
        text(0.1, -1.1 , expression(P(S[jk]==1)))
      else
        text(0.1, -1.1 , expression(P(R[jk]==1)))
    }
    else {
      if(block == "Intercept")
        text(0.1, -1.1 , expression(P(S[jk]==1~'|'~D)))
      else
        text(0.1, -1.1 , expression(P(R[jk]==1~'|'~D)))            
    }
    
  }
  
}

