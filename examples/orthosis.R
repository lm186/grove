library(wavethresh)
library(grove)
library(fda)

original.data = read.table("~/Dropbox/work/wavelets/data/orthosis.txt", header=TRUE)

condition.subset = c(1,2,3,4)
data.tmp = original.data[is.element(original.data$Condition,condition.subset),]


n_subjects = length(unique(data.tmp$Subject))
n_conditions = length(unique(data.tmp$Condition))
n_times = sum(data.tmp$Subject==1 & data.tmp$Replication==1 & data.tmp$Condition==condition.subset[1])
n_replicates = length(unique(data.tmp$Replication))



data = matrix(NA, nrow=n_subjects*n_conditions*n_replicates, ncol=n_times)

for(i in 1:n_conditions)
{
  for(j in 1:n_subjects)
  {
    indices = ((i-1)*n_replicates*n_subjects + (j-1)*n_replicates + 1 ):((i-1)*n_replicates*n_subjects + (j-1)*n_replicates + n_replicates)
    data[indices, ] = matrix(data.tmp$Moment[which(data.tmp$Subject==j & data.tmp$Condition==condition.subset[i])], nrow=n_replicates, byrow=TRUE)
    
  }
}

Y = data
W = DWT(Y)

D = W$D
C = W$C

XX = data.frame(data.tmp[data.tmp$Time==0.0020,c("Condition","Subject")])
XX$Condition = factor(XX$Condition,levels=intersect(c(1,2,3,4),condition.subset))
XX$Subject = factor(XX$Subject)

frml = formula( ~ 1 + Condition + Subject)


system.time({ans = groveEB.fixed.eta.and.gamma(W, frml, XX, alpha = 0.5, beta = 1, eta=c(0.3,0.3),gamma=c(0.4,0.4),nu = 5, n_samp = 500, transition_mode = "Markov", verbose = FALSE)})

ans$hyper


x11
par(mfrow=c(1,2))
plotStates(ans, block = "Condition", legend = TRUE)
plotStates(ans, block = "Subject", legend = TRUE)


contrastlist = list(c(1,0,0,0,0,0,0,0,0,0),c(0,1,0,0,0,0,0,0,0,0),c(0,0,1,-1,0,0,0,0,0,0),c(0,-1,1,0,0,0,0,0,0,0),c(0,-1,1,1,0,0,0,0,0,0))
contrastnames = c("Intercept","Orthosis vs Control","Spring 1 vs Spring 2","Spring 1 vs Orthosis","(Spring 1 + Spring 2) - (Orthosis + Control)")


x11()
par(mfrow=c(3,2))
for (i in c(2,3,5)) {
  aa = invDWT( ans, x = contrastlist[[i]], include_C = TRUE, sample_C=TRUE )
  bb = invDWT( ans, x = contrastlist[[i]], include_C = FALSE, sample_C=FALSE )
  ylim = range(c(aa,bb))

  for (include_C in c(FALSE,TRUE)) {
    aa = invDWT( ans, x = contrastlist[[i]], include_C = include_C, sample_C=include_C )
    plotFun(aa,main = paste(contrastnames[i],"--",c("including scaling coefficient","by mother coefficients")[2-as.integer(include_C)]),ylim=ylim)
    abline(h=0,col='red',lty="dashed")
  }
}




## visualize the data
x11()
condition.vec = c("Control","Orthosis","Spring 1","Spring 2")
par(mfcol=c(7,4))
for (condition in 1:4) {
  for (subject in 1:7) {
    for (replicate in 1:10) {
      if (replicate == 1) {xlim=c(1,256);ylim=range(Y[XX$Condition==condition & XX$Subject==subject,])}
      if (replicate > 1) par(new=TRUE)
      plot(1:256,Y[XX$Condition==condition & XX$Subject==subject,][replicate,],type='l',col='gray',xlim=xlim,ylim=ylim,xlab="t",ylab="Moment",main=paste(condition.vec[condition],"- Subject",subject))
      
    }
  }
}




