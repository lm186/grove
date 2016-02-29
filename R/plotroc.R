.plotroc = function(p, p0, fpr = seq(0,1,by=0.01), col="black", xlim=c(0,1), ylim=c(0,1), main="", lty=1,lwd=2) { ## helper function to plot ROC
  plot(fpr,ecdf(p)(quantile(p0,fpr)), type="l", xlab="False postive rate", ylab="True positive rate", col=col, xlim=xlim, ylim=ylim, main=main, lty=lty,lwd=lwd)
}
