plot.pvr<-function(x, trait = 1, ...){
	plot(x$PSR.curve.axis.x,x$PSR.curve.axis.y[trait,],ylim=c(0,1),xlim=c(0,1),xlab="Cumulative eigenvalues",ylab="R2")
	abline(0,1)	
}