#' @title Residuals from Mantel test
#' 
#' @description Function to extract residuals from Mantel test.
#' 
#' @encoding UTF-8
#' @importFrom stats lm residuals summary.lm as.dist
#' @aliases mantel.residuals
#' @param x Dissimilarity matrices of dist class.
#' @param y Dissimilarity matrices of dist class.
#' @return \item{statistic}{The Mantel statistic.} \item{residuals}{Residuals 
#' extracted from the Mantel test (Class dist).}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords daee
#' @examples
#' 
#' #require(vegan)
#' #A<-vegdist(matrix(sample(1:10, 20, replace  =TRUE), 5, 5))
#' #B<-vegdist(matrix(sample(1:10, 30, replace = TRUE), 5, 6))
#' #mantel(A, B)
#' #mantel.residuals(A, B)
#' 
#' @export
mantel.residuals<-function(y, x){
	res <- list(call = match.call())
	y<-stats::as.dist(y)
	x<-stats::as.dist(x)
	N<-attr(y,"Size")
	if(N != attr(x,"Size")){
		stop("\n Incompatible dimensions\n")
	}
	model<-stats::lm(as.vector(y)~as.vector(x))
	res.resid<-stats::residuals(model)
	RES<-matrix(NA, N, N)
	m<-0
	for(j in 1:N){
		for(i in j:N){
			if(i==j){
				RES[i,j]<-0
			} else {
				m<-m+1
				RES[i,j]<-res.resid[m]
			}
		}
	}
	R2<-sqrt(stats::summary.lm(model)$r.squared)
	if(model$coefficients[2]<0){
		R2<-R2*-1
	}
	res$statistic<-R2
	res$residuals<-stats::as.dist(RES, diag = TRUE)
return(res)
}