mantel.residuals<-function(A,B){
	a<-as.vector(as.dist(A))
	b<-as.vector(as.dist(B))
	modelo<-lm(a~b)
	RESIDUOS<-residuals(modelo)
	N<-dim(as.matrix(A))[1]
	RES<-matrix(NA,N,N)
	OR<-0
	for(j in 1:N){
		for(i in 1:N){
			if(i==j){
				RES[i,j]=0
			}
			if(i>j){
				OR<-OR+1
				RES[i,j]=RESIDUOS[OR]
			}	
		}	
	}
	R<-sqrt(summary(modelo)$r.squared)
	if(modelo$coefficients[2]<0){
		R<-R*-1
	}
	Dis_RES<-as.dist(RES,diag=T)
return(list(Mantel_statistic_r=R,Residuals=Dis_RES))
}