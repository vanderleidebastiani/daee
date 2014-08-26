PVR<-function(traits,dist,cumulative=0.99,VMoran=0.025,pMoran=0.05,check=TRUE){
	if((0<=cumulative & cumulative<=1)==FALSE){
		stop("\n Cumulative percentage must be higher than the cumulative percentage of the first two eigenvalues, and less than 1\n")
	}
	traits<-as.matrix(traits)
	dist<-as.matrix(dist)
	if (check==TRUE){
		if (is.null(rownames(traits))) {
            stop("\n Error in rownames of traits\n")
        }
        if (is.null(colnames(dist))) {
            stop("\n Error in colnames of dist\n")
        }
        if (is.null(rownames(dist))) {
            stop("\n Error in rownames of dist\n")
        }			
		chec<-rownames(traits)==colnames(dist)
		if (length(which(chec==FALSE))!=0){
			stop("\n Matrices labels do not match. Check data\n")
		}
	}
	ordination<-wcmdscale(sqrt(dist),eig=TRUE)
	values<-ordination$eig[which((ordination$eig>=0)==TRUE)]
	vectors<-ordination$points
	colnames(vectors)=colnames(vectors,do.NULL=FALSE,prefix="Axis.")
	relative<-values/sum(values)
	if(cumulative<sum(relative[1:2])){
		print(paste("Relative eigenvalue for axis 1 =",relative[1]))
		print(paste("Relative eigenvalue for axis 2 =",relative[2]))
		stop("\n Cumulative percentage must be higher than the cumulative percentage of the first two eigenvalues\n")
	}
	cumulative2<-as.vector(rep(NA,length(values)))
	for (i in 1:length(values)){
		cumulative2[i]<-sum((values/sum(values))[1:i])
	}
	use<-which((cumulative2<=cumulative)==TRUE)
	if(cumulative==1){
		use<-1:(length(values)-1)
	}
	values2<-cbind(values,relative,cumulative2)
	colnames(values2)=c("Eigenvalues","Relative_eig","Cumul_eig")
	rownames(values2)=1:length(values)
	x<-vectors[,use]
	vare<-dim(traits)[2]
	fac<-length(use)
	xnam <- paste("x[,", 1:fac,"]", sep="")
	result<-matrix(NA,nrow=vare,ncol=fac)
		for(i in 1:vare){
			y<-traits[,i]
			for (j in 1:fac){
				result[i,j]<-as.numeric(summary(lm(as.formula(paste("y ~ ", paste(xnam[1:j], collapse= "+")))))$r.squared)
				}	
		}
	colnames(result)=colnames(vectors)[1:fac]
	rownames(result)=colnames(traits,do.NULL=FALSE,prefix="Trait.")
	resFIN<-matrix(NA,vare,6)
		for(i in 1:vare){
			y<-traits[,i]
			MOR=NULL
			xnam2<-xnam
			possi<-1:fac
			possi2<-1:fac
			f<-character()		
			res.MORAN<-matrix(NA,fac,2)
				repeat{
					result.moran<-matrix(NA,2,fac)
					for (j in possi2){
						regre<-lm(as.formula(paste("y ~ ",paste(f,sep="",collapse=" + "),paste("+",xnam[j]),collapse=" + ")))
						moran<-Moran.I(regre$residuals,dist,scaled=TRUE)
						result.moran[1,j]<-moran$observed
						result.moran[2,j]<-moran$p.value
						}
					colnames(result.moran)=possi
					menor<-as.numeric(noquote(names(sort(result.moran[1,],decreasing=FALSE))[1]))
					xnam2<-xnam2[-which(xnam2==xnam[menor])]
					MOR<-result.moran[1,menor]
					pMOR<-result.moran[2,menor]
					f<-c(f,paste(xnam[menor]))
					res.MORAN[length(f),1]<-MOR
					res.MORAN[length(f),2]<-pMOR
					possi2<-possi2[-which(possi2==possi[menor])]			
					if (length(f)==fac) break
					if ((MOR<=VMoran)==TRUE){
						if (pMOR>=pMoran) break
					}
				}
			minMORAN<-length(f)
			regre2<-lm(as.formula(paste("y ~ ",paste(f[1:minMORAN],sep="",collapse=" + "))))
			resFIN[i,1]<-paste("y ~ ",paste(f[1:minMORAN],sep="",collapse=" + "))
			resFIN[i,2]<-minMORAN
			moran2<-Moran.I(regre2$residuals,dist,scaled=TRUE)
			resFIN[i,3]<-moran2$observed
			resFIN[i,4]<-moran2$p.value
			summary.regre2<-summary(regre2)
			resFIN[i,5]<-as.numeric(summary.regre2$r.squared)
			resFIN[i,6]<-as.numeric(pf(summary.regre2$fstatistic[1],summary.regre2$fstatistic[2],summary.regre2$fstatistic[3],lower.tail=FALSE))
		}
	colnames(resFIN)=c("Parameters","N.Parameters","Moran.Obs","Moran.p","R.Squared","p")
	rownames(resFIN)=colnames(traits,do.NULL=FALSE,prefix="Trait.")
	res.x<-t(as.matrix(values2[1:fac,3]))
	colnames(res.x)=colnames(res.x,do.NULL=FALSE,prefix="Axis.")
	rownames(res.x)="Cumulative"
	inf<-res.x[1,fac]
	RES<-list(values=values2,vectors=vectors,inf.cumulative=inf,n.axis.considered=fac,moran.less.than=VMoran,p.moran.greater.than=pMoran,PSR.curve.axis.x=res.x,PSR.curve.axis.y=result,minimum.moran=resFIN)
	class(RES)<-"pvr"
	return(RES)
}