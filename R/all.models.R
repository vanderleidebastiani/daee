all.models<-function(response,variables,subset,type ="AICc",only_intercept=FALSE,importance=TRUE){
	ImpValue=NULL
	response<-as.matrix(response)
	variables<-as.matrix(variables)
	colnames(variables)<-colnames(variables,do.NULL=FALSE,prefix="var_")
	dados<-as.data.frame(cbind(response,variables))
	colnames(dados)=c("y",colnames(variables))
	s<-1:subset
	m<-dim(variables)[2]
	nobs<-dim(variables)[1]
	possibilities<-factorial(m)/(factorial(s) * factorial(m - s))
	N_models<-sum(possibilities)
	RES <- vector("list", N_models)
	for(i in 1:subset){
		combinations <- combn(colnames(variables), i, simplify = TRUE)
		for (j in 1:possibilities[i]) {
			RES[[(j + sum(possibilities[1:i - 1]))]] <- glm(as.formula(paste("y ~ ", paste(combinations[, j], collapse= "+"))),data=dados)
		}
	}
	if(only_intercept==TRUE){
		RES[[N_models+1]]<-glm(y~1,data=dados)
		names(RES)=sprintf("Model_%.d",1:(N_models+1))
	}
	else{
		names(RES)=sprintf("Model_%.d",1:N_models)
	}
	TAB_IC<-ICtab(RES,type=type,nobs=nobs,base=T,weights=T,mnames=names(RES))
	class(TAB_IC)="data.frame"
	if(importance==TRUE){
		coef_uni<-colnames(variables)
		ImpValue<-matrix(0,length(coef_uni),1)
		rownames(ImpValue)=coef_uni
		colnames(ImpValue)="Relative_Importance"
		tab_ic<-ICtab(RES,type=type,nobs=nobs,weights=T,sort=FALSE,mnames=names(RES))
		for(i in 1:N_models){
			weight<-tab_ic$weight[i]
			coe<-names(RES[[i]]$coefficients[-1])
			equi<-match(coe,rownames(ImpValue))
			for(j in 1:length(coe)){
				ImpValue[equi[j],]=ImpValue[equi[j],]+weight
			}
		}
	}
	if(only_intercept==TRUE){
		N_models=N_models+1
	}
	Results<-list(call=match.call(),N_models=N_models,Models=RES,IC=TAB_IC,ImpValue=ImpValue)
	class(Results)<-"models"
	return(Results)
}