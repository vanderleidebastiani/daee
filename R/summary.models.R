summary.models<-function(object, ...){
	res<-list()
    res$call <- object$call
    res$N_models<-object$N_models
    res$IC<-object$IC
    res$ImpValue<-object$ImpValue
    return(res)
}