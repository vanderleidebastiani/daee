#' @rdname allmodels
#' @encoding UTF-8
#' @export
summary.allmodels<-function(object, ...){
	res<-list()
    res$call <- object$call
    res$N_models<-object$N_models
    res$Envir_class<-object$Envir_class
    res$IC<-object$IC
    res$ImpValue<-object$ImpValue
    return(res)
}