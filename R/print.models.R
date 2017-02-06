#' @rdname all.models
#' @encoding UTF-8
#' @export
print.models<-function(x, n = 10, ...){
	if(x$N_models<n){
		n<-x$N_models
	}
	cat("Call:\n")
	cat(deparse(x$call), "\n\n")
	cat("Number of models:\n")
	cat(deparse(x$N_models), "\n\n")
	cat("Best models:\n")
	print(as.matrix(x$IC[1:n,]))
	invisible(x)
}