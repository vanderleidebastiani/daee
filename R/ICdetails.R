#' @title Extract details of information criterion 
#' 
#' @description The function extract details of information criterion (IC, e.g. AIC, AICc) allowing 
#' calculate difference in IC from minimum - IC model (delta), the Akaike weights and the sum of 
#' Akaike weights by each explanatory variable (relative importance of each variable).
#' 
#' @encoding UTF-8
#' @aliases ICdetails
#' @param x A data frame or matrix containing the information criterion for each model. The 
#' information criterion in the first column and the models in rows. See examples.
#' @param variables A data frame or matrix containing the indication of all independent 
#' variables used in each model. The models in row and variables in columns. This must be 
#' binary, where 1 indicate presence of variable in respective model and 0 the ausence of 
#' the variable  in respective model. See examples (default variables = NULL).
#' @param order Logical argument (TRUE or FALSE) to specify if informtion criterion is sorted
#' (Default order = TRUE).
#' @return A list with: \item{call}{The arguments used.} \item{IC}{A data frame containing the 
#' information criterion, difference in IC from minimum - IC model (delta) and the weights IC.} 
#' \item{IValue}{Relative importance of each variable.} 
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @references 
#' Anderson, D.R. 2008. Model Based Inference in the Life Sciences: A Primer on Evidence. Springer-Verlag New York. 
#' @seealso \code{\link{ICtab}}
#' @keywords daee
#' @examples
#'
#' my.aic <- matrix(c(4,2,3,5),4,1)
#' colnames(my.aic) <- "AIC"
#' rownames(my.aic) <- rownames(my.aic, do.NULL = FALSE, prefix = "M.")
#' my.aic
#' ICdetails(my.aic)
#' 
#' my.models.by.variables<-matrix(c(1,0,0,0,0,1,1,0,0,1,0,1,0,0,0,1),4,4)
#' colnames(my.models.by.variables) <- colnames(my.models.by.variables, 
#'     do.NULL = FALSE, prefix = "Var.")
#' rownames(my.models.by.variables) <- rownames(my.aic)
#' my.models.by.variables
#' ICdetails(my.aic, variables = my.models.by.variables)
#' 
#' @export
ICdetails <- function(x, variables = NULL, order = TRUE){
	if(ncol(x)>1){
		stop("x must be only one column")
	}
	if(is.null(colnames(x))){
		colnames(x) <- "IC"
	}
	rownames(x) <- rownames(x, do.NULL = FALSE, prefix = "model.")
	x <- data.frame(x)
	res <- list(call =  match.call())
	x$delta <- x[,1]-min(x[,1])
	x$weights <- round(exp((-(1/2))*x$delta)/sum(exp((-(1/2))*x$delta)), 3)
	if(!is.null(variables)){
		if(nrow(x) != nrow(variables)){
			stop("number of rows in variables not equal number of rows in x")
		}
		if(!all(variables %in% c(0, 1))){
			stop("variables must be a binary matrix")
		}
		variables <- as.matrix(variables)
		colnames(variables) <- colnames(variables, do.NULL = FALSE, prefix = "v.")
		weight <- as.matrix(x$weight)
		colnames(weight) <- "relative.importance"
		res$IValue <- t(variables)%*%weight
	}
	if(order){
		x <- x[order(x[,1]),]
	}
	res$IC = x
	return(res)
}