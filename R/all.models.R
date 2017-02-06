#' Linear Models with all possible combinations of all variables
#' 
#' The function generates Linear Models (LM) with all possible combinations of all
#' variables included in the full model. Each model is a GLM of family gaussian with
#' the response variable modeled by one or more independent variables (Y= X1b1+e, 
#' Y= X1b1+X2b2+e, and so on).
#'
#' @encoding UTF-8
#' @import PCPS
#' @importFrom stats glm as.formula
#' @importFrom bbmle ICtab
#' @importFrom utils combn
#' @aliases all.models print.models summary.models
#' @param response A data frame containing the response variable.
#' @param variables A data frame containing all the independent variables for the models.
#' @param subset Maximum number of independent  variables to be considered in each model.
#' @param type Information criterion to be used "AIC","BIC","AICc","qAIC" and "qAICc"
#' (Default type = "AICc").
#' @param only_intercept Logical argument (TRUE or FALSE) to specify if a model containing 
#' only the intercept should be included (Default only_intercept = FALSE).
#' @param importance Logical argument (TRUE or FALSE) to specify if the relative importance
#' of variables should be calculated. (Default importance = TRUE).
#' @param object An object of class models.
#' @param x An object of class models.
#' @param n Number of model for print.
#' @param ... Other parameters for the respective functions. In all.models function the 
#' parameters are for the glm function.
#' @return \item{Envir_class}{The class of each variable.} \item{N_models}{The number of 
#' models.} \item{IC}{A data frame containing the information criterion, the number of 
#' parameters, difference in IC from minimum - IC model and the weights IC. The same that
#' function ICtab.} \item{Models}{A list with all models. Each model is the class glm.}
#' \item{ImpValues}{Relative importance of each variable.} 
#' @note If the model with only the intercept is include, this will be the last model.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{glm}} \code{\link{ICtab}}
#' @keywords daee
#' @examples
#'
#' #data(mite)
#' #data(mite.pcnm)
#' #response<-rowSums(mite)
#' #Res<-all.models(response,mite.pcnm,subset=3)
#' #Res
#' #summary(Res)
#' #Res$Models$Model_1157
#' 
#' @export
all.models<-function(response, variables, subset, type = "AICc", only_intercept = FALSE, importance = TRUE, ...){
	ImpValue<-NULL
	colnames(variables)<-colnames(variables, do.NULL = FALSE, prefix = "var_")
	dados<-data.frame(response, variables)
	colnames(dados)=c("y", colnames(variables))
	envir_class <- matrix(NA, dim(dados)[2], 1)
	rownames(envir_class) <- colnames(dados)
	colnames(envir_class) <- c("Class")
	for (j in 1:dim(dados)[2]) {
		envir_class[j,1]<-class(dados[,j])
	}
	s<-1:subset
	m<-dim(variables)[2]
	nobs<-dim(variables)[1]
	possibilities<-factorial(m)/(factorial(s) * factorial(m - s))
	N_models<-sum(possibilities)
	RES <- vector("list", N_models)
	for(i in 1:subset){
		combinations <- utils::combn(colnames(variables), i, simplify = TRUE)
		for (j in 1:possibilities[i]) {
			RES[[(j + sum(possibilities[1:i - 1]))]] <- stats::glm(stats::as.formula(paste("y ~ ", paste(combinations[, j], collapse = "+"))),data = dados, ...)
		}
	}
	if(only_intercept==TRUE){
		RES[[N_models+1]]<-stats::glm(y~1, data = dados, ...)
		names(RES)=sprintf("Model_%.d",1:(N_models+1))
	}
	else{
		names(RES)=sprintf("Model_%.d",1:N_models)
	}
	TAB_IC<-bbmle::ICtab(RES, type = type, nobs = nobs, base = TRUE, weights = TRUE, mnames = names(RES))
	class(TAB_IC)="data.frame"
	if(importance){
		coef_uni<-colnames(variables)
		ImpValue<-matrix(0,length(coef_uni),1)
		rownames(ImpValue)<-coef_uni
		colnames(ImpValue)<-"Relative_Importance"
		tab_ic<-bbmle::ICtab(RES, type = type, nobs = nobs, weights = TRUE, sort = FALSE, mnames = names(RES))
		for(i in 1:N_models){
			weight<-tab_ic$weight[i]
			coe<-names(RES[[i]]$coefficients[-1])
			equi<-match(coe,rownames(ImpValue))
			for(j in 1:length(coe)){
				ImpValue[equi[j],]<-ImpValue[equi[j],]+weight
			}
		}
	}
	if(only_intercept){
		N_models<-N_models+1
	}
	Results<-list(call = match.call(), Envir_class = envir_class, N_models = N_models, Models = RES, IC = TAB_IC, ImpValue = ImpValue)
	class(Results)<-"models"
	return(Results)
}