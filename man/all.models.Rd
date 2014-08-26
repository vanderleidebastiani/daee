\encoding{UTF-8}
\name{all.models}
\alias{all.models}
\alias{print.models}
\alias{summary.models}

\title{Linear Models with all possible combinations of all variables}

\description{

The function generates Linear Models (LM) with all possible combinations of all variables included in the full model. Each model is a GLM of family gaussian with the response variable modeled by one or more independent variables (Y= X1b1+e, Y= X1b1+X2b2+e, and so on).

}

\usage{
all.models(response, variables, subset, type = "AICc", only_intercept = FALSE, importance = TRUE)
\method{print}{models}(x, ...)
\method{summary}{models}(object, ...)
}

\arguments{
\item{response}{A data frame containing the response variable.}
\item{variables}{A data frame containing all the independent variables for the models.}
\item{subset}{Maximum number of independent  variables to be considered in each model.}
\item{type}{Information criterion to be  used "AIC","BIC","AICc","qAIC" and "qAICc" (Default type="AICc").}
\item{only_intercept}{Logical argument (TRUE or FALSE) to specify if a model containing only the intercept should be included (Default only_intercept = FALSE).}
\item{importance}{Logical argument (TRUE or FALSE) to specify if the relative importance of variables should be calculated. (Default importance = TRUE).}
\item{object, x}{An object of class models.}
\item{...}{Other parameters for the respective functions.}
}

\value{

\item{N_models}{The number of models.}
\item{IC}{A data frame containing the information criterion, the number of parameters, difference in IC from minimum - IC model and the weights IC. The same that function ICtab.}
\item{Models}{A list with all models. Each model is the class glm.}
\item{ImpValues}{Relative importance of each variable.}
}

\note{
If the model with only the intercept is include, this will be the last model.
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{glm}}, \code{\link{ICtab}} 
}

\examples{
data(mite)
data(mite.pcnm)
response<-rowSums(mite)
Res<-all.models(response,mite.pcnm,subset=3)
Res
summary(Res)
Res$Models$Model_1157
}

\keyword{daee}
