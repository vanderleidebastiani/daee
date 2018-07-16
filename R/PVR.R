#' @title Phylogenetic Eigenvector Regression
#' 
#' @description Phylogenetic Eigenvector Regression (PVR) and eigenvector selection. 
#' 
#' @details This function is based on a non-sequential approach, that uses the combination
#' of eigenvectors that minimizes the residual phylogenetic autocorrelation, 
#' measured by Moran I. The method can be used to measure the level of 
#' phylogenetic signal in ecological data and to study correlated evolution
#' (Diniz-Filho et al 2011).
#'
#' The phylogenetic distance matrix is double-centered and submitted to principal
#' coordinates analysis (PCoA). This method generates orthogonal eigenvectors that
#' summarize the phylogenetic structure (Diniz-Filho et al 2008).
#'
#' The sets of eigenvectors is selected with multiple regression model:
#' 
#' Y = a+Xb+e
#' 
#' where Y is a vector describing trait variation in the set of species, X contains
#' a set of k eigenvectors, a the intercept, b is the vector with regression coefficients estimated 
#' and e its residuals, the part of variation in Y that is not explained by the 
#' X (Diniz-Filho et al 1998). 
#'
#' The function use an iterative search for the eigenvector that reduces the autocorrelation
#' in the residuals. Primarily, the regression for all eigenvectors is calculated, obtaining
#' the residuas. Then, Moran I for each eigenvector is calculated for the residuals. The 
#' function select the eingenvector with the lowest Moran I, and then, as new eigenvectors 
#' are added to the model, residuals are updated and autocorrelation is reestimated. The search
#' stops when residual autocorrelation is reduced below threshold Moran I specified and when 
#' the statistical significance is reached (Diniz-Filho et al 2011).
#' 
#' @encoding UTF-8
#' @import ape
#' @importFrom stats lm residuals summary.lm as.dist as.formula pf
#' @importFrom graphics abline plot
#' @aliases PVR plot.pvr
#' @param traits Species described by continuous traits, with traits as columns 
#' and species as rows.
#' @param dist Phylogenetic distance matrix.
#' @param cumulative Percentage of variation in the phylogenetic distances 
#' considered in the analysis. Cumulative percentage must be higher than the 
#' cumulative percentage of the first two eigenvalues, and less than 1.
#' @param VMoran Stopping rule based on Moran I value (Absolute value, smaller than
#' the specified value).
#' @param pMoran Stopping rule based on the p-value of Moran I (Greater than the specified). 
#' check Logical argument (TRUE or FALSE) that checks  whether traits and phylogeny
#' taxa labels match. The sequence of species in the trait data must be the same as 
#' that in the phylogenetic distance matrix.
#' @param check Logical argument (TRUE or FALSE) to check if species sequence in the traits data
#' follows the same order as the one in the phylodist matrices (Default checkdata = TRUE).
#' @param x An object of class pvr.
#' @param trait Trait for plot.
#' @param ... Other parameters for the respective functions.
#' @return \item{values}{Eigenvalues, relative eigenvalues and cumulative eigenvalues 
#' for the PCoA of distance matrix.} \item{vectors}{The principal coordinates with 
#' positive eigenvalues.} \item{inf.cumulative}{Percentage of the variation in the 
#' phylogenetic distances considered in the analysis ( The result should be approximately
#' the specified cumulative value).} \item{n.axis.considered}{Number of axes considered.}
#' \item{moran.less.than}{Morans I value considered in the stopping rules (Absolute 
#' value).} \item{p.moran.greater.than}{Stopping rule for the p-value.} \item{PSR.curve.axis.x}{Values
#' for the PSR curve (Abscissa).} \item{PSR.curve.axis.y}{Values for the PSR curve
#' (Ordinate, for each traits).} \item{minimun.moran}{Parameters, number of parameters, 
#' observed Moran I, p-value for Moran I, R Squared and p-value for regression model
#' that minimize autocorrelation coefficients in the residuals for each trait.} 
#' @note The parameter in minimun.moran is shown as follows y ~  x[,3] + x[,5], 
#' in other words Trait~Axis_3+Axis_5 and so on.
#' @references Diniz-Filho, J. A. F., Santana, C. E. R., Bini, L. M. (1998). An 
#' eigenvector method for estimating phylogenetic inertia. Evolution, 52(5), 1247-1262.
#' @references Diniz-Filho, J.A.F., Bini, L. M., Rangel, T. F., Morales-Castilla, I., 
#' Olalla-Tarraga, M. A., Rodriguez, M. A. & Hawkins, B. A. (2011). On the selection of 
#' phylogenetic eigenvectors for ecological analyses. Ecography, 35(3), 239-249.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords daee
#' @examples
#'
#' data(flona)
#' Res<-PVR(flona$traits[,1:4], flona$phylo, VMoran = 0.01)
#' Res
#' plot(Res, trait = 1)
#' 
#' @export
PVR<-function(traits, dist, cumulative = 0.99, VMoran = 0.025, pMoran = 0.05, check = TRUE){
	if(0>cumulative | cumulative>1){
		stop("\n Cumulative percentage must be higher than the cumulative percentage of the first two eigenvalues, and less than 1\n")
	}
	traits<-as.matrix(traits)
	if (check){
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
	dist<-as.dist(dist)
	ord<-PCPS::wcmdscale.org(dist, squareroot = TRUE, eig = TRUE, correlations = FALSE)
	use<-which(ord$values[,3]<=cumulative)
	if(cumulative==1){
		use<-use[-(length(use))]
	}
	x<-ord$vectors[,use]
	vare<-ncol(traits)
	fac<-length(use)
	xnam <- paste("x[,", 1:fac,"]", sep="")
	result<-matrix(NA, nrow = vare, ncol = fac)
	for(i in 1:vare){
		y<-traits[,i]
		for (j in 1:fac){
			result[i,j]<-as.numeric(stats::summary.lm(stats::lm(stats::as.formula(paste("y ~ ", paste(xnam[1:j], collapse= "+")))))$r.squared)
		}	
	}
	colnames(result)<-colnames(ord$vectors)[1:fac]
	rownames(result)<-colnames(traits, do.NULL = FALSE, prefix = "Trait.")
	dist<-as.matrix(dist)
	resFIN<-matrix(NA, vare, 6)
		for(i in 1:vare){
			y<-traits[,i]
			MOR<-NULL
			xnam2<-xnam
			possi<-1:fac
			possi2<-1:fac
			f<-character()		
			res.MORAN<-matrix(NA, fac, 2)
				repeat{
					result.moran<-matrix(NA, 2, fac)
					for (j in possi2){
						regre<-stats::lm(stats::as.formula(paste("y ~ ",paste(f, sep = "", collapse = " + "), paste("+", xnam[j]), collapse = " + ")))
						moran<-ape::Moran.I(regre$residuals, dist, scaled = TRUE)
						result.moran[1,j]<-moran$observed
						result.moran[2,j]<-moran$p.value
					}
					colnames(result.moran)<-possi
					menor<-as.numeric(noquote(names(sort(result.moran[1,], decreasing = FALSE))[1]))
					xnam2<-xnam2[-which(xnam2==xnam[menor])]
					MOR<-result.moran[1,menor]
					pMOR<-result.moran[2,menor]
					f<-c(f,paste(xnam[menor]))
					res.MORAN[length(f),1]<-MOR
					res.MORAN[length(f),2]<-pMOR
					possi2<-possi2[-which(possi2==possi[menor])]			
					if (length(f)==fac) break
					if (MOR<=VMoran){
						if (pMOR>=pMoran) break
					}
				}
			minMORAN<-length(f)
			regre2<-stats::lm(stats::as.formula(paste("y ~ ", paste(f[1:minMORAN], sep = "", collapse = " + "))))
			resFIN[i,1]<-paste("y ~ ", paste(f[1:minMORAN], sep = "", collapse = " + "))
			resFIN[i,2]<-minMORAN
			moran2<-ape::Moran.I(regre2$residuals, dist, scaled = TRUE)
			resFIN[i,3]<-moran2$observed
			resFIN[i,4]<-moran2$p.value
			summary.regre2<-stats::summary.lm(regre2)
			resFIN[i,5]<-as.numeric(summary.regre2$r.squared)
			resFIN[i,6]<-as.numeric(stats::pf(summary.regre2$fstatistic[1], summary.regre2$fstatistic[2], summary.regre2$fstatistic[3], lower.tail = FALSE))
		}
	colnames(resFIN)<-c("Parameters", "N.Parameters", "Moran.Obs", "Moran.p", "R.Squared", "p")
	rownames(resFIN)<-colnames(traits, do.NULL = FALSE, prefix = "Trait.")
	res.x<-t(as.matrix(ord$values[1:fac,3]))
	colnames(res.x)<-colnames(res.x, do.NULL = FALSE, prefix = "Axis.")
	rownames(res.x)<-"Cumulative"
	inf<-res.x[1,fac]
	RES<-list(values = ord$values, vectors = ord$vectors, inf.cumulative = inf, n.axis.considered = fac, moran.less.than = VMoran, p.moran.greater.than = pMoran, PSR.curve.axis.x = res.x, PSR.curve.axis.y = result, minimum.moran = resFIN)
	class(RES)<-"pvr"
	return(RES)
}