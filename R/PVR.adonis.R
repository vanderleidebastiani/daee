#' @title Phylogenetic Eigenvector Regression (PVR) and eigenvector selection
#' 
#' @description Phylogenetic Eigenvector Regression (PVR) and eigenvector selection 
#' using analysis of variance with distance matrices (adonis).
#' 
#' @details The phylogenetic distance matrix is double-centered and submitted to principal
#' coordinates analysis (PCoA). This method generates orthogonal eigenvectors that
#' summarize the phylogenetic structure (Diniz-Filho et al 2008).
#' 
#' This function is similar the function \code{\link{PVR}} that use a non-sequential approach 
#' to perform the eigenvector selection, but the selection is based in multivariate
#' analysis of variance. The function search to combination of eigenvectors that maximize 
#' the F value in the analysis of variance with distance matrices using the \code{\link{adonis}} function. 
#' Primarily, an analysis for each eigenvectors is performed, obtaining the F values. Then, the function 
#' select the eingenvector with the higher F value, and then, new eigenvectors are added to the model, 
#' models are updated and F values are obtained. The search stops when all eigenvectors are included in 
#' the model. The subset of eigenvectors that maximize the global F value must be selected manually in the results.
#' 
#' @encoding UTF-8
#' @importFrom PCPS wcmdscale.org
#' @importFrom vegan adonis
#' @param traits Data matrix or a dissimilarity matrix (recommended), usually related to species traits. 
#' This will be passed to the left side of the formula in the adonis function. The sequence of species in 
#' the traits data matrix or dissimilarity matrix must be the same as that in the phylogenetic distance 
#' matrix. See details in \code{\link{adonis}}.
#' @param dist Phylogenetic distance matrix.
#' @param cumulative Percentage of variation in the phylogenetic distances 
#' considered in the analysis. Cumulative percentage must be higher than the 
#' cumulative percentage of the first two eigenvalues, and less than 1.
#' @return \item{values}{Eigenvalues, relative eigenvalues and cumulative eigenvalues 
#' for the PCoA of distance matrix.} \item{vectors}{The principal coordinates with 
#' positive eigenvalues.} \item{inf.cumulative}{Percentage of the variation in the 
#' phylogenetic distances considered in the analysis (The result should be approximately
#' the specified cumulative value).} \item{n.axis.considered}{Number of axes considered.}
#' \item{results.unique}{F value for each PVR axis} \item{results.sequential}{F value for 
#' sequential approach using all PVR axes (PVR 1,PVR 1 + PVR 2, ...).} 
#' \item{results.stepwise}{F value for non-sequential approach, that uses the combination 
#' of PVRs axes that maximize the F value. The selection finishes using all PVRs 
#' considered, but the max F value must be selected manually in the results.}
#' @references Diniz-Filho, J. A. F., Santana, C. E. R., Bini, L. M. (1998). An 
#' eigenvector method for estimating phylogenetic inertia. Evolution, 52(5), 1247-1262.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{PVR}}
#' @keywords daee
#' @examples
#'
#' # require(SYNCSA)
#' # require(vegan)
#' # data(flona)
#' # traits.dist <- vegdist(decostand(flona$traits[,c(1,3)], 
#' # 								 method = "standardize"), 
#' # 					   method = "euclidean")
#' # results <- PVR.adonis(traits.dist, flona$phylo, cumulative = 0.7)
#' # results
#' # plot(factor(results$results.unique$PVR, levels =results$results.unique$PVR), 
#' # 	 results$results.unique$F.value, 
#' # 	 xlab = "PVR", ylab = "F value", main = "results.unique")
#' # plot(factor(results$results.sequential$PVRs, levels = results$results.sequential$PVRs), 
#' # 	 results$results.sequential$F.value,
#' # 	 xlab = "PVRs", ylab = "F value", main = "results.sequential")
#' # plot(factor(results$results.stepwise$PVRs, levels = results$results.stepwise$PVRs), 
#' # 	 results$results.stepwise$F.value,
#' # 	 xlab = "PVRs", ylab = "F value", main = "results.stepwise")
#' 
#' @export
PVR.adonis <- function(traits, dist, cumulative = 0.99){
	RES <- list(call = match.call())
	if(0>cumulative | cumulative>1){
		stop("\n Cumulative percentage must be higher than the cumulative percentage of the first two eigenvalues, and less than 1\n")
	}
	# PVRs
	ord <- PCPS::wcmdscale.org(as.dist(dist), squareroot = TRUE, eig = TRUE, correlations = FALSE)
	if(cumulative<sum(ord$values$Relative_eig[1:2])){
		print(paste("Relative eigenvalue for axis 1 =", round(ord$values$Relative_eig[1], 3)))
		print(paste("Relative eigenvalue for axis 2 =", round(ord$values$Relative_eig[2], 3)))
		stop("\n Cumulative percentage must be higher than the cumulative percentage of the first two eigenvalues\n")
	}
	used <- which(ord$values[,3]<=cumulative)
	if(cumulative==1){
		used <- used[-(length(used))]
	}
	x <- ord$vectors[,used, drop = FALSE]
	fac <- length(used)
	# Unique and sequential
	result.unique <- matrix(NA, nrow = fac, ncol = 1)
	result.seq <- matrix(NA, nrow = fac, ncol = 1)
	for(i in 1:fac){
		mod.temp.uniq <- vegan::adonis(traits ~ x[, i, drop = FALSE], permutations = 1)
		mod.temp.seq <- vegan::adonis(traits ~ x[, 1:i, drop = FALSE], permutations = 1)
		result.unique[i,] <- mod.temp.uniq$aov.tab$F.Model[1]
		result.seq[i,] <- mod.temp.seq$aov.tab$F.Model[1]
	}
	RES.unique <- cbind.data.frame(used, result.unique)
	colnames(RES.unique) <- c("PVR", "F.value")
	RES.seq <- cbind.data.frame(paste0("1 to ", used), result.seq)
	colnames(RES.seq) <- c("PVRs", "F.value")
	# Stepwise
	result.step.pvr <- matrix(NA, nrow = fac, ncol = 1)
	result.step.f <- matrix(NA, nrow = fac, ncol = 1)
	remainder <- used
	included <- c()
	for(i in 1:fac){
		result.temp <- matrix(NA, nrow = 1, ncol = length(remainder))
		for(j in 1:length(remainder)){
			include.temp<-c(included, remainder[j])
			mod.temp <- vegan::adonis(traits ~ x[, include.temp, drop = FALSE], permutations = 1)
			result.temp[1, j] <- mod.temp$aov.tab$F.Model[1]
		}
		included <- c(included, remainder[which(result.temp==max(result.temp))])
		remainder <- setdiff(used, included)
		result.step.pvr[i,] <- paste(included, sep = "", collapse = " ")
		result.step.f[i,] <- max(result.temp)
	}
	RES.step <- cbind.data.frame(result.step.pvr, result.step.f)
	rownames(RES.step) <- paste0("Step.", used)
	colnames(RES.step) <- c("PVRs", "F.value")
	# Results
	RES$values <- ord$values
	RES$vectors <- ord$vectors
	RES$n.axis.considered <- fac
	RES$inf.cumulative <- ord$values[fac, 3]
	RES$results.unique <- RES.unique
	RES$results.sequential <- RES.seq
	RES$results.stepwise <- RES.step
	return(RES)
}
