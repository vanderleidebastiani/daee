#' @title Add tips at random to the tree
#' 
#' @description Function to add new tip at random to a tree with branch length. This function is a small change in the function \code{\link{add.random}}, see details.
#' 
#' @details This function is a small change in the function \code{\link{add.random}} allowing add only one tip at random only in ultrametric phylogenetic tree. 
#' Two methods are available. The method "length" add the tip with probability equivalent to edge length, long edge have more probability than short edge. 
#' The method "equal" each edge has equal probability to receive the new tip.
#' 
#' @encoding UTF-8
#' @importFrom phytools bind.tip
#' @importFrom stats runif
#' @param tree phylogeny as an object of class "phylo".
#' @param tip tip name for the added.
#' @param method The method to add the new tip, partial match to "equal" or "length". See details (Default method = "length").
#' @return The phylogeny as an object of class "phylo" with new tip.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{add.random}} \code{\link{add.taxa.phylo}}
#' @keywords daee
#' @examples
#' set.seed(10)
#' tree <- rtree(5)
#' tree <- compute.brlen(tree)
#' tree <- makeNodeLabel(tree)
#' plot.phylo(tree, show.node.label = TRUE)
#' res <- add.tip.random(tree, "NEW", method = "length")
#' res
#' plot.phylo(res, show.node.label = TRUE)
#' @export
add.tip.random<-function (tree, tip, method = "length") {
	METHOD <- c("length", "equal")
	method <- pmatch(method, METHOD)
	if (length(method) > 1) {
		stop("\n Only one argument is accepted in method \n")
	}
	if (is.na(method)) {
		stop("\n Invalid method \n")
	}
	if(length(tip)>1){
		stop("\n Only one tip is allowed \n")
	}
	randomPosn <- function(tree, method) {
		if(method == 1){
			cum.edge <- cumsum(tree$edge.length)
			random.number <- stats::runif(1)
			pos <- random.number * sum(tree$edge.length)
			edge <- 1
			while (pos > cum.edge[edge]) edge <- edge + 1
			res <- list(node = tree$edge[edge, 2], posn = cum.edge[edge] - pos)
		} 
		if(method == 2){
			edge <- sample(seq_len(length(tree$edge.length)), 1)
			posn <- stats::runif(1, min = 0, max = tree$edge.length[edge])
			res <- list(node = tree$edge[edge, 2], posn = posn)
		}
		return(res)
	}
	where <- randomPosn(tree, method)
	tree <- phytools::bind.tip(tree, tip, where = where$node, position = where$posn)
	return(tree)
}