#' @title Show information about a label in a phylogenetic tree
#'
#' @description Function for discover the edge position, edge length, edge height, max height of tree and type of the
#' label in a phylogenetic tree. This function is based in the function \code{\link{nodeHeights}}.
#'
#' @encoding UTF-8
#' @import ape
#' @importFrom phytools nodeHeights
#' @param tree phylogeny as an object of class "phylo".
#' @param label One tip label or one node label.
#' @return A data.frame with type of label, edge position, edge length, edge height and max height of tree.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{nodeHeights}}
#' @keywords daee
#' @examples
#' tree <- read.tree(text = "(C:32,(B:16,A:16)N2:16)N1;")
#' plot(tree, show.node = TRUE)
#' tree.label.info(tree, label = "A")
#' @export
tree.label.info <- function(tree, label){
	if(is.null(tree$edge.length) | any(is.na(tree$edge.length))){
		stop("\n tree$edge.length with any NA or NULL\n")
	}
	if(is.null(tree$tip.label) | is.null(tree$node.label)){
		stop("\n tree$tip.label and/or tree$node.label NULL\n")
	}
	n <- length(tree$tip.label)
	where <- which(c(tree$tip.label, tree$node.label) == label)
	if(length(where)>1){
		stop("\n Only one label accepted or label with multiple occurrences in tree\n")
	}
	if(length(where)!=1){
		stop("\n label not found\n")
	}
	H <- phytools::nodeHeights(tree)
	HM <- max(H)
	if(where<=n){
		edge <- which(tree$edge[,2] == where)
		edge.length <- tree$edge.length[edge]
		edge.height <- HM - H[edge, 2]
		type <- "tip"
	}else{
		if (where == (length(tree$tip.label) + 1)){
			edge<-which(tree$edge[,1] == where)[1]
			edge.length <- 0
			edge.height <- HM - H[edge, 1]
			edge <- NA
			type <- "root"
		} else {
			edge <- which(tree$edge[,2] == where)
			edge.length <-tree$edge.length[edge]
			edge.height <- HM - H[edge, 2]	
			type <- "node"
		}
	}
	res <- data.frame(edge, edge.length, edge.height, max.height = HM, type = type)
	rownames(res) <- label
	return(res)
}