#' @title Makes node labels
#'
#' @description This function makes node labels in a tree keeping already existing node names. This function is based in the function \code{\link{makeNodeLabel}}.
#'
#' @encoding UTF-8
#' @import ape
#' @param tree phylogeny as an object of class "phylo".
#' @param m One number to starting sequence (default m = 0).
#' @param prefix The prefix (default prefix = "NodE")
#' @return A list with: \item{call}{The call arguments.}\item{m.start}{The start m value.} \item{m.current}{The new m value, the difference between m.start and m.current is the number of node names created.} \item{prefix}{Prefix used.} \item{tree}{The tree, class "phylo".} 
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{makeNodeLabel}}
#' @keywords daee
#' @examples
#' tree <- read.tree(text = "(C:32,(B:16,A:16):16);")
#' plot(tree, show.node = TRUE)
#' tree <- node.tree(tree, m = 0, prefix = "N")$tree
#' plot(tree, show.node = TRUE)
#' @export
node.tree<-function(tree, m = 0, prefix = "NodE"){
	res <- list(call = match.call())
	res$m.start <- m
	if(is.null(tree$node.label)){
		tree <- ape::makeNodeLabel(tree)
		tree$node.label <- rep("", tree$Nnode)
	}
	for(i in 1:tree$Nnode){
		if(tree$node.label[i] == "NA" | tree$node.label[i] == ""){
			m <- m+1
			tree$node.label[i] <- paste(prefix, m, sep = "")
		}
	}
	res$m.current <- m
	res$prefix <- prefix
	res$tree <- tree
	return(res)
}