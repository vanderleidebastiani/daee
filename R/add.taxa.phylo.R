#' @title Add species in a phylogenetic tree
#'
#' @description Function add species in a phylogenetic tree. Three method are available, see details. This function is based in the functions \code{\link{bind.tree}} and \code{\link{bind.tip}}.
#'
#' @details  The function add.taxa.phylo.random add new tips at random. Taxa object is matrix with two columns. The first column is the node label to new tip be anchored. The tip will be anchored random inside this clade. The most recent common ancestors known can be used to define the point of anchorage. The second column is the label to new tip. Two methods are available, the method "length" add the tip with probability equivalent to edge length, long edge have more probability than short edge. The method "equal" each edge has equal probability to receive the new tip.
#' 
#' The function add.taxa.phylo allows to add the new tips according to more related species or node provided. The taxa object is matrix with three columns. The first column is the node label or tip label to new tip be anchored. The second column is the label to new tip. The third column is the terminal edge length for the added tips, can be NA. If anchored label is a node the new tip is inserted as polytomy and edge length is not used, the length is computed to keep the tree ultrametric. If anchored is a tip label and edge length was provided the new tip is inserted with length provide, and if edge length was NA the length is computed simply by split current edge length by two.
#' 
#' The function add.taxa.phylo.phylomatic allows to add the new tips according structure provided, this method is similar to software phylomatic. The taxa object is matrix with two columns. The first column is the node label or tip label to new tip be anchored. The second column is the structure to new tip. This structure is composed by labels separated by common slash (/), for example family/genus/sp. The anchored label is removed and the edge length is divided by the number of levels (number of slash) within each clade. Artificial tips (named TReMoVtIP) are inserted to mark the steps of splits, they can be removed.
#' 
#' @encoding UTF-8
#' @import ape
#' @aliases add.taxa.phylo add.taxa.phylo.random add.taxa.phylo.phylomatic
#' @importFrom phytools bind.tip add.random
#' @param tree as an object of class "phylo".
#' @param taxa Taxa structure to add species. See details.
#' @param random.order Logical argument (TRUE or FALSE) to specify if sequence order in input is keeped (FALSE) or random (TRUE) (default random.order = TRUE). 
#' @param renove.artificial.tip Logical argument (TRUE or FALSE) to specify if artificial tips are removed or not, only to add.taxa.phylo.phylomatic function (default renove.artificial.tip = TRUE). 
#' @param m One number to starting sequence (default m = 0).
#' @param prefix The prefix to new nodes labels (default "NeWNodEPhylO", "NeWNodERandoM" or "NeWNodEPhylomatiC").
#' @param method The method to add the new tip, partial match to "equal" or "length". See details (Default method = "length").
#' @return A list with: \item{call}{The call arguments.}\item{m.start}{The start m value.} \item{m.current}{The new m value, the difference between m.start and m.current is the number of node names created.} \item{prefix}{Prefix used.} \item{new.tips}{The new tips.}\item{tree}{The tree, class "phylo".} 
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords daee
#' @seealso \code{\link{bind.tree}}, \code{\link{bind.tip}}, \code{\link{add.random}}, \code{\link{add.tip.random}}
#' @examples
#' tree <- read.tree(text="(C:32,(B:16,A:16)N2:16)N1;")
#' tree <- compute.brlen(tree)
#' plot(tree, show.node=TRUE)
#' 
#' # add.taxa.phylo without edge length
#' taxa <- matrix(c("A", "C", "N2", "N2", "A1","C1", "in_N2.1", "in_N2.2",rep(NA,4)),4,3)
#' taxa
#' plot(add.taxa.phylo(tree, taxa)$tree, show.node.label = TRUE)
#' 
#' # add.taxa.phylo taxa with edge length
#' taxa <- matrix(c("A", "C", "A1","C1", 0.01, 0.9), 2,3)
#' taxa
#' plot(add.taxa.phylo(tree, taxa)$tree, show.node.label = TRUE)
#' 
#' # add.taxa.phylo.random
#' taxa <- matrix(c("N1","N2","Anywhere","Only_in_N2"), 2, 2)
#' taxa
#' plot(add.taxa.phylo.random(tree, taxa)$tree, show.node.label = TRUE)
#' plot(add.taxa.phylo.random(tree, taxa)$tree, show.node.label = TRUE)
#' 
#' # add.taxa.phylo.phylomatic
#' taxa <- matrix(c("B","B","D/D1/DD1","D/D2/DD2/DDD1"),2,2)
#' taxa
#' plot(add.taxa.phylo.phylomatic(tree, taxa, renove.artificial.tip = TRUE)$tree, 
#' 	 show.node.label = TRUE)
#' plot(add.taxa.phylo.phylomatic(tree, taxa, renove.artificial.tip = FALSE)$tree, 
#' 	 show.node.label = TRUE)
#' @export
add.taxa.phylo <- function(tree, taxa, m = 0, prefix = "NeWNodEPhylO"){
	res <- list(call = match.call())
	res$m.start <- m
	if (!inherits(tree, "phylo")){ 
		stop("\n Object tree is not of class phylo \n")
	}
	if(is.null(tree$node.label)){
		stop("\n tree$node.label is NULL. Use the function node.tree\n")
	}
	tree.labels <- c(tree$tip.label, tree$node.label)
	match.names <- match(taxa[,1], tree.labels)
	if(any(is.na(match.names))){
		print("Some nodes/tips in taxa are not present in tree$node.label/tree$tip.label:")
		mNA <- is.na(match.names)
		for(i in (1:nrow(taxa))[mNA]){
			print(as.character(taxa[i, 1]))
		}
		stop("\n Check tree and/or taxa")
	}
	show.warning <- FALSE
	for(i in 1:nrow(taxa)){
		if(length(which(c(tree.labels) == taxa[i, 1]))>1){
			stop(paste("\n The label", taxa[i,1], "with multiple occurrences in tree$node.label"))
		}
		info.temp <- tree.label.info(tree, taxa[i, 1])
		if(!is.na(taxa[i,3]) & (info.temp$type == "node" | info.temp$type == "root")){
			taxa[i,3] <- NA
			show.warning <- TRUE
		}
		if(!is.na(taxa[i,3]) & info.temp$edge.length<taxa[i,3]){
			stop(paste("\n The edge length to larger to tip", taxa[i,2]))
		}
	}
	if(show.warning){
		warning("\n Edge length not used in internal node anchor")
	}
	u.taxa <- unique(taxa[,1])
	for(i in 1:length(u.taxa)){
		if(length(unique(as.character(taxa[which(taxa[,1] == u.taxa[i]),3])))>1)
			stop("\n Edge length in a sigle anchor tip must be equal")
	}
	edges.length <- as.numeric(taxa[, 3, drop = FALSE])
	control<-matrix(NA, nrow(taxa), 2)
	control[,1] <- taxa[,1]
	for(i in 1:nrow(taxa)){
		label <- taxa[i, 1]
		if(length(which(label == tree$tip.label)) > 0){
			control.temp<-control[which(control[,1] == label), 2]
			if(all(is.na(control.temp))){
				where <- which(tree$tip.label == label)
				if(is.na(edges.length[i])){
					edges.length[i] <- tree.label.info(tree, label)[1, 2]/2 #edge.length
				}
				tree <- phytools::bind.tip(tree, taxa[i, 2], where = where, position = edges.length[i])
			} else {
				control.tips <- c(control.temp, label)
				control.tips <- control.tips[!is.na(control.tips)]
				node.name <- c(tree$tip.label, tree$node.label)[ape::getMRCA(tree, control.tips)]
				where <- which(c(tree$tip.label, tree$node.label) == node.name)
				tree <- phytools::bind.tip(tree, taxa[i, 2], where = where, position = 0)
			}
			control[i, 2] <- taxa[i, 2]
		} else {
			where <- which(c(tree$tip.label, tree$node.label) == label)
			tree <- phytools::bind.tip(tree, taxa[i, 2], where = where, position = 0)
		}
		tree.temp <- node.tree(tree, m = m, prefix = prefix)
		tree <- tree.temp$tree
		m <- tree.temp$m.current
	}
	res$m.current <- m
	res$prefix <- prefix
	res$new.tips <- taxa[,2]
	res$tree <- tree
	return(res)
}