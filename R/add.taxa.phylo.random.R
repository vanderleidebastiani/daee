#' @rdname add.taxa.phylo
#' @encoding UTF-8
#' @export
add.taxa.phylo.random <- function(tree, taxa, random.order = TRUE, m = 0, prefix = "NeWNodERandoM"){
	res <- list(call = match.call())
	res$m.start <- m
	if (!inherits(tree, "phylo")){ 
		stop("\n Object tree is not of class phylo \n")
	}
	if(random.order){
		taxa <- taxa[sample(1:nrow(taxa)), ,drop = FALSE]
	}
	if(is.null(tree$node.label)){
		stop("\n tree$node.label is NULL. Use the function node.tree\n")
	}
	tree.labels <- c(tree$node.label)
	match.names <- match(taxa[,1], tree.labels)
	if(any(is.na(match.names))){
		print("Some nodes in taxa are not present in tree$node.label:")
		mNA <- is.na(match.names)
		for(i in (1:nrow(taxa))[mNA]){
			print(as.character(taxa[i, 1]))
		}
		stop("\n Check tree and/or taxa")
	}
	for(i in 1:nrow(taxa)){
		if(length(which(c(tree.labels) == taxa[i, 1]))>1){
			stop(paste("\n The label", taxa[i,1], "with multiple occurrences in tree$node.label"))
		}
	}
	for(i in 1:nrow(taxa)){
		tree.clade <- ape::extract.clade(tree, node = taxa[i, 1])
		tree.random <- phytools::add.random(tree.clade, tips = taxa[i, 2])
		if(!length(tree$tip.label) == length(tree.clade$tip.label)){
			tree <- ape::drop.tip(tree, tree.clade$tip.label, trim.internal = FALSE, subtree = TRUE)
			tree$tip.label[which(tree$tip.label == paste("[", length(tree.clade$tip.label), "_tips]", sep = ""))] <- taxa[i, 1]
			where <- which(c(tree$tip.label, tree$node.label) == taxa[i, 1])
			tree <- ape::bind.tree(tree, tree.random, where = where)
		} else{
			tree <- tree.random
		}
		tree.temp <- node.tree(tree, m, prefix = prefix)
		tree <- tree.temp$tree
		m <- tree.temp$m.current
	}
	res$m.current <- m
	res$prefix <- prefix
	res$new.tips <- taxa[,2]
	res$tree <- tree
	return(res)
}