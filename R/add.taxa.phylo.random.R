#' @rdname add.taxa.phylo
#' @encoding UTF-8
#' @export
add.taxa.phylo.random <- function(tree, taxa, random.order = TRUE, m = 0, prefix = "NeWNodERandoM", method = "length"){
	# method = c("m1", "m2")
	res <- list(call = match.call())
	res$m.start <- m
	if (!inherits(tree, "phylo")){ 
		stop("\n Object tree is not of class phylo \n")
	}
	# METHOD <- c("m1", "m2")
	# method <- pmatch(method, METHOD)
	# if (length(method) > 1) {
	# 	stop("\n Only one argument is accepted in method \n")
	# }
	# if (is.na(method)) {
	# 	stop("\n Invalid method \n")
	# }
	if(random.order){
		taxa <- taxa[sample(1:nrow(taxa)), ,drop = FALSE]
	}
	if(is.null(tree$node.label)){
		stop("\n tree$node.label is NULL. Use the function node.tree\n")
	}
	if(method==1){
		tree.labels <- c(tree$node.label) # only node label
	} else{
		tree.labels <- c(tree$node.label, tree$tip.label)
	}
	match.names <- match(taxa[,1], tree.labels)
	if(any(is.na(match.names))){
		print("Some nodes in taxa are not present in labels:")
		mNA <- is.na(match.names)
		for(i in (1:nrow(taxa))[mNA]){
			print(as.character(taxa[i, 1]))
		}
		stop("\n Check tree and/or taxa")
	}
	for(i in 1:nrow(taxa)){
		if(length(which(c(tree.labels) == taxa[i, 1]))>1){
			stop(paste("\n The label", taxa[i,1], "with multiple occurrences in labels"))
		}
	}
	for(i in 1:nrow(taxa)){
		edge.info <- tree.label.info(tree, taxa[i, 1])
		if(edge.info$type=="root"){
			# tree <- phytools::add.random(tree, tips = taxa[i, 2])
			tree <- add.tip.random(tree, tip = taxa[i, 2], method = method)
		}
		if(edge.info$type=="node"){
			tree.clade <- ape::extract.clade(tree, node = taxa[i, 1])
			# if(method==1){
			# 	tree.random <- phytools::add.random(tree.clade, tips = taxa[i, 2])
			# 	tree <- ape::drop.tip(tree, tree.clade$tip.label, trim.internal = FALSE, subtree = TRUE)
			# 	tree$tip.label[which(tree$tip.label == paste("[", length(tree.clade$tip.label), "_tips]", sep = ""))] <- taxa[i, 1]
			# 	plot.phylo(tree, show.node.label = T)
			# 	where <- which(c(tree$tip.label, tree$node.label) == taxa[i, 1])
			# 	tree <- ape::bind.tree(tree, tree.random, where = where)
			# } else{
			tree.clade.tips <- tree.clade$tip.label
			tip <- list(edge = matrix(c(3, 3, 1, 2), 2, 2), tip.label = c("TReMoVtIP1TemP", "TReMoVtIP2TemP"), Nnode = 1, node.label = "NTReMoV1TemP", edge.length = c(edge.info$edge.length, edge.info$edge.length+edge.info$edge.height))
			class(tip) <- "phylo"
			tree.clade <- ape::bind.tree(tip, tree.clade, where = 1, position = 0)
			go <- TRUE
			while(go){
				# tree.random <- phytools::add.random(tree.clade, tips = taxa[i, 2])
				tree.random <- add.tip.random(tree.clade, tip = taxa[i, 2], method = method)
				node.check <- getMRCA(tree.random, c("TReMoVtIP2TemP", taxa[i, 2]))
				if(c(tree.random$tip.label, tree.random$node.label)[node.check]!="NA"){
					go <- FALSE
				}
			}
			node.random1 <- ape::getMRCA(tree.random, c(taxa[i, 2], tree.clade.tips))
			node.random2 <- ape::getMRCA(tree.random, c(tree.clade.tips))
			if(node.random1!=node.random2){
				tree.random$node.label[which(tree.random$node.label==taxa[i, 1])] <- "NTReMoV2TemP"
				tree.random$node.label[which(tree.random$node.label=="NA")] <- taxa[i, 1]
				tree.random$node.label[which(tree.random$node.label=="NTReMoV2TemP")] <- "NA"
			}
			tree <- ape::drop.tip(tree, tree.clade$tip.label, trim.internal = FALSE, subtree = TRUE)
			tree$edge.length[edge.info$edge] <- 0
			tree$tip.label[which(tree$tip.label == paste("[", length(tree.clade$tip.label)-1, "_tips]", sep = ""))] <- taxa[i, 1]
			where <- which(c(tree$tip.label, tree$node.label) == taxa[i, 1])
			tree <- ape::bind.tree(tree, tree.random, where = where)
			tree <- ape::drop.tip(tree, "TReMoVtIP2TemP")
			# }
		}
		if(edge.info$type=="tip"){
			tip <- list(edge = matrix(c(3, 3, 1, 2), 2, 2), tip.label = c(taxa[i, 1], "TReMoVtIP2TemP"), Nnode = 1, node.label = "XX", edge.length = c(edge.info$edge.length, edge.info$edge.length))
			class(tip) <- "phylo"
			go <- TRUE
			while(go){
				# tree.random <- phytools::add.random(tip, tips = taxa[i, 2])
				tree.random <- add.tip.random(tip, tip = taxa[i, 2], method = method)
				node.check <- getMRCA(tree.random, c("TReMoVtIP2TemP", taxa[i, 2]))
				if(c(tree.random$tip.label, tree.random$node.label)[node.check]!="NA"){
					go <- FALSE
				}
			}
			tree <- ape::drop.tip(tree, taxa[i, 1], trim.internal = T, subtree = TRUE)
			tree$edge.length[edge.info$edge] <- 0
			where <- which(c(tree$tip.label, tree$node.label) == paste("[", 1, "_tip]", sep = ""))
			tree <- ape::bind.tree(tree, tree.random, where = where)
			tree <- ape::drop.tip(tree, "TReMoVtIP2TemP")
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