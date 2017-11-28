#' @rdname add.taxa.phylo
#' @encoding UTF-8
#' @export
add.taxa.phylo.phylomatic<-function(tree, taxa, renove.artificial.tip = TRUE,  m = 0, prefix = "NeWNodEPhylomatiC"){
	res <- list(call = match.call())
	res$m.start <- m
	if (!inherits(tree, "phylo")){ 
		stop("\n Object tree is not of class phylo \n")
	}
	if(is.null(tree$node.label)){
		stop("\n tree$node.label is NULL. Use the function node.tree\n")
	}
	spp <- strsplit(cbind(paste(taxa[,1], taxa[,2], sep = "/")), "/")
	SPP <- matrix(NA, length(spp), max(sapply(spp, length)))
	for(i in 1:length(spp)){
		SPP[i, 1:length(spp[[i]])] <- spp[[i]]
	}
	tree.labels <- c(tree$tip.label, tree$node.label)
	match.names <- match(taxa[,1], tree.labels)
	if(any(is.na(match.names))){
		print("Some species in taxa are not present in the tip labels or node label:")
		mNA<-is.na(match.names)
		for(i in 1:nrow(taxa)){
			if(mNA[i]){
				print(as.character(taxa[i,1]))
			}
		}
		stop("\n Check tree and/or taxa")
	}
	for(i in 1:nrow(taxa)){
		if(length(which(c(tree.labels) == taxa[i, 1]))>1){
			stop(paste("\n The label", taxa[i,1], "with multiple occurrences in tree$tip.label or tree$node.label"))
		}
	}
	show.warning <- FALSE
	for(i in 1:nrow(taxa)){
		for(j in 2:length(spp[[i]])){
			match.names<-match(spp[[i]][j], tree.labels)
			if(!is.na(match.names)){
				show.warning <- TRUE
				spp[[i]][1] <- spp[[i]][j]
				spp[[i]][j] <- NA
			}
		}
	}
	if(show.warning){
		warning("\n structure of taxa was reduced, because some intermediaries nodes were found in the tree")
		for(i in 1:nrow(taxa)){
			if(length(which(c(tree.labels) == spp[[i]][1]))>1){
				stop(paste("\n The label", spp[[i]][1], "with multiple occurrences in tree$tip.label or tree$node.label"))
			}
		}
	}
	for(i in 1:length(spp)){
		spp[[i]] <- spp[[i]][!is.na(spp[[i]])]
	}
	SPP <- SPP[sapply(spp, length)>1, , drop = FALSE]
	spp <- spp[sapply(spp, length)>1]
	if(nrow(SPP)==0){
		stop("\n No tips to add, because all intermediaries nodes were found in the tree")
	}
	CONTROL <- as.data.frame(matrix(NA, length(spp), 7))
	colnames(CONTROL) <- c("base.node", "base.edge.length", "type", "n.split", "new.tip", "edge.length", "max.split")
	for(i in 1:length(spp)){
		info.temp <- tree.label.info(tree, spp[[i]][1])
		CONTROL[i,1]<-spp[[i]][1]
		CONTROL[i,2]<-info.temp$edge.length
		CONTROL[i,3]<-as.character(info.temp$type)
		if(CONTROL[i,3]=="root"){
			CONTROL[i,2]<-info.temp$max.height
		}
		CONTROL[i,4]<-length(spp[[i]])
		CONTROL[i,5]<-SPP[i,(CONTROL[i,4])]
		if(CONTROL[i,3]!="tip"){
			CONTROL[i,4] <- CONTROL[i,4]-1
		}
	}
	for(i in 1:length(unique(SPP[,1]))){
		CONTROL[which(SPP[, 1] == unique(SPP[,1])[i]), 6] <- CONTROL[which(SPP[, 1] == unique(SPP[, 1])[i]), 2]/max(CONTROL[which(SPP[, 1] == unique(SPP[,1])[i]), 4])
		CONTROL[which(SPP[, 1] == unique(SPP[,1])[i]), 7] <- max(CONTROL[which(SPP[, 1] == unique(SPP[,1])[i]), 4])
	}
	for(k in 2:ncol(SPP)){
		new.tips <- unique(SPP[,k])
		new.tips <- new.tips[!is.na(new.tips)]
		for(i in 1:length(new.tips)){
			label <- unique(SPP[which(SPP[,k] == new.tips[i]), k-1])
			edge.length <- CONTROL[which(SPP[, k] == new.tips[i]), 6][1]
			if(length(which(label==tree$tip.label)) > 0){
				tip <- list(edge = matrix(c(3, 3, 1, 2), 2, 2), tip.label = c(new.tips[i], "TReMoVtIP"), Nnode = 1, node.label = label, edge.length = c(edge.length, edge.length))
				if(any(new.tips[i] == CONTROL[, 5])){
					ntemp <- CONTROL[which(SPP[, k] == new.tips[i]), 7]
					if((CONTROL[which(SPP[,k] == new.tips[i]), 4] < ntemp)){
						NN <- ntemp - CONTROL[which(SPP[, k] == new.tips[i]), 4]
						tip$edge.length <- c((edge.length*(NN+1)), (edge.length*(NN+1)))
					}
				}
				class(tip) <- "phylo"
				tree$edge.length[tree.label.info(tree, label)[1, 1]] <- edge.length
				where <- which(c(tree$tip.label, tree$node.label) == label)
				tree <- ape::bind.tree(tree, tip, where = where, position = 0)
			} else {
				where <- which(c(tree$tip.label, tree$node.label)==label)
				if(any(new.tips[i] == CONTROL[, 5])){
					ntemp <- CONTROL[which(SPP[, k] == new.tips[i]), 7]
					if((CONTROL[which(SPP[,k] == new.tips[i]), 4] < ntemp)){
						NN <- ntemp-CONTROL[which(SPP[, k] == new.tips[i]), 4]
						edge.length <- edge.length*(NN+1)
					}
				}
				tree<-phytools::bind.tip(tree,new.tips[i], edge.length = edge.length, where = where, posicao = 0)
			}
			tree.temp <- node.tree(tree, m = m, prefix = prefix) 
			tree <- tree.temp$tree
			m <- tree.temp$m.current
		}
	}
	if(renove.artificial.tip){
		tree <- drop.tip(tree, "TReMoVtIP")
	}
	res$m.current <- m
	res$prefix <- prefix
	res$new.tips <- CONTROL[, 5]
	res$tree <- tree
	return(res)
}