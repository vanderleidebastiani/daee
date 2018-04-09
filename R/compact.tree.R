#' @title Compact a phylogenetic tree
#' 
#' @description Funtion to compact a phylogenetic tree. See details.
#' 
#' @details This function compacts the node of a phylogenetic tree. Inside each node the function removes tips
#' using the function \code{\link{drop.tip}}. Only two tips are keeped: the longhest and shortest edge length within each 
#' node. This can be useful only for the function \code{\link{plotcollapse.phylo}}.
#' 
#' @encoding UTF-8
#' @import ape
#' @importFrom phylobase descendants phylo4 ancestor
#' @param tree phylogeny as an object of class "phylo".
#' @param nodes a vector with node label to compact the edges. Only two tips are keeps, most short a most long within each nodes
#' @return The phylogeny as an object of class "phylo" with nodes compacted. 
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{plotcollapse.phylo}} \code{\link{graphical.node.patterns}}
#' @keywords daee
#' @examples
#' set.seed(10)
#' tree <- rtree(15)
#' tree <- makeNodeLabel(tree)
#' plot.phylo(tree, show.node.label = TRUE)
#' nodes <- c("Node2", "Node9")
#' res <- compact.tree(tree, nodes)
#' res
#' plot.phylo(res, show.node.label = TRUE)
#' @export
compact.tree <- function(tree, nodes){
  if (!inherits(tree, "phylo")) {
    stop("\n Object tree is not of class phylo \n")
  }
  tree <- node.tree(tree)$tree
  tree1 <- phylobase::phylo4(tree)
  check.nodes <- unlist(lapply(nodes, function(node) names(phylobase::descendants(node, phy = tree1, type = "all"))))
  if(any(check.nodes %in% nodes)){
    stop("nodes cannot be collapsed one inside other")
  }
  node.inf <- lapply(nodes, tree.label.info, tree = tree)
  node.inf <- as.numeric(organize.list(node.inf)[,1])
  if(any(is.na(node.inf))){
    stop("node to collapse must not be the root") 
  }
  # new.order <- order(node.inf, decreasing = TRUE)
  # nodes <- nodes[new.order]
  for(i in 1:length(nodes)){
    CCC.org <- organize.list(lapply(c(names(phylobase::descendants(tree1, nodes[i], "tips"))), tree.label.info, tree = tree))
    tip.order <- rownames(CCC.org)[(order(as.numeric(CCC.org[,3]), decreasing = TRUE))]
    if(length(tip.order)>2){
      tree <- drop.tip(tree, tip.order[c(-1, -length(tip.order))])
      tree1 <- phylobase::phylo4(tree)
      tree$node.label[which(tree$node.label  %in% names(phylobase::ancestor(tree1, tip.order[1])))] <- nodes[i]
    }
  }
  return(tree)
}