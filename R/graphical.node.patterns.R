#' @title Define graphical parameters to plot phylogenetic tree.
#' 
#' @description Function to define graphical parameters to each edge when draw a phylogenetic tree. See details.
#' 
#' @details This function can be used to especify diferent graphical parameters (e.g. color, width and 
#' line types) for specific nodes when draw a phylogenetic tree. First, the basicpattern argument is 
#' defined for all edges of phylogenetic tree and when the basic pattern is changed in all edge within 
#' of each node, following the nodespatterns specify. The argument force.order specify if changes following 
#' order in nodes arguments step by step (this case, some changes may have no effect) or change are done 
#' from root to tips.
#' 
#' @encoding UTF-8
#' @import ape
#' @importFrom phylobase descendants phylo4
#' @param tree phylogeny as an object of class "phylo".
#' @param nodes a vector with node label to search the edges to change the basic graphical parameter.
#' @param basicpattern the basic pattern for graphical parameter. This is apply for all edges.
#' @param nodespatterns a vector with new graphical parameter for each node label. This change the basic graphical parameter in each node.
#' @param include.node logical argument (TRUE or FALSE) to specify if edge of each node is include in change (default include.node = TRUE).
#' @param force.order logical argument (TRUE or FALSE) to specify if force the search as according to edges (default force.order = TRUE).
#' @return A vector with the new graphical parameters for each edge in phylogenetic tree.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{plot.phylo}} \code{\link{plotcollapse.phylo}}
#' @keywords daee
#' @examples
#' set.seed(10)
#' tree <- rtree(15)
#' tree <- makeNodeLabel(tree)
#' plot.phylo(tree, show.node.label = TRUE)
#'
#' edge.col <- graphical.node.patterns(tree, nodes = c("Node2", "Node8"),
#'    basicpattern = "black", nodespatterns = c("red","blue"))
#' edge.col # Color vector for each edge
#' plot.phylo(tree, show.node.label = TRUE, edge.color = edge.col)
#'
#' edge.width <- graphical.node.patterns(tree, nodes = c("Node11","Node3"), 
#'    basicpattern = 1, nodespatterns = 5, include.node = FALSE)
#' edge.width # width vector for each edge
#' plot.phylo(tree, show.node.label = TRUE, edge.width = edge.width)
#' 
#' tree <- rtree(250)
#' tree <- makeNodeLabel(tree)
#' plot(tree, show.tip.label = FALSE)
#' edge.col <- graphical.node.patterns(tree, nodes = tree$node.label, 
#'    basicpattern = "black", nodespatterns = rainbow(length(tree$node.label)))
#' edge.col
#' plot.phylo(tree, edge.color = edge.col, show.tip.label = FALSE)
#' 
#' @export
graphical.node.patterns <- function(tree, nodes, basicpattern, nodespatterns, include.node = TRUE, force.order = TRUE){
  if (!inherits(tree, "phylo")) {
    stop("\n Object tree is not of class phylo \n")
  }
  tree <- node.tree(tree)$tree
  if(length(nodes)!=length(nodespatterns)){
    nodespatterns <- rep(nodespatterns, length(nodes))
  }
  tree1 <- phylobase::phylo4(tree)
  res <- rep(basicpattern, nrow(tree$edge))
  if(force.order){
    node.inf <- lapply(nodes, tree.label.info, tree = tree)
    new.order <- order(as.numeric(organize.list(node.inf)[,1]), na.last = FALSE)
    nodes <- nodes[new.order]
    nodespatterns <- nodespatterns[new.order]
  }
  for(i in 1:length(nodes)){
    names.descendants <- names(phylobase::descendants(tree1, nodes[i], "all"))
    if(include.node){
      names.descendants <- c(nodes[i], names.descendants)
    }
    edge.inf <- lapply(names.descendants, tree.label.info, tree = tree)
    res[as.numeric(organize.list(edge.inf)[,1])] <- nodespatterns[i]
  }
  return(res)
}