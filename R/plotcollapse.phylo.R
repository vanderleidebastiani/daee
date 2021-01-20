#' @title Plot phylogenetic tree with nodes collapsed.
#' 
#' @description This function plot phylogenetic tree with nodes collapsed. All arguments used are of the 
#' function \code{\link{plot.phylo}}. See details.
#' 
#' @details The function use the function \code{\link{plot.phylo}} to plot the phylonenies. This function only 
#' hidden all edge, node labels and tips labels of each collapsed node and draw a polygon that represent a 
#' collapsed node. Type "phylogram" and "cladogram" work properly, other type of phyloneny will be implemented 
#' in future (Tip labels are not draw in correct position). The function \code{\link{compact.tree}} can be used
#' to compact tips of phylogenetic tree before plot. The function \code{\link{graphical.node.patterns}} can be
#' used to define graphical parameters (e.g. color, width and line types) for each edge in phylogenetic tree.
#' 
#' This function require node labels. The function \code{\link{makeNodeLabel}} or \code{\link{node.tree}} can be
#' used to make node labels before plot, anyway the the function \code{\link{node.tree}} is apply automatically. 
#' Make sure that that there no conflict with node labels.
#' 
#' @encoding UTF-8
#' @import ape
#' @importFrom phylobase phylo4 descendants
#' @importFrom graphics polygon text
#' @param tree phylogeny as an object of class "phylo".
#' @param nodes a vector with node label to collapse the edges.
#' @param show.tip.label a logical indicating whether to show the tip labels on the phylogeny. The same used in the function plot.phylo (default show.tip.label = TRUE).
#' @param show.node.label a logical indicating whether to show the node labels on the phylogeny. The same used in the function plot.phylo (default show.node.label = FALSE).
#' @param edge.color a vector of colours giving the colours used to draw the branches of the plotted phylogeny. The same used in the function plot.phylo (default edge.color = "black").
#' @param edge.width a numeric vector giving the width of the branches of the plotted phylogeny. The same used in the function plot.phylo (default edge.width = 1).
#' @param edge.lty a numeric vector giving the line types of the branches of the plotted phylogeny. The same used in the function plot.phylo (default edge.lty = 1).
#' @param polygon.color a vector of colours for each node collapsed, this is used to filling each polygon (default polygon.color = "gray").
#' @param density a numeric vector with density of shading lines (lines per inch) for each node collapsed, this is used to filling each polygon (default density = NULL).
#' @param border a vector of colours for each node collapsed, this is used to draw the border in each polygon (default border = NULL).
#' @param nhchar number of hidden characters in tip.labels. This adjusts the margin of tip or nodes names are not shown completely (default nhchar = 2).
#' @param text.nodes.color a vector of colours to text label in each nodes collapsed. If NULL color of text is same that tip.label (default text.nodes.color = NULL).
#' @param ... other arguments to function plot.phylo.
#' @return The plot of phylogenetic trees with nodes collapsed and returns invisibly a list with components which values are those used for the plot. See details in \code{\link{plot.phylo}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{plot.phylo}} \code{\link{compact.tree}} \code{\link{graphical.node.patterns}}
#' @keywords daee
#' @examples
#' 
#' set.seed(10)
#' tree<-rtree(15)
#' tree <- makeNodeLabel(tree)
#' plot.phylo(tree, show.node.label = TRUE)
#' nodes <- c("Node6", "Node13")
#' 
#' plot.phylo(tree, show.node.label = TRUE)
#' plotcollapse.phylo(tree, nodes)
#' plotcollapse.phylo(tree, nodes, nhchar = 5) # nhchar ajust the margin
#' plotcollapse.phylo(tree, nodes, show.tip.label = FALSE)
#' plotcollapse.phylo(tree, nodes, direction = "rightwards", polygon.color = c("red","blue"))
#' plotcollapse.phylo(tree, nodes, direction = "rightwards", polygon.color = c("red","blue"), 
#'    text.nodes.color = c("red", "blue"))
#' plotcollapse.phylo(tree, nodes, direction = "upwards", polygon.color = c("red","blue"))
#' plotcollapse.phylo(tree, nodes, direction = "downwards", polygon.color = c("red","blue"))
#' plotcollapse.phylo(tree, nodes, direction = "leftwards", polygon.color = c("red","blue"))
#' plotcollapse.phylo(tree, nodes, direction = "leftwards", polygon.color = c("red","blue"), 
#'    nhchar = 6)
#' 
#' set.seed(20)
#' tree <- rtree(15)
#' tree <- makeNodeLabel(tree)
#' plot.phylo(tree, show.node.label = TRUE)
#' tree.compacted <- compact.tree(tree, "Node6")
#' plot.phylo(tree.compacted, show.node.label = TRUE) # plot tree compacted
#' nodes <- c("Node5","Node7")
#' plotcollapse.phylo(tree, nodes, nhchar = 4)
#' plotcollapse.phylo(tree.compacted, nodes, nhchar = 4)
#' 
#' e.wid <- graphical.node.patterns(tree, nodes, basicpattern = 1, 
#'     nodespatterns = 3, include.node = TRUE)
#' e.wid
#' plotcollapse.phylo(tree, nodes, show.node.label = TRUE, edge.width = e.wid)
#' e.col <- graphical.node.patterns(tree, nodes, basicpattern = "black", 
#'     nodespatterns = c("red","blue"), include.node = TRUE)
#' e.col
#' plotcollapse.phylo(tree, nodes, show.node.label = TRUE, edge.width = e.wid, edge.color = e.col)
#' 
#' tree <- rtree(50)
#' tree <- compute.brlen(tree)
#' tree <- makeNodeLabel(tree)
#' plot.phylo(tree, type ="fan", show.tip.label = FALSE)
#' nodes <- "Node30"
#' plotcollapse.phylo(tree, nodes, type ="fan", show.tip.label = FALSE)
#' 
#' @export
plotcollapse.phylo <- function(tree, nodes, show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", edge.width = 1, edge.lty = 1, polygon.color = "gray", density = NULL, border = NULL, nhchar = 2, text.nodes.color = NULL, ...){
  nnodes <- length(nodes)
  if (!inherits(tree, "phylo")) {
    stop("\n Object tree is not of class phylo \n")
  }
  tree <- node.tree(tree)$tree
  tree1 <- phylobase::phylo4(tree)
  check.nodes <- unlist(lapply(nodes, function(node) names(phylobase::descendants(node, phy = tree1, type= "all"))))
  if(any(check.nodes %in% nodes)){
    stop("nodes cannot be collapsed one inside other")
  }
  node.inf <- lapply(nodes, tree.label.info, tree = tree)
  node.inf <- as.numeric(organize.list(node.inf)[,1])
  if(any(is.na(node.inf))){
    stop("node to collapse must not be the root") 
  }
  if(length(edge.color)!=nrow(tree$edge)){
    edge.color <- rep(edge.color, nrow(tree$edge))
  }
  if(length(edge.width)!=nrow(tree$edge)){
    edge.width <- rep(edge.width, nrow(tree$edge))
  }
  if(length(edge.lty)!=nrow(tree$edge)){
    edge.lty <- rep(edge.lty, nrow(tree$edge))
  }
  if(length(polygon.color)!=nnodes){
    polygon.color <- rep(polygon.color, nnodes)
  }
  if(length(density)!=nnodes){
    density <- rep(density, nnodes)
  }
  if(length(border)!=nnodes){
    border <- rep(border, nnodes)
  }
  if(nhchar<0){
    stop("nhchar must not be negative")
  }
  if(!nhchar%%1==0){
    stop("\n nhchar must an integer \n")
  }
  if(!is.null(text.nodes.color)){
    if(length(text.nodes.color)!=nnodes){
      text.nodes.color <- rep(text.nodes.color, nnodes)
    }
  }
  tree.temp <- tree
  for(i in 1:length(nodes)){
    tips.node <- names(phylobase::descendants(tree1, nodes[i], "all"))
    edge.width[as.numeric(organize.list(lapply(tips.node, tree.label.info, tree = tree))[,1])] <- 0
    tree$node.label[which(tree$node.label %in% c(nodes[i], tips.node))] <- ""
    tips.node <- tips.node[tips.node %in% tree$tip.label]
    ncn <- nchar(nodes[i])
    nct <- nchar(tips.node)
    space <- sapply(nct, function(x) paste(rep(" ", max(x,ncn)+nhchar), collapse = ""))
    tree$tip.label[which(tree$tip.label  %in% tips.node)] <- space
  }
  last_plot.phylo <- plot.phylo(tree, show.tip.label = show.tip.label, show.node.label = show.node.label, edge.color = edge.color, edge.width = edge.width, edge.lty = edge.lty, ...)
  last_phylo_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  for(i in 1:length(nodes)){
    edge.group <- which(tree.temp$tip.label %in% names(descendants(tree1, nodes[i], "tips")))
    edge.node <- which(c(tree.temp$tip.label,tree.temp$node.label) %in% nodes[i])
    x.node <- last_phylo_plot$xx[edge.node]
    y.node <- last_phylo_plot$yy[edge.node]
    if(last_phylo_plot$type == "phylogram" || last_phylo_plot$type == "cladogram"){
      if (last_phylo_plot$direction == "rightwards"){
        x.coords <- c(x.node, max(last_phylo_plot$xx[edge.group]), max(last_phylo_plot$xx[edge.group]))
        y.coords <- c(y.node, max(last_phylo_plot$yy[edge.group]), min(last_phylo_plot$yy[edge.group]))
      }
      if (last_phylo_plot$direction == "upwards"){
        x.coords <- c(x.node, max(last_phylo_plot$xx[edge.group]), min(last_phylo_plot$xx[edge.group]))
        y.coords <- c(y.node, max(last_phylo_plot$yy[edge.group]), max(last_phylo_plot$yy[edge.group]))
      }
      if (last_phylo_plot$direction == "leftwards"){
        x.coords <- c(x.node, min(last_phylo_plot$xx[edge.group]), min(last_phylo_plot$xx[edge.group]))
        y.coords <- c(y.node, max(last_phylo_plot$yy[edge.group]), min(last_phylo_plot$yy[edge.group]))
      }
      if (last_phylo_plot$direction == "downwards"){
        x.coords <- c(x.node, max(last_phylo_plot$xx[edge.group]), min(last_phylo_plot$xx[edge.group]))
        y.coords <- c(y.node, min(last_phylo_plot$yy[edge.group]), min(last_phylo_plot$yy[edge.group]))
      }
    } else{
      x.coords <- c(x.node, last_phylo_plot$xx[edge.group])
      y.coords <- c(y.node, last_phylo_plot$yy[edge.group])
    }
    graphics::polygon(x.coords, y.coords, col = polygon.color[i], density = density[i], border = border[i])
    if(show.tip.label){
      if(is.null(text.nodes.color)){
        text.nodes.color <- last_phylo_plot$tip.color
      }
      if (last_phylo_plot$direction == "rightwards"){
        graphics::text(max(last_phylo_plot$xx[edge.group]), mean(c(max(last_phylo_plot$yy[edge.group]), min(last_phylo_plot$yy[edge.group]))), nodes[i], srt = 0, font = last_phylo_plot$font, cex = last_phylo_plot$cex, adj = 0, col = text.nodes.color[i])
      }
      if (last_phylo_plot$direction == "upwards"){
        graphics::text(mean(c(max(last_phylo_plot$xx[edge.group]), min(last_phylo_plot$xx[edge.group]))), max(last_phylo_plot$yy[edge.group]), nodes[i], srt = 90, font = last_phylo_plot$font, cex = last_phylo_plot$cex, adj = 0, col = text.nodes.color[i])
      }
      if (last_phylo_plot$direction == "leftwards"){
        graphics::text(min(last_phylo_plot$xx[edge.group]), mean(c(max(last_phylo_plot$yy[edge.group]), min(last_phylo_plot$yy[edge.group]))), nodes[i], srt = 0, pos = 2, font = last_phylo_plot$font, cex = last_phylo_plot$cex, adj = 0, col = text.nodes.color[i])
      }
      if (last_phylo_plot$direction == "downwards"){
        graphics::text(mean(c(max(last_phylo_plot$xx[edge.group]), min(last_phylo_plot$xx[edge.group]))), min(last_phylo_plot$yy[edge.group]), nodes[i], srt = 270, font = last_phylo_plot$font, cex = last_phylo_plot$cex, adj = 0, col = text.nodes.color[i])
      }
    }
  }
  invisible(last_plot.phylo)
}