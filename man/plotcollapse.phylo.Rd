% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotcollapse.phylo.R
\encoding{UTF-8}
\name{plotcollapse.phylo}
\alias{plotcollapse.phylo}
\title{Plot phylogenetic tree with nodes collapsed.}
\usage{
plotcollapse.phylo(
  tree,
  nodes,
  show.tip.label = TRUE,
  show.node.label = FALSE,
  edge.color = "black",
  edge.width = 1,
  edge.lty = 1,
  polygon.color = "gray",
  density = NULL,
  border = NULL,
  nhchar = 2,
  text.nodes.color = NULL,
  ...
)
}
\arguments{
\item{tree}{phylogeny as an object of class "phylo".}

\item{nodes}{a vector with node label to collapse the edges.}

\item{show.tip.label}{a logical indicating whether to show the tip labels on the phylogeny. The same used in the function plot.phylo (default show.tip.label = TRUE).}

\item{show.node.label}{a logical indicating whether to show the node labels on the phylogeny. The same used in the function plot.phylo (default show.node.label = FALSE).}

\item{edge.color}{a vector of colours giving the colours used to draw the branches of the plotted phylogeny. The same used in the function plot.phylo (default edge.color = "black").}

\item{edge.width}{a numeric vector giving the width of the branches of the plotted phylogeny. The same used in the function plot.phylo (default edge.width = 1).}

\item{edge.lty}{a numeric vector giving the line types of the branches of the plotted phylogeny. The same used in the function plot.phylo (default edge.lty = 1).}

\item{polygon.color}{a vector of colours for each node collapsed, this is used to filling each polygon (default polygon.color = "gray").}

\item{density}{a numeric vector with density of shading lines (lines per inch) for each node collapsed, this is used to filling each polygon (default density = NULL).}

\item{border}{a vector of colours for each node collapsed, this is used to draw the border in each polygon (default border = NULL).}

\item{nhchar}{number of hidden characters in tip.labels. This adjusts the margin of tip or nodes names are not shown completely (default nhchar = 2).}

\item{text.nodes.color}{a vector of colours to text label in each nodes collapsed. If NULL color of text is same that tip.label (default text.nodes.color = NULL).}

\item{...}{other arguments to function plot.phylo.}
}
\value{
The plot of phylogenetic trees with nodes collapsed and returns invisibly a list with components which values are those used for the plot. See details in \code{\link{plot.phylo}}.
}
\description{
This function plot phylogenetic tree with nodes collapsed. All arguments used are of the 
function \code{\link{plot.phylo}}. See details.
}
\details{
The function use the function \code{\link{plot.phylo}} to plot the phylonenies. This function only 
hidden all edge, node labels and tips labels of each collapsed node and draw a polygon that represent a 
collapsed node. Type "phylogram" and "cladogram" work properly, other type of phyloneny will be implemented 
in future (Tip labels are not draw in correct position). The function \code{\link{compact.tree}} can be used
to compact tips of phylogenetic tree before plot. The function \code{\link{graphical.node.patterns}} can be
used to define graphical parameters (e.g. color, width and line types) for each edge in phylogenetic tree.

This function require node labels. The function \code{\link{makeNodeLabel}} or \code{\link{node.tree}} can be
used to make node labels before plot, anyway the the function \code{\link{node.tree}} is apply automatically. 
Make sure that that there no conflict with node labels.
}
\examples{

set.seed(10)
tree<-rtree(15)
tree <- makeNodeLabel(tree)
plot.phylo(tree, show.node.label = TRUE)
nodes <- c("Node6", "Node13")

plot.phylo(tree, show.node.label = TRUE)
plotcollapse.phylo(tree, nodes)
plotcollapse.phylo(tree, nodes, nhchar = 5) # nhchar ajust the margin
plotcollapse.phylo(tree, nodes, show.tip.label = FALSE)
plotcollapse.phylo(tree, nodes, direction = "rightwards", polygon.color = c("red","blue"))
plotcollapse.phylo(tree, nodes, direction = "rightwards", polygon.color = c("red","blue"), 
   text.nodes.color = c("red", "blue"))
plotcollapse.phylo(tree, nodes, direction = "upwards", polygon.color = c("red","blue"))
plotcollapse.phylo(tree, nodes, direction = "downwards", polygon.color = c("red","blue"))
plotcollapse.phylo(tree, nodes, direction = "leftwards", polygon.color = c("red","blue"))
plotcollapse.phylo(tree, nodes, direction = "leftwards", polygon.color = c("red","blue"), 
   nhchar = 6)

set.seed(20)
tree <- rtree(15)
tree <- makeNodeLabel(tree)
plot.phylo(tree, show.node.label = TRUE)
tree.compacted <- compact.tree(tree, "Node6")
plot.phylo(tree.compacted, show.node.label = TRUE) # plot tree compacted
nodes <- c("Node5","Node7")
plotcollapse.phylo(tree, nodes, nhchar = 4)
plotcollapse.phylo(tree.compacted, nodes, nhchar = 4)

e.wid <- graphical.node.patterns(tree, nodes, basicpattern = 1, 
    nodespatterns = 3, include.node = TRUE)
e.wid
plotcollapse.phylo(tree, nodes, show.node.label = TRUE, edge.width = e.wid)
e.col <- graphical.node.patterns(tree, nodes, basicpattern = "black", 
    nodespatterns = c("red","blue"), include.node = TRUE)
e.col
plotcollapse.phylo(tree, nodes, show.node.label = TRUE, edge.width = e.wid, edge.color = e.col)

tree <- rtree(50)
tree <- compute.brlen(tree)
tree <- makeNodeLabel(tree)
plot.phylo(tree, type ="fan", show.tip.label = FALSE)
nodes <- "Node30"
plotcollapse.phylo(tree, nodes, type ="fan", show.tip.label = FALSE)

}
\seealso{
\code{\link{plot.phylo}} \code{\link{compact.tree}} \code{\link{graphical.node.patterns}}
}
\author{
Vanderlei Julio Debastiani <vanderleidebastiani@yahoo.com.br>
}
\keyword{daee}
