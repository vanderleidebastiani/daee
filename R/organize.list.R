#' @title Organize a list with data.frame/matrix in a single matrix
#'
#' @description Function organize a list of data.frame or matrix in a single matrix. The function combine the data.frame/matrix by
#' columns, thus all objects in the list must be equal number of columns
#'
#' @encoding UTF-8
#' @param x A list with data.frame or matrix
#' @param force.names Logical argument (TRUE or FALSE) to specify if new row names are created, using list name as prefix.
#' @return A matrix combined.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @keywords daee
#' @examples
#' A<-matrix(1,10,6)
#' rownames(A) <- rownames(A, do.NULL = FALSE, prefix = "A.")
#' B<-matrix(2,3,6)
#' rownames(B) <- rownames(B, do.NULL = FALSE, prefix = "B.")
#' C<-matrix(3,15,6)
#' rownames(C) <- rownames(C, do.NULL = FALSE, prefix = "C.")
#' res<-list(as.data.frame(A), as.data.frame(B), as.data.frame(C))
#' organize.list(RES)
#' organize.list(RES, force.names = TRUE)
#' @export
organize.list <- function(x, force.names = FALSE){
	if(!inherits(x, "list")){
		stop("\n x must be a list \n")
	}
	n.row.list <- sapply(x, nrow)
	n.col.list <- sapply(x, ncol)
	d <- cumsum(c(0, n.row.list))
	if(length(unique(n.col.list))!=1){
		stop("\n all data.frame or matrix in x must be equal number of columns \n")
	}
	res <- matrix(NA, sum(n.row.list), n.col.list[1])
	colnames(res) <- colnames(x[[1]])
	rownames(res) <- rownames(res, do.NULL = FALSE)
	if(is.null(names(x))){
		names(x) <- paste("list", 1:length(x), sep = ".")
	}
	for(i in 1:length(x)){
		if(n.row.list[i]>0){
			res[(d[i]+1):d[i+1],] <- as.matrix(x[[i]])
			rownames(res)[(d[i]+1):d[i+1]] <- rownames(x[[i]], do.NULL = FALSE, prefix = paste(names(x)[i], ".", sep = ""))
			if(force.names){
				rownames(res)[(d[i]+1):d[i+1]] <- paste(names(x)[i], ".", rownames(res)[(d[i]+1):d[i+1]], sep = "")
			}
		}
	}
	return(res)
}