\encoding{UTF-8}
\name{mantel.residuals}
\alias{mantel.residuals}

\title{Residuals from Mantel test}

\description{

Function to extract residuals from Mantel test.

}

\usage{
mantel.residuals(A,B)
}

\arguments{
\item{A,B}{Dissimilarity matrices or dist objects.}
}

\value{
\item{Mantel_statistic_r}{The Mantel statistic.}
\item{Residuals}{Residuals extracted from the Mantel test (Class dist).}
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{
\code{\link{mantel}}
}

\examples{
X<-vegdist(matrix(sample(1:10,20,replace=TRUE),5,5))
Y<-vegdist(matrix(sample(1:10,30,replace=TRUE),5,6))
mantel(X,Y)
mantel.residuals(X,Y)
}

\keyword{daee}
