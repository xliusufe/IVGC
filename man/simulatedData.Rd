\name{simulatedData}
\alias{simulatedData}
\docType{data}
\title{
Genetical Genomics Data Incorporating Network Structures
}
\description{
This is a simulated data on gene expression. Certain correlation structures on genes are also imposed to show the impact of network structures on variable selection and estimation.
}
\usage{
data("simulatedData")
}
\details{
We simulated a total of p variables of gene expression and considered three group structures on coefficients with five variables in each group. Within each group, the variables are correlated. The strength of the correlation was controlled by the number of effective SNPs they had in common.

Y: the response, a vector of size 200

X: the expression data dimension \eqn{200\times 100} 

gene: the SNP variable with dimension \eqn{200\times 100}

graph: the graph structure of X with dimension \eqn{100\times 100}
}

\references{
Gao, B., Liu, X., Li, H. and Cui, Y. (2019+) Integrative Analysis of Genetical Genomics Data Incorporating Network Structures. Manuscript.
}
\examples{
data(simulatedData)
y <- simulatedData$Y
x <- simulatedData$X
S <- simulatedData$gene
G <- simulatedData$graph
}
\keyword{datasets}
