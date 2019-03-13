\name{main}
\alias{main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit an IV regression model with lasso and graph-constrained regularization
}
\description{
This function is a two-step procedure to achieve variable selection and estimation under an IV regression framework by incorporating gene regulatory network structures. Because of high correlations between coefficients, in the first step, the LASSO-type regression is applied to estimate the coefficient matrix. And in the second step, two penalties (LASSO and graph-constrained penalties) are applied to encourage the shrinking coefficients of a pair of connected variable in a network, thus achieving a grouping effect.
}
\usage{
main(y, x, S, G,is.graph="IVGC",nlambda=50,nalpha=50,nfold=5,epsilon=1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The response, a vector of size \eqn{n}
}
  \item{x}{
The expression data dimension \eqn{n\times p} 
}
  \item{S}{
The SNP variable with dimension \eqn{n\times q}
}
  \item{G}{
The graph structure of X with dimension \eqn{p\times p}
}
  \item{is.graph}{
The method including "IVGC", "IV", "both". "IVGC" denotes IV regression incorporating network structures ( Gao et al. 2019); "IV" denotes IV regression without network structures; "both" outputs the results by "IVGC" and "IV". Default is "IVGC".
  }
  \item{nlambda}{
The length of tuning parameter \eqn{\lambda} which is selected automatically. Default is 50.
}
  \item{nalpha}{
The length of tuning parameter \eqn{\alpha} which is selected automatically. Default is 50.
}
  \item{nfold}{
The number of folds for the tuning selection by CV. Default is 5.
}
  \item{epsilon}{
The perturbation when the proposed algorithm is used. Default is epsilon=1e-6.
}
}
\details{
The simulated data set is saved in simulatedData.RData, which is list with four values Y, X, gene (S) and graph (G). 
}
\value{
\item{beta_ivgc}{the beta estimator by the proposed method}
\item{beta_iv}{the beta estimator by IV method without graph structure}
\item{Gamma_hat}{the Gamma estimator by grouped lasso}
}
\references{
Gao, B., Liu, X., Li, H. and Cui, Y. (2019+) Integrative Analysis of Genetical Genomics Data Incorporating Network Structures. Manuscript.
}
\author{
Gao Bin, Xu Liu (liu.xu@sufe.edu.cn), Hongzhe Li, and Yuehua Cui
}


\examples{
data(simulatedData)
y <- simulatedData$Y
x <- simulatedData$X
S <- simulatedData$gene
G <- simulatedData$graph
fit <- main(y, x, S, G)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IV regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ network structure }% __ONLY ONE__ keyword per line
