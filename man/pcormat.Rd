\name{Partial correlation matrix}
\alias{pcormat}
\title{
Partial correlation matrix
}

\description{
Partial correlation matrix.
}

\usage{
pcormat(x, type = "pearson")
}

\arguments{
\item{x}{
A numerical matrix with two columns.
}
\item{type}{
The type of the partial correlation matrix to compute, "pearson" or "spearman".
}
}

\details{
Given a correlation matrix the function will caclulate the partial correlation between variables indx and indy conditioning on variable(s) indz. The logarithm of the p-value is also returned.
}

\value{
The partial correlation coefficient and the p-value for the test of no association.
}

%\references{
%
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{partialcor}, \link{partialcor2}
}
}

\examples{
pcormat( iris[, 1:4] )
}


