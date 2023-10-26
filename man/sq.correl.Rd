\name{Squared multivariate correlation between two sets of variables}
\alias{sq.correl}
\title{
Squared multivariate correlation between two sets of variables
}

\description{
Squared multivariate correlation between two sets of variables.
}

\usage{
sq.correl(y, x)
}

\arguments{
\item{y}{
A numerical matrix.
}
\item{x}{
A numerical matrix.
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
\code{ \link{correl}
}
}

\examples{
sq.correl( iris[, 1:2], iris[, 3:4] )
}


