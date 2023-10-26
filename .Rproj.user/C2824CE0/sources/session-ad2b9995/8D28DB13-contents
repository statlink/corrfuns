\name{Asymptotic p-value for a correlation coefficient}
\alias{correl}
\title{
Asymptotic p-value for a correlation coefficient
}

\description{
Asymptotic p-value a correlation coefficient.
}

\usage{
correl(y, x, type = "pearson", rho = 0, alpha = 0.05)
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A numerical vector.
}
\item{type}{
The type of correlation coefficient to compute, "pearson" or "spearman".
}
\item{rho}{
The hypothesized value of the true partial correlation.
}
\item{alpha}{
The significance level.
}
}

\details{
Given a correlation matrix the function will caclulate the partial correlation between variables indx and indy conditioning on variable(s) indz. The logarithm of the p-value is also returned.
}

\value{
The partial correlation coefficient and the logged p-value for the test of no association.
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
\code{ \link{correl}, \link{permcorrels}
}
}

\examples{
y <- rnorm(40)
x <- matrix(rnorm(40 * 1000), ncol = 1000)
a <- correls(y, x )
}


