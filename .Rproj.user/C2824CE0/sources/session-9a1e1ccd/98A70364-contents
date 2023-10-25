\name{Hypothesis test for equality of two correlation coefficients}
\alias{correls2.test}
\title{
Hypothesis test for equality of two correlation coefficients
}

\description{
Hypothesis test for equality of two correlation coefficients.
}

\usage{
correls2.test(r1, r2, n1, n2, type = "pearson")
}

\arguments{
\item{r1}{
The value of the first correlation coefficient.
}
\item{r2}{
The value of the second correlation coefficient.
}
\item{n1}{
The sample size of the first sample from which the first correlation coefficient was computed.
}
\item{n2}{
The sample size of the second sample from which the first correlation coefficient was computed.
}
\item{type}{
The type of correlation coefficients, "pearson" or "spearman".
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
\code{ \link{correl}, \link{correls}
}
}

\examples{
y <- rnorm(40)
x <- matrix(rnorm(40 * 1000), ncol = 1000)
a <- correls(y, x )
}


