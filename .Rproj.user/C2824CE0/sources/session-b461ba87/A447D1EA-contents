\name{Empirical and exponential empirical likelihood test for a correlation coefficient}
\alias{el.cor.test}
\alias{eel.cor.test}
\title{
Empirical and exponential empirical likelihood test for a correlation coefficient
}

\description{
Empirical and exponential empirical likelihood test for a correlation coefficient.
}

\usage{
el.cor.test(y, x, rho = 0, tol = 1e-07)
el.cor.test(y, x, rho = 0, tol = 1e-07)
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A numerical vector.
}
\item{rho}{
The hypothesized value of the true partial correlation.
}
\item{tol}{
The tolerance vlaue to terminate the Newton-Raphson algorithm.
}
}

\details{

}

\value{
A list including:
\item{iters}{
The number of iterations required by the Newton-Raphson.
}
\item{info}{
A vector with three values, the value of \eqn{\lambda}, the test statistic and its associated asymptotic p-value.
}
\item{p}{
The probabilities of the EL or of the EEL.
}
}

\references{
Efron B. and Tibshirani R.J. (1993). An introduction to the bootstrap. Chapman & Hall/CRC.

Fieller E.C., Hartley H.O. and Pearson E.S. (1957). Tests for rank correlation coefficients.
I. Biometrika, 44(3/4): 470--481.

Fieller E.C. and Pearson E.S. (1961). Tests for rank correlation coefficients: II. Biometrika,
48(1/2): 29--40.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{correls}, \link{permcorrels}
}
}

\examples{
a <- correl( iris[, 1], iris[, 2] )
}


