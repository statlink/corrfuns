\name{Permutation based EL and EEL test for a correlation coefficient}
\alias{perm.elcortest}
\alias{perm.eelcortest}
\title{
Permutation based EL and EEL test for a correlation coefficient
}

\description{
Permutation based EL and EEL test for a correlation coefficient.
}

\usage{
perm.elcortest(y, x, tol = 1e-07, B = 999)
perm.eelcortest(y, x, tol = 1e-07, B = 999)
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A numerical vector.
}
\item{tol}{
The tolerance vlaue to terminate the Newton-Raphson algorithm.
}
\item{B}{
The numer of permutations to perform.
}
}

\details{
These functions compute the permutation based p-value for the hypothesis test that the Pearson correlation coefficient is equal to zero using the empirical likelihood or the exponential empirical likelihood.
}

\value{
A vector with two values, the test statistic and its associated permutation based p-value.
}

\references{
Efron B. (1981) Nonparametric standard errors and confidence intervals. Canadian Journal of
Statistics, 9(2): 139--158.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{el.cor.test}, \link{permcor}
}
}

\examples{
perm.eelcortest( iris[, 1], iris[, 2])
}


