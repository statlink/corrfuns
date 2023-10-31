\name{Confidence interval for the percentage of variance retained by the first k components}
\alias{eigci}
\title{
Confidence interval for the percentage of variance retained by the first \eqn{\kappa} components
}

\description{
Confidence interval for the percentage of variance retained by the first \eqn{\kappa} components.
}

\usage{
eigci(x, k, alpha = 0.05, B = 1000, graph = TRUE)
}

\arguments{
\item{x}{
A numerical matrix with more rows than columns.
}
\item{k}{
The number of principal components to use.
}
\item{alpha}{This is the significance level. Based on this, an \eqn{(1-\alpha)\%} confidence interval will be computed.
}
\item{B}{
The number of bootstrap samples to generate.
}
\item{graph}{
Should the plot of the bootstrap replicates appear? Default value is TRUE.
}
}

\details{
The algorithm is taken by Mardia Kent and Bibby (1979, pg. 233--234). The percentage retained by the fist \eqn{\kappa} principal components denoted by \eqn{\hat{\psi}} is equal to
\deqn{
\hat{\psi}=\frac{ \sum_{i=1}^{\kappa}\hat{\lambda}_i }{\sum_{j=1}^p\hat{\lambda}_j },
}
where \eqn{\hat{\psi}} is asymptotically normal with mean \eqn{\psi} and variance
\deqn{
\tau^2 & = & \frac{2}{\left(n-1\right)\left(tr\pmb{\Sigma} \right)^2}\left[ \left(1-\psi\right)^2\left(\lambda_1^2+...+\lambda_k^2\right)+
\psi^2\left(\lambda_{\kappa+1}^2+...\lambda_p^2\right) \right] \\
& = & \frac{2tr\pmb{\Sigma}^2}{\left(n-1\right)\left(tr\pmb{\Sigma} \right)^2}\left(\psi^2-2a\psi+a\right),
}
where
\eqn{a=\left( \lambda_1^2+...+\lambda_k^2\right)/\left( \lambda_1^2+...+\lambda_p^2\right)}
and \eqn{\text{tr}\pmb{\Sigma}^2=\lambda_1^2+...+\lambda_p^2}

The bootstrap version provides an estimate of the bias, defined as \eqn{\hat{\psi}_{boot}-\hat{\psi}} and confidence intervals calculated via the percentile method and via the standard (or normal) method Efron and Tibshirani (1993). The funciton gives the option to perform bootstrap.
}

\value{
A list including:
\item{res}{
If B=1 (no bootstrap) a vector with the esimated percentage of variance due to the first \eqn{k} components, \eqn{\hat{\psi}} and its associated asymptotic \eqn{(1-\alpha)\%} confidence interval.
If B>1 (bootstrap) a vector with: the esimated percentage of variance due to the first \eqn{k} components, \eqn{\hat{\psi}}, its bootstrap estimate and its bootstrap estimated bias.
}
\item{ci}{
This appears if B>1 (bootstrap). The standard bootstrap and the empirical bootstrap \eqn{(1-\alpha)\%} confidence interval for \eqn{\psi}.
}

Futher, if B>1 and "graph" was set equal to TRUE, a histogram with the bootstrap \eqn{\hat{\psi}} values, the observed \eqn{\hat{\psi}} value and its bootstrap estimate.
}

\references{
Mardia K.V., Kent, J.T. and Bibby, J.M. (1979). Multivariate Analysis.
London: Academic
Press.

Efron B. and Tibshirani R. J. (1993). An introduction to the bootstrap.
Chapman & Hall/CRC.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pc.choose}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
eigci(x, k = 2, B = 1)
}

