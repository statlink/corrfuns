partialcor <- function(R, indx, indy, indz, n) {
  ## R is a correlation matrix
  ## i and j denote the two variables whose conditional correlation is to be estimated
  ## k denotes the set of conditioning variables
  d <- length(indz)
  if ( d == 1 ) {
    a1 <- R[indx, indy]     ;    a2 <- R[indx, indz]
    a3 <- R[indy, indz]
    r <- (a1 - a2 * a3) / sqrt( (1 - a3^2) * (1 - a2^2) )

  } else if ( d > 1 ) {
    rho <- solve( R[c(indx, indy, indz), c(indx, indy, indz)] )
    r <-  - rho[1, 2] / sqrt(rho[1, 1] * rho[2, 2])
  }
  if ( abs(r) > 1 ) r <- 0.99999
  stat <- 0.5 * log( (1 + r) / (1 - r) ) * sqrt( n - d - 3 )
  pvalue <- 2 * pt( abs(stat), n - d - 3, lower.tail = FALSE )
  res <- c(r, pvalue)
  names(res) <- c("partial cor", "p-value")
  res
}

