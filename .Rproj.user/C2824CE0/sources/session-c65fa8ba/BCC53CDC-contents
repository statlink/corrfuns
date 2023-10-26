bootcor2 <- function(x, B = 999) {
  x <- as.matrix(x)
  s <- cov(x)
  n <- dim(x)[1]
  eig <- eigen(s, symmetric = TRUE)
  lam <- eig$values
  vec <- eig$vectors
  ## The next row makes the correlation matrix
  ## equal to the identity matrix, thus rho = 0
  z <- x %*% vec %*% ( t(vec) / sqrt(lam) )
  rb <- numeric(B)
  r <- cor(x[, 1], x[, 2])
  stat <- 0.5 * log( (1 + r)/(1 - r) )  ## the test statistic
  for (i in 1:B) {
    nu <- sample(1:n, replace = TRUE)
    rb[i] <- cor( z[nu, 1], z[nu, 2] )
  }

  tb <- 0.5 * log( (1 + rb)/(1 - rb) )
  pvalue <- ( sum( abs(tb) >  abs(stat) ) + 1)/(B + 1)  ## bootstrap p-value

  result <- c(r, pvalue)
  names(result) <- c('correlation', 'p-value')
  result

}
