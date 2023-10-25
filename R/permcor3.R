permcor3 <- function(x, B = 999) {
  ## x is a 2 column matrix containing the data
  ## type can be either "pearson" or "spearman"
  ## R is the number of permutations
  x1 <- x[, 1]      ;     x2 <- x[, 2]
  n <- length(x1)
  m1 <- sum(x1)     ;     m12 <- sum(x1^2)
  m2 <- sum(x2)     ;     m22 <- sum(x2^2)
  up <-  m1 * m2 / n
  down <- sqrt( (m12 - m1^2 / n) * (m22 - m2^2 / n) )
  r <- ( sum(x1 * x2) - up) / down
  test <- log( (1 + r) / (1 - r) )  ## test statistic
  sxy <- numeric(B)
  for (i in 1:B) {
    y1 <- sample(x1, n, replace = FALSE)
    sxy[i] <- sum(y1 * x2)
  }
  rb <- (sxy - up) / down
  tb <- log( (1 + rb) / (1 - rb) )  ## test statistic
  pvalue <- ( sum( abs(tb) > abs(test) ) + 1 ) / (B + 1)  ## permutation p-value
  res <- c( r, pvalue )
  names(res) <- c('correlation', 'p-value')
  res
}
