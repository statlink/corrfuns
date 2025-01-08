permcorrels <- function(y, x, B = 999) {

  p <- dim(x)[2]
  n <- length(y)
  r <- as.vector( cor(y, x) )
  test <- log( (1 + r) / (1 - r) )  ## test statistic
  m1 <- sum(y)         ;     m12 <- sum(y^2)
  m2 <- Rfast::colsums(x)     ;     m22 <- Rfast::colsums(x^2)
  up <-  m1 * m2 / n
  down <- sqrt( (m12 - m1^2 / n) * (m22 - m2^2 / n) )
  sxy <- matrix(0, p, B)
  for (i in 1:B) {
    y1 <- sample(y, n)
    sxy[, i] <- Rfast::colsums(y1 * x)
  }
  rb <- (sxy - up) / down
  tb <- log( (1 + rb)/(1 - rb) )  ## test statistic
  pvalue <- ( Rfast::rowsums( abs(tb) > abs(test) ) + 1 ) / (B + 1)  ## permutation p-value
  res <- cbind( r, pvalue )
  colnames(res) <- c('correlation', 'p-value')
  if ( is.null(colnames(x)) ) {
    rownames(res) <- paste("X", 1:dim(x)[2], sep = "")
  } else  rownames(res) <- colnames(x)
  res
}
