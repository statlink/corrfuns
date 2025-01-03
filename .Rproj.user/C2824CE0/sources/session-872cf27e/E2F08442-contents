perm.eelcortest <- function(y, x, tol = 1e-07, B = 999) {

  funa <- function(z, z2, n, tol) {
    lam1 <- 0
    f <- mean(z)
    der <- ( n * sum(z2) - sum(z)^2 ) / n^2
    lam2 <- lam1 - f / der
    i <- 2
    while (abs(lam1 - lam2) > tol) {
      i <- i + 1
      lam1 <- lam2
      com <- exp( lam1 * z )
      up <- sum( z * com )
      down <- sum(com)
      f <- up / down
      der <- ( sum(z2 * com) * down - up^2 ) / down^2
      lam2 <- lam1 - f / der
    }
    list(iters = i, lam2 = lam2, p = com / down )
  }

  x <- ( x - mean(x) ) / Rfast::Var(x, std = TRUE)
  y <- ( y - mean(y) ) / Rfast::Var(y, std = TRUE)
  z <- x * y   ;   z2 <- z^2

  n <- length(z)

  res <- try( funa(z, z2, n, tol), silent = TRUE )
  if ( identical(class(res), "try-error") ) {
    stat <- 10^5
    pval <- 0
  } else {
    stat <-  -2 * sum( log(n * res$p) )
    tb <- numeric(B)
    for ( i in 1:B ) {
      x1 <- Rfast2::Sample(x, n)
      z <- x1 * y   ;   z2 <- z^2
      res <- try( funa(z, z2, n, tol), silent = TRUE )
      if ( identical( class(res), "try-error") ) {
        tb[i] <- 10^5
      } else  tb[i] <-  -2 * sum( log(n * res$p) )
    }
    pval <- ( sum(tb > stat) + 1 ) / (B + 1)
  }

  res <- c(stat, pval)
  names(res) <- c("statistic", "p-value")
  res
}
