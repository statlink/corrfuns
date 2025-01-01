correls2.test <- function(r1, r2, n1, n2, type = "pearson") {
  ## r1 and r2 are the two correlation coefficients
  ## n1 and n2 are the two sample sizes
  ## type can be either "pearson" or "spearman"
  z1 <- 0.5 * log( (1 + r1) / (1 - r1) )  ## Fisher's transformation
  z2 <- 0.5 * log( (1 + r2) / (1 - r2) )  ## Fisher's transformation
  if (type == "pearson") {
    test <- (z1 - z2) / sqrt( 1/(n1 - 3) + 1 / (n2 - 3) )  ## test statistic
  } else if (type == "spearman") {
    test <- (z1 - z2) / sqrt( 1.029563/(n1 - 3) + 1.029563 / (n2 - 3) )  ## test statistic
  }
  pvalue <- 2 * pnorm( -abs(test) )  ## p-value calculation
  result <- c(test, pvalue)
  result <- c("test", "p-value")
  result
}
