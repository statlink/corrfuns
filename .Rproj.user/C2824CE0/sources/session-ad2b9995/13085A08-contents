partialcor2 <- function(y, x, z, type = "pearson", rho = 0, alpha = 0.05) {

  n <- length(y)  ## sample size
  d <- dim(z)[2]  ## dimensionality of z
  res <- resid( lm( cbind(y, x) ~., data = as.data.frame(z) ) )  ## residuals of y and x on z
  r <- cor(res, method = type)[1, 2]  ## partial correlation of y and x conditioning on z
  zh0 <- 0.5 * log( (1 + rho) / (1 - rho) )  ## Fisher's transform for Ho
  zh1 <- 0.5 * log( (1 + r) / (1 - r) )  ## Fisher's transform for H1
  if (type == "pearson") {
    se <- 1/sqrt(n - d - 3)  ## standard error for Fisher's transform under Ho
  } else if ( type == "spearman" ){
    se <-  1.029563 / sqrt(n - d - 3)  ## standard error for Fisher's transform under Ho
  }

  test <- (zh1 - zh0)/se  ## test statistic
  pvalue <- 2 * pt( abs(test), n - d - 3, lower.tail = FALSE )  ## p-value
  zL <- zh1 - qt(1 - alpha/2, n - d - 3) * se
  zH <- zh1 + qt(1 - alpha/2, n - d - 3) * se
  fishL <- (exp(2 * zL) - 1)/(exp(2 * zL) + 1)  ## lower confidence limit
  fishH <- (exp(2 * zH) - 1)/(exp(2 * zH) + 1)  ## upper confidence limit
  ci <- c(fishL, fishH)
  names(ci) <- paste( c( alpha/2 * 100, (1 - alpha/2) * 100 ), "%", sep = "" )

  result <- c(r, pvalue)
  names(result) <- c('correlation', 'p-value')
  list(result = result, ci = ci)
}
