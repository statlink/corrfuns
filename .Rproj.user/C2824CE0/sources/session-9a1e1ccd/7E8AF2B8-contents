pcormat <- function(x, type = "pearson") {
  ## x is a matrix with data
  ## type can be either "pearson" or "spearman"
  ## it can of course be "kendall" but I have not examined it
  ## in other functions
  r <- cor(x, method = type) ## correlation matrix of x
  r2 <-  - chol2inv( chol(r) )
  diag(r2) <-  - diag(r2)
  cov2cor(r2)
}
