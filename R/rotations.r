

# The code for varimax is modified from the base R package's code:
# `t(X) %*%` replaced by crossprod(x,)
varimax <- function(x, normalize = TRUE, eps = 1e-5)
{
  nc <- ncol(x)
  if(nc < 2) return(x)
  if(normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  p <- nrow(x)
  TT <- diag(nc)
  d <- 0
  for(i in 1L:1000L) {
    z <- x %*% TT
    B  <- crossprod(x, z^3 - z %*% diag(drop(rep(1, p) %*% z^2))/p)
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    if(d < dpast * (1 + eps)) break
  }
  z <- x %*% TT
  if(normalize) z <- z * sc
  dimnames(z) <- dimnames(x)
  class(z) <- "loadings"
  list(loadings = z, rotmat = TT)
}


# The code for promax is modified from the base R package's code:
# U %*% diag(sqrt(d)) changed so that diag(sqrt(d)) is not created.
promax <- function(x, m = 4)
{
  if(ncol(x) < 2) return(x)
  dn <- dimnames(x)
  xx <- varimax(x)
  x <- xx$loadings
  Q <- x * abs(x)^(m-1)
  U <- lm.fit(x, Q)$coefficients
  d <- diag(solve(crossprod(U)))
  U <- .postmdiag(U , sqrt(d))
  dimnames(U) <- NULL
  z <- x %*% U
  U <- xx$rotmat %*% U
  dimnames(z) <- dn
  class(z) <- "loadings"
  list(loadings = z, rotmat = U)
}

