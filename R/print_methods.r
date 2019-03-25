#' Print Loadings in Factor Analysis using fad.
#' @rdname loadings
#' @aliases loadings print.loadings print.fad
#' @description Extract or print loadings in factor analysis from package fad.
#' @param x An object of class "fad" or the loadings component of a "fad" object.
#' @param digits number of decimal places to use in printing uniquenesses and loadings.
#' @param loadings smaller than this (in absolute value) are suppressed.
#' @param sort logical. If true, the variables are sorted by their importance on each factor.
#' @param cutoff loadings smaller than this (in absolute value) are suppressed.
#' Each variable with any loading larger than 0.5 (in modulus) is assigned to
#' the factor with the largest loading, and the variables are printed in the order of
#' the factor they are assigned to, then those unassigned.
#' @param \dots further arguments for other methods, ignored for loadings.


#' @export loadings
#' @rdname loadings
loadings <- function(x,...)
{
  x$loadings
}

#' @rdname loadings
#' @export
print.loadings <- function(x, digits = 3L, cutoff = 0.1, sort = FALSE, ...)
{
  Lambda <- unclass(x)
  p <- nrow(Lambda)
  factors <- ncol(Lambda)
  if (sort) {
    mx <- max.col(abs(Lambda))
    ind <- cbind(1L:p, mx)
    mx[abs(Lambda[ind]) < 0.5] <- factors + 1
    Lambda <- Lambda[order(mx, 1L:p),]
  }
  cat("\nLoadings:\n")
  fx <- setNames(format(round(Lambda, digits)), NULL)
  nc <- nchar(fx[1L], type="c")
  fx[abs(Lambda) < cutoff] <- strrep(" ", nc)
  print(fx, quote = FALSE, ...)
  vx <- colSums(x^2)
  varex <- rbind("SS loadings" = vx)
  if(is.null(attr(x, "covariance"))) {
    varex <- rbind(varex, "Proportion Var" = vx/p)
    if(factors > 1)
      varex <- rbind(varex, "Cumulative Var" = cumsum(vx/p))
  }
  cat("\n")
  print(round(varex, digits))
  invisible(x)
}


#' @export
#' @rdname loadings
print.fad <- function(x, digits = 3, ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Uniquenesses:\n")
  print(round(x$uniquenesses, digits), ...)
  print(x$loadings, digits = digits, ...)
  if (!is.null(x$rotmat)){

    tmat <- solve(x$rotmat)
    R <- tmat %*% t(tmat)
    factors <- x$factors
    rownames(R) <- colnames(R) <- paste0("Factor", 1:factors)

    if (TRUE != all.equal(c(R), c(diag(factors)))){
      cat("\nFactor Correlations:\n")
      print(R, digits=digits, ...)
    }


  }

  if(!is.na(x$BIC))    cat("\nThe BIC is: ",x$BIC)

  invisible(x)
}
