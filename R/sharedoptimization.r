grad_share <- function(fn_gr) {
  env <- new.env()
  env$x_last <- NULL
  env$fn <- function(x) {
    if (is.null(env$x_last) || any(env$x_last != x)) {
      out <- fn_gr(x)
      env$x_last <- x
      env$fn_val <- out[[1]]
      env$gr_val <- out[[2]]
    }
    env$fn_val
  }
  env$gr <- function(x = NULL) {
    if (is.null(env$x_last) || any(env$x_last != x)) {
      out <- fn_gr(x)
      env$x_last <- x
      env$fn_val <- out[[1]]
      env$gr_val <- out[[2]]
    }
    env$gr_val
  }
  env
}


fngr <- function(func, evalForNewX=TRUE,
                 recalculate_indices = c(),
                 check_all=FALSE
                 ) {
  if (evalForNewX == F & length(recalculate_indices) == 0) {
    stop("Values will never be calculated")
  }
  env <- new.env()
  env$x_last <- NULL
  env$f <- function(i, evalForNewX_=evalForNewX,
                    recalculate = any(i==recalculate_indices),
                    check=check_all
                    ) {
    function(x=NULL, check_now=check, recalculate_now=recalculate,
             evalForNewX_now=evalForNewX_
             ) {
      if (recalculate_now ||
          (evalForNewX_now &&
           ((is.null(env$x_last) && !is.null(x)) ||
            (!is.null(env$x_last) && !is.null(x) && any(x != env$x_last))
            ))
        ) {
        out <- func(x)
        env$x_last <- x
        env$out <- out
        out[[1]]
      } else {
        # Can check if evaluated at same value, but will only slow it down
        if (check_now) {
          if (!is.null(x) &&
              !any(is.nan(x)) &&
              (is.null(env$x_last) || any(x != env$x_last))) {
            warning("gr called at different x than fn")
          }
        }
      }
      env$out[[i]]
    }
  }
  env$f
}

optim_share <- function(par, fngr, ...) {
  env <- grad_share(fngr)
  optim(par=par, fn=env$fn, gr=env$gr, ...)
}

