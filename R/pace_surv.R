#' Calculate shape of survival over age
#'
#' Calculates a 'shape' value of survival lifespan inequality by comparing the 
#' area under a survival curve (over age) with the area under a constant 
#' survival function.
#'
#' @param surv Either 1) a numeric vector describing a survival curve (lx), or
#'   2) a \code{data.frame} / \code{list} with one column / element titled 'lx'
#'   describing a survival curve, optionally a column / element 'x' containing
#'   age classes (each element a number representing the age at the start of the
#'   class). 
#'
#'   If x is not supplied, the function will assume age classes starting at 0
#'   with time steps of 1 unit. If x begins at 0 then \code{lx[1]} should equal
#'   1. If x ends at maximum longevity, then \code{lx[which.max(x)]} should
#'   equal 0; however it is possible to supply partial survivorship curves.
#' @param survTable Logical to determine whether a survival table should be
#'   returned (see details).
#'
#' @return the pace of survival, equal to life expectancy at birth. If 
#'   \code{survTable = TRUE}, a list containing the pace of survival, and 
#'   a survival table which contains the parameters required to calculate pace:
#'   Lx (midpoint survivorshop), Tx (total number of age categories left for
#'   individuals surviving to beginning of age interval), ex (life expectancy
#'   at the start of the age category). Note that life expectancy at other ages 
#'   (especially at age of maturity) may also be considered as measures of 
#'   the pace of survival.
#' 
#' @author Iain Stott <iainmstott@@gmail.com>
#' 
#' @examples
#' # exponential decline in lx yields shape = 0
#' lx <- 0.7^(0:20)
#' pace_surv(lx)
#' 
#' @export pace_surv
pace_surv <- function(surv, survTable = FALSE) {
  if(class(surv) %in% "numeric") {
    lx <- surv
    x <- seq_along(lx) - 1
  }
  if(class(surv) %in% c("list", "data.frame")) {
    if(!all(c("x", "lx") %in% names(surv))) {
      stop("'surv' is a data.frame or list and doesn't contain both x and lx")
    }
    x <- surv$x
    lx <- surv$lx
    if(length(x) != length(lx)) {
      stop("x and lx must be the same length")
    }
  }
  xmin <- min(x)
  xmax <- max(x)
  if (!(xmin %in% 0)) stop("x must begin with 0")
  if (!(lx[1] %in% 1) | !(lx[xmax + 1] %in% 0)) stop("lx must start with 1 and end with 0")
  if(any(diff(x) <= 0)) stop("much as we'd like to reverse aging, x must all be ascending")
  if(any(diff(lx) > 1e-7)) stop("please don't bring organisms back from the dead (check lx)")
  if(length(x) <= 2) {
    stop("must have > 2 values of lx to calculate pace")
  }
  Lx <- (lx + lx[seq(2, xmax + 2, 1)]) / 2
  Lx[xmax + 1] <- 0
  Tx <- numeric(xmax + 1)
  Tx[1] <- sum(Lx)
  i <- 1
  while (i < (xmax + 1)) {
    Tx[i + 1] <- Tx[i] - Lx[i]
    i <- i + 1
  }
  Tx[xmax + 1] <- 0
  ex <- Tx / lx
  ex[xmax + 1] <- NA
  pace <- ex[1]
  if (!survTable) return(pace)
  if (survTable) {
    survTable <- data.frame(x = x,
                            lx = lx,
                            Lx = Lx,
                            Tx = Tx,
                            ex = ex)
    return(list(pace = pace, survTable = survTable))
  }
}
