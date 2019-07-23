#' Calculate shape of reproduction over age
#'
#' Calculates a 'shape' value of distribution of reproduction over age by
#' comparing the area under a cumulative reproduction curve (over age) with the
#' area under a cumulative function describing constant reproduction.
#'
#' @param rep Either 1) a numeric vector describing reproduction over age (mx),
#'   or 2) a \code{data.frame} / \code{list} with one column / element titled
#'   'mx' describing a reproduction over age, optionally a column / element 'x'
#'   containing age classes (each element a number representing the age at the
#'   start of the class).
#'   
#'   If x is not supplied, the function will assume age classes starting at 0
#'   with time steps of unit. If x ends at maximum longevity,
#'   \code{mx[which.max(x)]} should equal 0; however it is possible to supply
#'   partial reproduction schedules.
#' @param xmin,xmax The minimum and maximum age repectively over which to
#'   evaluate shape. If not given, these default to \code{min(x)} and
#'   \code{max(x)} respectively.
#'
#' @return a shape value describing symmetry of reproduction over age by
#'   comparing the area under a cumulative reproduction curve over age with the
#'   area under constant reproduction. May take any real value between -0.5 and
#'   +0.5. A value of 0 indicates negligible aging (neither generally increasing
#'   nor generally decreasing reproduction with age); positive values indicate
#'   senescence (generally decreasing reproduction with age); negative values
#'   indicate negative senescence (generally increasing reproduction with age).
#'   A value of +0.5 indicates that (hypothetically) all individuals are born to
#'   individuals of age 0; a value of -0.5 indicates that all individuals are
#'   born at the age of maximum longevity.
#' 
#' @author Iain Stott <iainmstott@@gmail.com>
#' 
#' @examples
#' # increasing mx yields negative shape
#' mx <- c(0, 0, 0.3, 0.4, 0.5, 0.6)
#' shape_rep(mx)
#' 
#' # decreasing mx yields positive shape
#' mx <- c(1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4)
#' shape_rep(mx)
#' 
#' # constant mx yields shape = 0
#' mx <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' shape_rep(mx)
#'
#' @export shape_rep
shape_rep <- function(rep, xmin = NULL, xmax = NULL, fertTable = FALSE) {
  if(class(rep) %in% "numeric") {
    mx <- rep
    x <- seq_along(mx) - 1
  }
  if(class(rep) %in% c("list", "data.frame")) {
    if(!all(c("x", "mx") %in% names(rep))) {
      stop("'rep' is a data.frame or list and doesn't contain both x and mx")
    }
    x <- rep$x
    mx <- rep$mx
    if(length(x) != length(mx)) {
      stop("x and mx must be the same length")
    }
  }
  if(any(duplicated(x))) stop("all x must be unique values")
  if(any(diff(x) <= 0)) stop("much as we'd like to reverse aging, x must all be ascending")
  if(any(mx[!is.na(mx)] < 0)) stop("You appear to have minus-babies (check mx)")
  if(any(length(xmin) > 1, length(xmax) > 1)){
    stop("xmin and xmax must have length 1 or NULL")
  }
  ltdim <- length(x) 
  if(is.null(xmin)) xmin_fix <- x[min(which(mx > 0))]
  if(!is.null(xmin)) xmin_fix <- xmin
  if(is.null(xmax)) { 
    if(is.na(mx[ltdim])) {
      x_fix <- x
      ltdim_fix <- length(x_fix)
      xmax_fix <- max(x_fix)
      mx_fix <- mx
      x_sub <- x_fix[x_fix >= xmin_fix & x_fix <= xmax_fix]
      ltdim_sub <- length(x_sub)
      mx_sub <- mx_fix[x >= xmin_fix & x <= xmax_fix]
    }
    if(!is.na(mx[ltdim])) {
      x_fix <- c(x, x[ltdim] + (x[ltdim] - x[ltdim - 1]))
      ltdim_fix <- length(x_fix)
      xmax_fix <-  max(x_fix)
      mx_fix <- c(mx, NA)
      x_sub <- x_fix[x_fix >= xmin_fix & x_fix <= xmax_fix]
      ltdim_sub <- length(x_sub)
      mx_sub <- mx_fix[x_fix >= xmin_fix & x_fix <= xmax_fix]
    }
  }
  if(!is.null(xmax)){ 
    if(is.na(mx[which(x == xmax)])){
      x_fix <- x
      ltdim_fix <- length(x_fix)
      xmax_fix <- xmax
      mx_fix <- mx
      x_sub <- x_fix[x_fix >= xmin_fix & x_fix <= xmax_fix]
      ltdim_sub <- length(x_sub)
      mx_sub <- mx_fix[x_fix >= xmin_fix & x_fix <= xmax_fix]
    }
    if(!is.na(mx[which(x == xmax)])){
      x_fix <- x
      ltdim_fix <- length(x_fix)
      xmax_fix <- xmax
      mx_fix <- mx
      x_sub <- x_fix[x_fix >= xmin_fix & x_fix <= xmax_fix]
      ltdim_sub <- length(x_sub)
      mx_sub <- mx_fix[x >= xmin_fix & x <= xmax_fix]
      mx_sub[xmax_fix] <- NA
    }
  }
  if(ltdim_sub <= 2 ) {
    stop("must have > 2 values of mx to calculate shape")
  }
  lt_sub_int <- diff(x_sub)
  Bx_sub <- c(0, cumsum(mx_sub[seq(1, ltdim_sub-1, 1)]) * lt_sub_int)
  B <- max(Bx_sub)
  x_std <- (x_sub - xmin_fix) / (xmax_fix - xmin_fix)
  # standardised mx has mean of 1
  # last "fix" class is NA as no more offspring past the start of the class
  mx_std <- (mx_sub / B) * (xmax_fix - xmin_fix)
  Bxmin <- Bx_sub[which.min(x_std)]
  Bxmax <- Bx_sub[which.max(x_std)]
  Bx_std <- (Bx_sub - Bxmin) / (Bxmax - Bxmin) 
  auc_std <- area_under_curve(x_std, Bx_std)
  auc_flat <- 0.5
  shape <- auc_std - auc_flat
  if(!fertTable) return(shape)
  if(fertTable) { 
    fertTable <- data.frame(x = x_sub,
                            mx = mx_sub,
                            Bx = Bx_sub,
                            xStd = x_std,
                            mxStd = mx_std,
                            BxStd = Bx_std)
    return(list(shape = shape, fertTable = fertTable))
  }
}
