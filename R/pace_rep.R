
#' Calculate shape of reproduction over age
#'
#' Calculates a 'shape' value of distribution of reproduction over age by comparing
#' the area under a cumulative reproduction curve (over age) with the area 
#' under a cumulative function describing constant reproduction.
#'
#' @param rep Either 1) a numeric vector describing reproduction over age (mx); 2) a 
#'   \code{data.frame} / \code{list} with one column / element titled 'mx' 
#'   describing a reproduction over age, optionally a column / element 'x' containing 
#'   age classes (each element a number representing the age at the start of the
#'   class); 3) a list containing two elements: 'matU', a U matrix (the survival 
#'   component of a matrix population model, i.e. a square projection matrix reflecting 
#'   survival-related transitions, e.g. progression, stasis, and retrogression) and
#'   'matF', and F matrix(the reproduction component of a matrix projection model, 
#'   i.e. a square projection matrix of the same dimension as matU reflecting 
#'   those transitions describing sexual reproduction); 4) a \code{CompadreMat} 
#'   object (\code{RCompadre-package})containing a matrix population model in 
#'   the format described in the \code{CompadreMat} class.
#'   In the case of 1 and 2 where x is not supplied, the function will assume
#'   age classes starting at 0 with steps of 1 year.
#'   In the case of 3 and 4, an age-based reproduction schedule will be generated from 
#'   a stage-based matrix using the \code{makeLifeTable} function of 
#'   \code{RCompadre-package}. 
#'   In all cases where x ends at maximum longevity, \code{mx[which.max(x)]} 
#'   should equal 0, however it is possible to supply partial reproduction schedules.
#' @param xmin,xmax The minimum and maximum age repectively over which to evaluate
#'   shape. If not given, these default to \code{min(x)} and \code{max(x)} 
#'   respectively.
#' @param fertTable logical determining whether to return the fertility table 
#'   used to calculate pace, including standardised measures of age, fertility
#'   and standardised fertility.
#'
#' @return the pace of reproduction, equal to the age of the mother at which an
#'   average baby is born. This is the birth equivalent to life expectancy, 
#'   calculated using the cumulative reproduction.
#' 
#' @author Iain Stott <iainmstott@@gmail.com>
#' 
#' @examples
#' mx <- c(0, 0, 0.3, 0.4, 0.5, 0.6)
#' pace_rep(mx)
#'
#' @export pace_rep
pace_rep <- function(rep, xmin = NULL, xmax = NULL, fertTable = FALSE){
    if (class(rep) %in% "numeric") {
        mx <- rep
        x <- seq_along(mx) - 1
    }
    if (class(rep) %in% c("list", "data.frame")) {
        if (!all(c("x", "mx") %in% names(rep))) {
            stop("'rep' is a data.frame or list and doesn't contain both x and mx")
        }
        x <- rep$x
        mx <- rep$mx
        if (length(x) != length(mx)) {
            stop("x and mx must be the same length")
        }
    }
    shape <- shape_rep(rep = data.frame(x = x, mx = mx), 
                       xmin = xmin, xmax = xmax, fertTable = TRUE)
    shape_xmin <- min(shape$fertTable$x)
    shape_xmax <- max(shape$fertTable$x)
    pace <- (shape_xmax - shape_xmin) * (0.5 - shape$shape)
    if (!fertTable) return(pace)
    if (fertTable) return(list(pace = pace, fertTable = shape$fertTable))
}
