#' Mean and Confidence Interval
#'  
#' @param x numeric vector
#' @param level confidence level
#' 
#' @importFrom stats qt sd
#'  
#' @note Internal function. Not supposed to be called directly
#'
#' @examples 
#'      stopifnot(all.equal(ci_mean(NA_real_), c(NA_real_, NA_real_, NA_real_)))
#'      stopifnot(all.equal(ci_mean(1), c(lower = NA_real_, mean = 1, upper = NA_real_)))
#'      stopifnot(all.equal(
#'          ci_mean(1:9, 0.95),
#'          c(lower = 2.934942, mean = 5.000000, upper = 7.065058), 
#'          tolerance = 0.0001)
#'      )
#'      
#' @export
ci_mean <- function(x, level = 0.90) {
    if (all(is.na(x))) {
        return(c(NA_real_, NA_real_, NA_real_))
    }
    n <- sum(!is.na(x))
    m <- mean(x, na.rm = TRUE)
    s <-   sd(x, na.rm = TRUE) / sqrt(n)
    t <- qt(p = 0.5 * (1 + level), df = n)
    c(lower = m - t * s, mean = m, upper = m + t * s)
}


#' Construct a Text Representation of a Weight Vector
#'  
#' @param x numeric or character vector
#'  
#' @note Internal function. Not supposed to be called directly
#'
#' @export
toString_weights <- function(x) {
    if (is.null(x)) {
        return("optimized")
    } else {
        if(attr(x, "normalized")) {
             return(paste(toString(formatC(x, format = "f", digits = 3)), "(normalized and fixed)"))
        } else {
            if (is.null(attr(x, "character"))) {
                return(paste(toString(formatC(x, format = "f", digits = 3)), "(fixed)"))
            } else {
                return(paste(toString(attr(x, "character")), "(fixed)"))
            }
        }
    }
}