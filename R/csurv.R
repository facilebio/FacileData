#####################################################################################
### Make a type called cSurv that is a character representation of survival::surv ###
### Serves as a more reliable way to use Surv objects as data.frame columns       ###
#####################################################################################

#' @family cSurv
eav_encode_Surv <- function(x, ...) {
    stopifnot(is(x, "Surv"))
    out <- as.character(x)
    class(out) <- "cSurv"
    out
}

#' @family cSurv
eav_decode_Surv <- function(x, attrname = character(), def = list(), ...) {
    x = as.character(x) # Both check type and drop attributes
    stopifnot(all(is.na(x) | grepl("\\d[\\+ ]*$", x)))
    status = ifelse(endsWith(x,"+"), 0, 1)
    Surv(as.numeric(gsub("[\\+ ]*$", "", x)), status)
}


setAs(
    from = "cSurv",
    to = "Surv",
    def = function(x) {

    })
