################################################################################
# Make a type called cSurv that is a character representation of survival::surv
################################################################################

#' @name cSurv
#' @title cSurv is a character representation of survival::surv ###
#'
#' @description cSurv serves as a more reliable way to use Surv objects as data.frame columns. A
#' data.frame is supposed to be able to hold Surv columns. There are multiple special
#' cases written into base for this. It seems the implementation is incomplete as
#' subsetting the DF breaks the Surv object. cSurv cannot do anything but get subset
#' and become a Surv again. In the FacileVerse we hold Surv objects as cSurv, which
#' allows us to survive a round-trip through an EAV sample metadata table. Survival
#' analyses can convert cSurv to Surv as needed. It is assumed that all Surv censoring
#' is right-censored.
#' @importFrom survival Surv
#' @examples
#' library(survival)
#' x = Surv(c(14,12,3), event = c(1,0,1))
#' y = as(x,"cSurv")
#' z = as(y, "Surv")
#' x2 = as.character(x)
#' z2 = as(x2, "Surv")
#' a = as(x, "cSurv")
#' b = as(a, "character")
#' c = as(b, "cSurv")
#' d = as(c, "Surv")
NULL

setOldClass("Surv")
setOldClass("cSurv")

#' @rdname cSurv
#' @family cSurv
#' @export
as_cSurv <- function(from) {
  structure(as.character(from), class = "cSurv")
}

#' @rdname cSurv
#' @family cSurv
#' @export
as_Surv <- function(from) {
  ns <- tryCatch(loadNamespace("survival"), error = function(e) NULL)
  if (is.null(ns)) stop("survival package required")

  from <- as.character(from) # Both check type and drop attributes
  stopifnot(all(is.na(from) | grepl("\\d[\\+ ]*$", from)))
  status <- ifelse(endsWith(from, "+"), 0, 1)
  ns$Surv(as.numeric(gsub("[\\+ ]*$", "", from)), status)
}

#' @family cSurv
setAs(
  from = "Surv",
  to = "cSurv",
  def = as_cSurv
)

#' @family cSurv
setAs(
  from = "cSurv",
  to = "Surv",
  def = as_Surv
)

#' @family cSurv
setAs(
  from = "character",
  to = "Surv",
  def = as_Surv
)

#' @family cSurv
setAs(
  from = "character",
  to = "cSurv",
  def = function(from) {
    as_cSurv(as_Surv(from))
  }
)
