##' Enumerate the types of feature stored in a FacileDataSet
##'
##' @export
##' @param x A \code{FacileDataSet}
feature_types <- function(x) {
  stopifnot(is.FacileDataSet(x))
  ## Damn, can't do distinct on sqlite
  feature_info_tbl(x) %>%
    distinct(feature_type) %>%
    collect(n=Inf) %$%
    feature_type
}

##' Test if a given feature type is stored in a FacileDataSet
##'
##' @export
##' @param x A \code{FacileDataSet}
##' @param feature_type a character vector of potential feature types
##' @return logical vector indicating whether or not a given \code{feature_type}
##'   is stored in \code{x}
has_feature_type <- function(x, feature_type) {
  stopifnot(is.FacileDataSet(x))
  assert_character(feature_type)
  feature_type %in% feature_types(x)
}
