#' Filter against the sample_covariate_tbl as if it were wide.
#'
#' This feature is implemented so poorly. It's only really meant to be used
#' interactively, and with extreme caution ... programatically specifying
#' column names in feature table, for instance, does not work right now.
#' @export
#' @param x A \code{FacileDataSet}
#' @param ... NSE claused to use in \code{\link[dplyr]{filter}} expressions
#' @family API
filter_features <- function(x, ...) {
  feature_info_tbl(x) %>% filter(...) %>% set_fds(x)
}
