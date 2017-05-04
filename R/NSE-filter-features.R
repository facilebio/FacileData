##' Filter against the sample_covariate_tbl as if it were wide.
##'
##' This feature is implemented so poorly. It's only really meant to be used
##' interactively, and with extreme caution ... programatically specifcying
##' column names in feature table, for instance, does not work right now.
##'
##' TODO: Professionaly implement this.
##'
##' @export
##' @param x A \code{FacileDataSet}
##' @param
filter_features <- function(x, ...) {
  feature_info_tbl(x) %>% filter(...) %>% set_fds(x)
}
