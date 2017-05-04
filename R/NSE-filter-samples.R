##' Filter against the sample_covariate_tbl as if it were wide.
##'
##' This feature is implemented so poorly. It's only really meant to be used
##' interactively, and with extreme caution ... programatically specifcying
##' the covariates, for instance, does not work right now.
##'
##' TODO: Professionaly implement this.
##'
##' @export
##' @param x A \code{FacileDataSet}
##' @param
filter_samples <- function(x, ..., with_covariates=FALSE) {
  stopifnot(is.FacileDataSet(x))
  dots <- lazyeval::lazy_dots(...)
  cov.table <- .create_wide_covariate_table(x, dots)
  out <- dplyr::filter_(cov.table, .dots=dots)
  if (!with_covariates) {
    out <- select(out, dataset, sample_id)
  }
  set_fds(out, x)
}

.create_wide_covariate_table <- function(x, dots) {
  stopifnot(is.FacileDataSet(x))
  sc <- sample_covariate_tbl(x)
  qvars <- .parse_filter_vars(x, dots)
  fetch_sample_covariates(x, covariates=qvars) %>% spread_covariates
}

.parse_filter_vars <- function(x, dots) {
  stopifnot(is.FacileDataSet(x))
  stopifnot(is(dots, 'lazy_dots'))

  all.vars <- sample_covariate_tbl(x) %>%
    distinct(variable) %>%
    collect(n=Inf)
  all.vars <- all.vars$variable

  dot.exprs <- names(lazyeval::auto_name(dots))
  hits <- sapply(all.vars, function(var) any(grepl(var, dot.exprs)))
  out <- names(hits)[hits]
  if (length(out) == 0) {
    stop("No sample covariates found in query: ",
         paste(dot.exprs, collapse=';'))
  }
  out
}
