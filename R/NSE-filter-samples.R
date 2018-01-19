#' Filter against the sample_covariate_tbl EAV table as if it were wide.
#'
#' This allows the user to query the `FacileDataSet` as if it were a wide
#' `pData` `data.frame` of all its covariates.
#'
#' This feature is implemented so poorly. It's only really meant to be used
#' interactively, and with extreme caution ... programatically specifcying
#' the covariates, for instance, does not work right now.
#'
#' TODO: Professionaly implement this. You will minimally want to use `tidyeval`
#'
#' @md
#' @export
#' @importFrom lazyeval lazy_dots auto_name
#'
#' @param x A \code{FacileDataSet}
#' @param ... NSE claused to use in \code{\link[dplyr]{filter}} expressions
#' @return a sample-descriptor `data.frame` that includes the dataset,sample_id
#'   pairs that match the virtual `filter(covaries, ...)` clause executed here.
#'
#' @examples
#' fds <- exampleFacileDataSet()
#'
#' # To identify all samples that are of "CMS3" or "CMS4" subtype(
#' # stored in the "subtype_crc_cms" covariate:
#' crc.34 <- filter_samples(fds, subtype_crc_cms %in% c("CMS3", "CMS4"))
#' eav.query <- sample_covariate_tbl(fds) %>%
#'   filter(variable == "subtype_crc_cms", value %in% c("CMS3", "CMS4")) %>%
#'   collect
#' setequal(crc.34$sample_id, eav.query$sample_id)
filter_samples <- function(x, ..., with_covariates=FALSE) {
  stopifnot(is.FacileDataSet(x))
  dots <- lazy_dots(...)
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

  dot.exprs <- names(auto_name(dots))
  hits <- sapply(all.vars, function(var) any(grepl(var, dot.exprs)))
  out <- names(hits)[hits]
  if (length(out) == 0) {
    stop("No sample covariates found in query: ",
         paste(dot.exprs, collapse=';'))
  }
  out
}
