#' Filter against the sample_covariate_tbl EAV table as if it were wide.
#'
#' This allows the user to query the `FacileDataSet` as if it were a wide
#' `pData` `data.frame` of all its covariates.
#'
#' This feature is only really meant to be
#' used interactively, and with extreme caution ... programatically specifying
#' the covariates, for instance, does not work right now.
#'
#' TODO: Implement using `tidyeval`
#'
#' @export
#' @family API
#'
#' @param x A `FacileDataSet`
#' @param ... NSE claused to use in [dplyr::filter()] expressions
#' @return a sample-descriptor `data.frame` that includes the dataset,sample_id
#'   pairs that match the virtual `filter(covaries, ...)` clause executed here.
#'
#' @examples
#' fds <- exampleFacileDataSet()
#'
#' # To identify all samples that are of "CMS3" or "CMS4" subtype(
#' # stored in the "subtype_crc_cms" covariate:
#' crc.34 <- filter_samples(fds, subtype_crc_cms %in% c("CMS3", "CMS4"))
#' eav.query <- fds %>%
#'   fetch_sample_covariates(covariates = "subtype_crc_cms") %>%
#'   filter(value %in% c("CMS3", "CMS4")) %>%
#'   collect()
#' setequal(crc.34$sample_id, eav.query$sample_id)
filter_samples.FacileDataSet <- function(x, ...,
                                         custom_key = Sys.getenv("USER"),
                                         with_covariates = FALSE) {
  # cov.table <- .create_wide_covariate_table(x, dots)
  # out <- dplyr::filter_(cov.table, .dots=dots)
  cov.table <- .create_wide_covariate_table(x, ..., custom_key = custom_key)
  out <- filter(cov.table, ...)
  if (!with_covariates) {
    out <- select(out, dataset, sample_id)
  }
  as_facile_frame(out, x)
}

#' @noRd
#' @importFrom lazyeval lazy_dots
.create_wide_covariate_table <- function(x, ...,
                                         custom_key = Sys.getenv("USER")) {
  assert_facile_data_store(x)
  sc <- fetch_sample_covariates(x, custom_key = custom_key)
  dots <- lazy_dots(...)
  qvars <- .parse_filter_vars(x, dots)
  filter(sc, variable %in% !!qvars) %>% spread_covariates()
}

#' @noRd
#' @importFrom lazyeval auto_name
.parse_filter_vars <- function(x, dots) {
  assert_facile_data_store(x)
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
