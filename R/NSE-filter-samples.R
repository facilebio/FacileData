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
#' eav.query <- fds |>
#'   fetch_sample_covariates(covariates = "subtype_crc_cms") |>
#'   filter(value %in% c("CMS3", "CMS4")) |>
#'   collect()
#' setequal(crc.34$sample_id, eav.query$sample_id)
#'
#' # You can keep filtering a filtered dataset
#' crc.34.male <- filter_samples(crc.34, sex == "m")
filter_samples.FacileDataSet <- function(x, ..., samples. = samples(x),
                                         custom_key = Sys.getenv("USER"),
                                         with_covariates = FALSE) {
  # cov.table <- .create_wide_covariate_table(x, dots)
  # out <- dplyr::filter_(cov.table, .dots=dots)

  force(samples.)
  assert_sample_subset(samples.)

  cov.table <- .create_wide_covariate_table(x, samples., ...,
                                            custom_key = custom_key)
  out <- filter(cov.table, ...)
  if (!with_covariates) {
    out <- select(out, dataset, sample_id)
  }
  if (nrow(out) == 0L) {
    warning("All samples have been filtered out", immediate. = TRUE)
  }
  as_facile_frame(out, x)
}

#' @noRd
#' @export
filter_samples.facile_frame <- function(x, ...,
                                        custom_key = Sys.getenv("USER"),
                                        with_covariates = FALSE) {
  .fds <- assert_facile_data_store(fds(x))
  assert_sample_subset(x)
  filter_samples(.fds, ..., samples. = x, custom_key = custom_key,
                 with_covariates = with_covariates)
}

#' @noRd
#' @importFrom lazyeval lazy_dots
.create_wide_covariate_table <- function(x, samples, ...,
                                         custom_key = Sys.getenv("USER")) {
  assert_facile_data_store(x)
  assert_sample_subset(samples)

  out <- fetch_sample_covariates(x, samples = samples, custom_key = custom_key)
  dots <- lazy_dots(...)
  qvars <- .parse_filter_vars(x, dots)

  # TODO: check if any of the query variables are dataset or sample_id, then
  # fiter `out` on the dataset or sample_id columns, THEN play with the
  # other sample covariates (sc)
  pk.vars <- intersect(qvars, c("dataset", "sample_id"))
  # if (length(pk.vars)) {
  #   out <- filter(out, pk.part.of.query)
  # }

  sc.vars <- setdiff(qvars, c("dataset", "sample_id"))
  if (length(sc.vars)) {
    out <- filter(out, variable %in% !!qvars)
  }
  out |>
    spread_covariates() |>
    distinct(dataset, sample_id, .keep_all = TRUE)
}

#' @noRd
#' @importFrom lazyeval auto_name
.parse_filter_vars <- function(x, dots) {
  assert_facile_data_store(x)
  stopifnot(is(dots, 'lazy_dots'))

  all.vars <- sample_covariate_tbl(x) |>
    distinct(variable) |>
    collect(n=Inf)
  all.vars <- c(all.vars$variable, "dataset", "sample_id")

  dot.exprs <- names(auto_name(dots))
  hits <- sapply(all.vars, function(var) any(grepl(var, dot.exprs)))
  out <- names(hits)[hits]
  if (length(out) == 0) {
    stop("No sample covariates found in query: ",
         paste(dot.exprs, collapse=';'))
  }
  out
}
