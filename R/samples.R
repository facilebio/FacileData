
##' Creates tbl that associates a sample with all of its indication/subtypes
##'
##' @export
##' @param db a \code{FacileDb} objectg
##' @return a \code{tbl} with indication and subtype information for all samples
##'   in the database
subtype_map <- function(db) {
  if (missing(db)) on.exit(dbDisconnect(db$con))
  sample.map <- sample_covariate_tbl(db) %>%
    filter(class == 'tumor_classification') %>%
    with_sample_covariates('sample_type', db, db$cov.def)

  main <- sample.map %>%
    filter(variable == 'indication') %>%
    transmute(indication=value, subtype='all', dataset=dataset,
              sample_id=sample_id, sample_type=sample_type)

  subs <- filter(sample.map, variable != 'indication') %>%
    transmute(subtype=value, dataset=dataset, sample_id=sample_id,
              sample_type=sample_type) %>% ## add indication
    left_join(select(main, dataset, sample_id, indication),
              by=c('dataset', 'sample_id'))

  bind_rows(main, subs) %>%
    arrange(indication, subtype, dataset) %>%
    set_fdb(db)
}

##' Fetches a sample descriptor that matches filter criterion over covariates.
##'
##' @export
##' @param db A \code{FacileDb}
##' @param ... filter clause to apply to \code{sample_covariate_tbl(db)}
##' @return a \code{tbl} sample-descriptor with a \code{FacileDb} connection
##' @examples
##' fetch_samples(db, indication %in% c('BRCA', 'COAD'))
fetch_samples <- function(db, ...) {
  sample_covariate_tbl(db) %>%
    filter(...) %>%
    select(dataset, sample_id) %>%
    set_fdb(db)
}

##' Filters the samples down in a dataset to ones specified
##'
##' Tables like \code{expression} and \code{sample_covariate} house different
##' datapoints per sample, and we often want to only retreive data points over
##' a subset of samples.
##'
##' @export
##' @param x likely a \code{tbl_sqlite} object, but a \code{tbl_df}-like
##'   object should work as well.
##' @param samples a sample descriptor \code{tbl_df}-like object (likely a
##'   \code{tbl_sqlite} object) that has \code{"dataset"} and \code{"samle_id"}
##'   columns.
##' @param semi if \code{TRUE}, appropximates a semi-join on the \code{samples},
##'   otherwise does an inner_join between \code{x} and \code{samples}
##'   (default \code{FALSE}).
##' @return joined result between \code{x} and \code{samples}
join_samples <- function(x, samples=NULL, semi=FALSE, distinct.samples=FALSE) {
  if (is.null(samples)) {
    return(x)
  }
  assert_sample_subset(samples)

  ## TODO: Rethink the internalization choice here. If the "external" dataset
  ##       is huge, then copying it into the database to do the join will be
  ##       painful
  internalize <- !same_src(samples, x)

  ## I think I should be using `semi_join` here, but that is so slow I might
  ## as well be looking things up by hand
  extra.cols <- setdiff(colnames(samples), c('dataset', 'sample_id'))
  if (semi && length(extra.cols) > 0L) {
    samples <- select(samples, dataset, sample_id)
  }

  inner_join(x, samples, by=c('dataset', 'sample_id'),
             copy=internalize, auto_index=internalize) %>%
    set_fdb(fdb(x))
}

##' Filter x down to specific samples
##'
##' @export
##' @param x something like a \code{tbl_sqlite} object
##' @param samples a sample descriptor
##' @return filtered version of \code{x} that only has the desired samples
filter_samples <- function(x, samples=NULL) {
 join_samples(x, samples, semi=TRUE)
}
