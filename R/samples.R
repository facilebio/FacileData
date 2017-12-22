
##' Creates tbl that associates a sample with all of its indication/subtypes
##'
##' @export
##' @param x a \code{FacileDataSet} object
##' @return a \code{tbl} with indication and subtype information for all samples
##'   in the database
subtype_map <- function(x) {
  stopifnot(is.FacileDataSet(x))
  sample.map <- sample_covariate_tbl(x) %>%
    filter(type == 'tumor_classification') %>%
    with_sample_covariates('sample_type', .fds=x)

  main <- sample.map %>%
    filter(variable == 'indication') %>%
    transmute(indication=value, subtype='all', dataset=dataset,
              sample_id=sample_id, sample_type=sample_type)

  subs <- sample.map %>%
    filter(variable != 'indication') %>%
    transmute(subtype=value, dataset=dataset, sample_id=sample_id,
              sample_type=sample_type) %>% ## add indication
    left_join(select(main, dataset, sample_id, indication),
              by=c('dataset', 'sample_id'))

  bind_rows(main, subs) %>%
    arrange(indication, subtype, dataset) %>%
    set_fds(x)
}

# ##' Fetches a sample descriptor that matches filter criterion over covariates.
# ##'
# ##' @export
# ##' @param x A \code{FacileDataSet} object
# ##' @param ... filter clause to apply to \code{sample_covariate_tbl(x)}
# ##' @return a facile sample descriptor with a \code{FacileDataSet} connection.
# ##' @examples
# ##' exampleFacileDataSet() %>%
# ##'   fetch_samples(indication %in% c('BRCA', 'COAD'))
# fetch_samples <- function(x, ...) {
#   stopifnot(is.FacileDataSet(x))
#   sample_covariate_tbl(x) %>%
#     filter(...) %>%
#     collect(n=Inf) %>%
#     distinct(dataset, sample_id) %>%
#     set_fds(x)
# }

##' Fetches a sample descriptor that matches the filter criterion.
##'
##' Use \code{...} as if this is a dplyr::filter call, and our
##' sample_covariate_tbl was "wide".
##'
##' This is experimental, so each "term" in the filter criteria should be
##' just one boolean operation. Multiple terms passed into \code{...} will be
##' "AND"ed together.
##'
##' @export
##' @importFrom lazyeval lazy_dots
##' @param x A \code{FacileDataRepository}
##' @param ... the NSE boolean filter criteria
##' @return a facile sample descriptor
fetch_samples <- function(x, samples=NULL, assay="rnaseq", ...) {
  stopifnot(is.FacileDataSet(x))
  dots <- lazy_dots(...)
  if (length(dots)) {
    stop("Currently rethinking how to make fetching samples intuitive, ie. ",
         "see fetch_samples.old")
  }
  if (is.null(samples)) samples <- sample_info_tbl(x)
  samples <- assert_sample_subset(samples) %>% collect(n=Inf)

  if (!missing(assay)) {
    assert_string(assay)
    if (is.character(samples$assay)) {
      warning("assay specified in parameter, but already exists in sample ",
              "descriptor. The paremter will override value in sample ",
              "descriptor", immediate.=TRUE)
    }
    samples$assay <- assay
    fds.tbl <- 'assay_sample_info'
  } else {
    fds.tbl <- 'sample_info'
  }

  copy <- !is(samples, 'tbl_dbi')
  pk <- primary_key(x, fds.tbl)

  tbl(x$con, fds.tbl) %>%
    semi_join(samples, by=pk, copy=copy, auto_index=copy) %>%
    set_fds(x)
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
    samples <- distinct(samples, dataset, sample_id)
  }

  inner_join(x, samples, by=c('dataset', 'sample_id'),
             copy=internalize, auto_index=internalize) %>%
    set_fds(fds(x))
}

## Filter x down to specific samples
##
## @export
## @param x something like a \code{tbl_sqlite} object
## @param samples a sample descriptor
## @return filtered version of \code{x} that only has the desired samples
# filter_samples <- function(x, samples=NULL) {
#  join_samples(x, samples, semi=TRUE)
# }

## FacileExplorer's filter_active_samples ======================================

##' @export
retrieve_samples_in_memory <- function(criteria, cov.table=NULL) {
  if(length(criteria)==0){
    # case where no filters have been defined
    indiv.results <- list(cov.table %>%
                            select(dataset, sample_id) %>%
                            collect)
  } else {
    indiv.results <- lapply(criteria, function(crit) {
      if (!is.null(crit)) {
        dots <- parse_sample_criterion(variable=crit$variable, value=crit$value)
        cov.table %>%
          filter_(.dots=dots) %>%
          select(dataset, sample_id) %>%
          collect
      }
    })
  }

  ## Now we go back and take the intersection of the dataset,sample_id pairs
  ## found in each element of the list above
  if (length(indiv.results) == 1L) {
    out <- indiv.results[[1]]
  } else {
    rf <- function(x,y) semi_join(x, y, by=c('dataset', 'sample_id'))
    out <- Reduce(rf, indiv.results[-1], init=indiv.results[[1]])
  }

  out %>%
    distinct(dataset, sample_id) %>%
    arrange(dataset, sample_id)
}

##' Creates a filter expression to select samples based on value of a covariate
##'
##' This leverages dplyr's standard (vs non-standard) evaluation mojo. There is
##' likely a cleaner way to do this, but to be honest I still find the
##' \code{\link[lazyeval]{interp}} stuff rather confusing
##'
##' @seealso \href{https://cran.r-project.org/web/packages/dplyr/vignettes/nse.html}{dplyr non-standard evaluation}
##'
##' @importFrom lazyeval interp
##' @importFrom stats formula
##' @param variable the name of the variable to look for in the sample_covariate
##'   \code{variable} column
##' @param value \code{character} vector of values for the \code{variable} that
##'   you want your samples to have.
##' @return a
parse_sample_criterion <- function(variable, value) {
  stopifnot(is.character(variable) && length(variable) == 1L)
  stopifnot(is.character(value) && length(value) >= 1)

  vals <- paste(sprintf("'%s'", value), collapse=',')
  if (length(value) == 1L) {
    crit <- paste0("~ variable == dbvar & value == ", vals)
  } else if (length(value) > 1) {
    crit <- paste0('~ variable == dbvar & value %in% c(', vals, ')')
  } else {
    stop("length(value) <= 0")
  }

  dots <- interp(formula(crit), dbvar=variable)
  dots
}
