#' @export
#' @noRd
fetch_feature_info.FacileDataSet <- function(x, feature_type,
                                             feature_ids = NULL, ...) {
  ftype <- assert_choice(feature_type, feature_types(x))
  out <- filter(feature_info_tbl(x), feature_type == ftype)
  if (!is.null(feature_ids)) {
    assert_character(feature_ids, min.len = 1)
    if (length(feature_ids) == 1L) {
      out <- filter(out, feature_id == feature_ids)
    } else {
      out <- filter(out, feature_id %in% feature_ids)
    }
  }
  out
}

#' @noRd
#' @export
with_feature_info.facile_frame <- function(x, covariates = NULL, ...,
                                           .fds = fds(x)) {
  x <- collect(x, n = Inf)
  .fds <- assert_facile_data_store(.fds)
  NextMethod(x, .fds = .fds)
}

#' @noRd
#' @export
with_feature_info.tbl <- function(x, covariates = NULL, ..., .fds = NULL) {
  with_feature_info.data.frame(collect(x, n = Inf), covariates = covariates,
                               ..., .fds = .fds)
}

#' NOTE: Would be nice to do with_feature_info(x, symbol = name)
#'
#' @noRd
#' @export
#' @method with_feature_info data.frame
#' @examples
#' efds <- exampleFacileDataSet()
#' some <- fetch_feature_info(efds, "entrez") %>%
#'   select(feature_id, feature_type) %>%
#'   collect(n = 5)
#' with_feature_info(some, c("name", "meta"))
#' with_feature_info(some, c(symbol = "name", biotype = "meta"))
with_feature_info.data.frame <- function(x, covariates = NULL, ...,
                                         .fds = NULL) {
  fkey <- c("feature_id", "feature_type")
  .fds <- assert_facile_data_store(.fds)
  assert_true(has_columns(x, fkey))

  # Possible things to ask for
  fattribs <- colnames(fetch_feature_info(.fds, feature_types(.fds)[1L]))
  fattribs <- setdiff(fattribs, fkey)

  if (is.null(covariates)) covariates <- fattribs
  # unique removes
  covariates <- assert_subset(covariates[!duplicated(covariates)], fattribs)
  covariates <- nameit(covariates)

  new_info <- lapply(unique(x[["feature_type"]]), function(ftype) {
    # I can do this lazily, but ...
    fi <- fetch_feature_info(.fds, ftype)
    fi <- select(fi, !!fkey, !!covariates)
    collect(fi, n = Inf)
  })
  new_info <- bind_rows(new_info)

  out <- left_join(x, new_info, by = fkey, suffix = c("", ".infds"))
  as_facile_frame(out, .fds)
}

#' Enumerate the types of feature stored in a FacileDataSet
#'
#' @export
#' @param x A \code{FacileDataSet}
feature_types <- function(x) {
  stopifnot(is.FacileDataSet(x))
  ## Damn, can't do distinct on sqlite
  feature_info_tbl(x) %>%
    distinct(feature_type) %>%
    collect(n=Inf) %>%
    pull(feature_type)
}

#' Test if a given feature type is stored in a FacileDataSet
#'
#' @export
#' @param x A \code{FacileDataSet}
#' @param feature_type a character vector of potential feature types
#' @return logical vector indicating whether or not a given \code{feature_type}
#'   is stored in \code{x}
has_feature_type <- function(x, feature_type) {
  stopifnot(is.FacileDataSet(x))
  assert_character(feature_type)
  feature_type %in% feature_types(x)
}

#' Returns table of names and aliases for features.
#'
#' #dropme
#'
#' @export
#' @param x \code{FacileDataSet}
#' @param feature_type a character vector specifying the feature type
#' @return a tibble with \code{feature_id, name, type} columns, where type
#'   is "primary" or "alias"
feature_name_map <- function(x, feature_type) {
  requireNamespace("AnnotationDbi") || stop("Failed to require AnnotationDbi.")
  stopifnot(has_feature_type(x, feature_type))
  ## http://jira.gene.com/jira/browse/FACILEDATA-64 will put this in database
  ## but for now we have flat files.
  finfo <- feature_info_tbl(x) %>%
    filter(feature_type == feature_type) %>%
    select(feature_id, name) %>%
    collect(n=Inf) %>%
    mutate(type='primary')
  if (organism(x) == 'Homo sapiens') {
    if (FALSE) {
      requireNamespace("org.Hs.eg.db") || stop("Failed to require org.Hs.eg.db")
      alias <- org.Hs.eg.db %>%
        AnnotationDbi::select(finfo$feature_id, c('ENTREZID', 'ALIAS')) %>%
        transmute(feature_id=ENTREZID, name=ALIAS, type='alias') %>%
        anti_join(finfo, by=c('feature_id', 'name'))
      write.csv(alias, 'inst/extdata/feature-alias-map.human.csv', row.names=FALSE)
    }
    alias <- system.file('extdata', 'feature-alias-map.human.csv', package='FacileData')
    alias <- read.csv(alias, colClasses='character')
  } else if (organism(x) == 'Mus musculus') {
      if (FALSE) {
          requireNamespace("org.Mm.eg.db") || stop("Failed to require org.Mm.eg.db")
      alias <- org.Mm.eg.db %>%
        AnnotationDbi::select(finfo$feature_id, c('ENTREZID', 'ALIAS')) %>%
        transmute(feature_id=ENTREZID, name=ALIAS, type='alias') %>%
        anti_join(finfo, by=c('feature_id', 'name')) %>%
        filter(!is.na(name))
      write.csv(alias, 'inst/extdata/feature-alias-map.mouse.csv', row.names=FALSE)
    }
    alias <- system.file('extdata', 'feature-alias-map.mouse.csv', package='FacileData')
    alias <- read.csv(alias, colClasses='character')
  } else {
    stop("Unsupported organism for now")
  }

  finfo %>%
    bind_rows(alias) %>%
    filter(!is.na(name)) %>%
    arrange(feature_id, name)
}
