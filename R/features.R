#' Enumerate the types of feature stored in a FacileDataSet
#'
#' @export
#' @param x A \code{FacileDataSet}
feature_types <- function(x) {
  stopifnot(is.FacileDataSet(x))
  ## Damn, can't do distinct on sqlite
  feature_info_tbl(x) %>%
    distinct(feature_type) %>%
    collect(n=Inf) %$%
    feature_type
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
