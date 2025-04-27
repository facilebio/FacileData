#' @export
#' @noRd
fetch_feature_info.FacileDataStore <- function(x, feature_type,
                                               feature_ids = NULL, ...) {
  stopifnot(has_feature_type(x, feature_type))
  out <- features(x) |> 
    filter(.data$feature_type == .env$feature_type) |> 
    collect(n = Inf)
  
  if (!is.null(feature_ids)) {
    if (is.character(feature_ids)) {
      feature_ids <- dplyr::tibble(feature_id = feature_ids)
    }
    assert_multi_class(feature_ids, c("data.frame", "tbl"))
    assert_choice("feature_id", names(feature_ids))
    out <- semi_join(out, feature_ids, by = "feature_id")
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
#' some <- fetch_feature_info(efds, "entrez") |>
#'   select(feature_id, feature_type) |>
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
  assert_class(x, "FacileDataStore")
  assay_info(x) |> 
    distinct(feature_type) |> 
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
  assert_class(x, "FacileDataStore")
  assert_character(feature_type)
  is.element(feature_type, feature_types(x))
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
  .Defunct(paste(
    "your own symbol <> alias <> id mapping function,",
    "you might take a look at {babelgene}"))
  reqpkg("AnnotationDbi")
  stopifnot(has_feature_type(x, feature_type))
  finfo <- features(x) |>
    filter(.data$feature_type == .env$feature_type) |>
    select(feature_id, name) |>
    collect(n = Inf) |>
    mutate(type = "primary")
  if (organism(x) == 'Homo sapiens') {
    if (FALSE) {
      requireNamespace("org.Hs.eg.db") || stop("Failed to require org.Hs.eg.db")
      alias <- org.Hs.eg.db |>
        AnnotationDbi::select(finfo$feature_id, c('ENTREZID', 'ALIAS')) |>
        transmute(feature_id=ENTREZID, name=ALIAS, type='alias') |>
        anti_join(finfo, by=c('feature_id', 'name'))
      write.csv(alias, 'inst/extdata/feature-alias-map.human.csv', row.names=FALSE)
    }
    alias <- system.file('extdata', 'feature-alias-map.human.csv', package='FacileData')
    alias <- read.csv(alias, colClasses='character')
  } else if (organism(x) == 'Mus musculus') {
      if (FALSE) {
          requireNamespace("org.Mm.eg.db") || stop("Failed to require org.Mm.eg.db")
      alias <- org.Mm.eg.db |>
        AnnotationDbi::select(finfo$feature_id, c('ENTREZID', 'ALIAS')) |>
        transmute(feature_id=ENTREZID, name=ALIAS, type='alias') |>
        anti_join(finfo, by=c('feature_id', 'name')) |>
        filter(!is.na(name))
      write.csv(alias, 'inst/extdata/feature-alias-map.mouse.csv', row.names=FALSE)
    }
    alias <- system.file('extdata', 'feature-alias-map.mouse.csv', package='FacileData')
    alias <- read.csv(alias, colClasses='character')
  } else {
    stop("Unsupported organism for now")
  }

  finfo |>
    bind_rows(alias) |>
    filter(!is.na(name)) |>
    arrange(feature_id, name)
}

#' Creates a feature descriptor for interactive ease
#'
#' Creates a data.frame of features and assays they come from
#' 
#' @export
#' @param x FacileDataSet
#' @param features a character string of feature ids (requires assay_name)
#'   or a data.frame with feature_id column.
#' @param assay_name the assay to get the featurespace from. if this is provided,
#'   it will trump an already existing assay_name column in \code{features}
#' @return a feature descriptor with feature_id and assay_name, which can be
#'   used to absolutely find features
create_assay_feature_descriptor <- function(
    x, 
    features = NULL,
    assay_name = NULL) {
  # TODO: Refactor the code inside `fetch_assay_data` to use this.
  assert_facile_data_store(x)
  if (is.null(assay_name)) assay_name <- default_assay(x)
  if (is.factor(features)) features <- as.character(features)
  
  if (is.character(features) || is.null(features) || is(features, 'tbl_sql')) {
    assert_string(assay_name)
    assert_choice(assay_name, assay_names(x))
  }
  
  ainfo <- assay_info(x, assay_name)
  
  if (is.null(features)) {
    features <- collect(features(x, assay_name), n = Inf)
  } else if (is.character(features)) {
    features <- tibble(
      feature_id = features, 
      assay = assay_name
    )
  } else if (is(features, 'tbl_sql')) {
    features <- mutate(collect(features, n = Inf))
    if (is.null(features[["assay"]])) {
      features[["assay"]] <- assay_name
    }
  }
  
  add.cols <- setdiff(c("assay", "assay_type", "feature_type"), colnames(features))  
  for (add in add.cols) {
    features[[add]] <- ainfo[[add]]
  }

  features <- dplyr::distinct(
    features, 
    .data$feature_id, 
    .data$assay,
    .data$assay_type,
    .data$feature_type)
  
  assert_assay_feature_descriptor(features, x)
  features
}