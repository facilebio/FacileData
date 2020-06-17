#' @noRd
#' @export
biocbox.FacileDataStore <- function(x, class = NULL, assay_name = NULL,
                                    features = NULL, sample_covariates = NULL,
                                    feature_covariates = NULL,
                                    normalized = FALSE, with_fds = FALSE,
                                    custom_key = Sys.getenv("USER"), ...) {
  xs <- samples(x)
  biocbox(xs, class = class, assay_name = assay_name, features = features,
          normalized = normalized, with_fds = with_fds, custom_key = custom_key,
          ...)
}

#' Assembles a Bioconductor data container for a given assay.
#'
#' There is a default bioc class provided for different assay types, however
#' the class type can be overrided by the `class` parameter. This function
#' simply puts the assay data requested into the container. There is no
#' sepcial functionality that happens downstream of that (for instance,
#' DGEList lib.size calculated from the data that made its way into the DGEList)
#'
#' @noRd
#' @export
#'
#' @param sample_covariates If `NULL` (default), all sample covariates will
#'   be included over samples in x. If a data.frame, we will treat the
#'   extra columns as custom covariates, and include them in the outgoing
#'   box, along with the internal ones.
biocbox.facile_frame <- function(x, class = NULL, assay_name = NULL,
                                 features = NULL, sample_covariates = NULL,
                                 feature_covariates = NULL,
                                 normalized = FALSE, with_fds = FALSE,
                                 custom_key = Sys.getenv("USER"), ...) {
  if (!is.null(feature_covariates)) {
    warning("We currently are not trying to merge extra feature_covariates.",
            immediate. = TRUE)
  }
  assert_sample_subset(x)
  extra.covs <- setdiff(colnames(x), c("dataset", "sample_id"))
  fds. <- assert_facile_data_store(fds(x))

  if (is.null(assay_name)) {
    assay_name <- default_assay(fds.)
  } else {
    assert_choice(assay_name, assay_names(fds.))
  }

  ainfo <- assay_info(fds., assay_name)
  bb.info <- biocbox_assay_class(ainfo[["assay_type"]], class)
  assert_flag(normalized)
  if (normalized && bb.info[["class"]] != "list") {
    stop("Can't return normalized version of data in bioc container")
  }

  # setup pData
  if (isFALSE(sample_covariates)) {
    samples. <- x
  } else if (is.null(sample_covariates)) {
    samples. <- with_sample_covariates(x, custom_key = custom_key)
  } else {
    if (test_sample_subset(sample_covariates)) {
      samples. <- left_join(x, sample_covariates,
                            by = c("dataset", "sample_id"),
                            suffix = c("", ".extra"))
      samples. <- with_sample_covariates(samples., custom_key = custom_key)
    } else {
      assert_character(sample_covariates)
      samples. <- with_sample_covariates(x, sample_covariates,
                                         custom_key = custom_key)
    }
  }

  samples. <- mutate(samples., samid = paste(dataset, sample_id, sep = "__"))

  # setup fData
  fids <- NULL
  if (!is.null(features)) {
    if (is.data.frame(features)) {
      features <- features[["feature_id"]]
    }
    if (is.factor(features)) features <- as.character(features)
    fids <- assert_character(features)
  }

  features. <- features(fds., assay_name = assay_name, feature_ids = fids)
  if (is.data.frame(features)) {
    # User provided more metadata in the features data.frame passed in
    take <- setdiff(colnames(features), colnames(features.))
    if (length(take)) {
      take <- c("feature_id", take)
      features. <- left_join(features., select(features, {{take}}),
                             by = "feature_id")
    }
  }
  axe.f.cols <- c("assay", "hd5_index", "assay_type")
  features. <- select(features., -any_of(axe.f.cols))
  features. <- as.data.frame(features.)
  if (!"symbol" %in% colnames(features.)) {
    features.[["symbol"]] <- features.[["name"]]
  }
  rownames(features.) <- features.[["feature_id"]]

  A <- fetch_assay_data(fds., fids, samples., assay_name = assay_name,
                        normalized = normalized, as.matrix = TRUE, ...)

  fxref <- match(rownames(A), features.[["feature_id"]])
  if (any(is.na(fxref))) {
    stop("There are rows in the assay matrix that we don't have feature info ",
         "for, this should not have happend")
  }
  features. <- features.[fxref,,drop = FALSE]

  sxref <- match(colnames(A), samples.[["samid"]])
  if (any(is.na(sxref))) {
    stop("There are columns in the assay matrix that we don't have sample ",
         "info for, this should not have happend")
  }
  samples. <- samples.[sxref,,drop = FALSE]
  samples. <- as.data.frame(samples.)
  rownames(samples.) <- samples.[["samid"]]
  samples.[["samid"]] <- NULL

  out <- bb.info[["fn"]](A, features., samples., fds.,
                         update_libsizes = update_libsizes,
                         update_normfactors = update_normfactors, ...)
  if (!with_fds) {
    out <- set_fds(out, NULL)
  }

  out
}

# Bioconductor Creation Helpers ------------------------------------------------

#' Map assay_types to default bioc containers
#' @noRd
.biocboxes <- tribble(
  ~assay_type,     ~class,                   ~package,
  "rnaseq",        "DGEList",                "edgeR",
  "rnaseq",        "SummarizedExperiment",   "SummarizedExperiment",
  "rnaseq",        "ExpressionSet",          "Biobase",
  "pseudobulk",    "DGEList",                "edgeR",
  "lognorm",       "EList",                  "limma")

#' Returns the class info to use to instantiate the right biocbox class
#'
#' @noRd
#' @param assay_type the type of assay to stuff in a box
#' @param class the specific class requested. If `NULL` (default), a default
#'   class will be used.
biocbox_assay_class <- function(assay_type, class = NULL) {
  assert_string(assay_type)
  if (!is.null(class) && class == "list") {
    return(list(assay_type = assay_type, class = "list", fn = bb.list))
  }

  if (!is.null(class)) assert_choice(class, .biocboxes$class)

  if (!assay_type %in% .biocboxes$assay_type) {
    if (is.null(class)) {
      warning("assay_type (", assay_type, ") not found in class list, and no ",
              "class was provided, returning an ExpressionsSet")
      class <- "ExpressionSet"
    } else {
      warning("assay_type (", assay_type, ") not found in class list, ",
              "giving you what you asked for: ", class)
    }
  } else {
    choices <- filter(.biocboxes, .data$assay_type == .env$assay_type)
    if (is.null(class)) {
      class <- choices[["class"]][1L]
    } else {
      if (!class %in% choices[["class"]]) {
        warning("asking for a `", class, "` container, which is not compatible",
                " with the assay_type: ", assay_type)
      }
    }
  }

  fn.name <- paste0("bb.", class)
  fn <- getFunction(fn.name)
  list(assay_type = assay_type, class = class, fn = fn)
}

#' @noRd
bb.list <- function(amatrix, fdat, pdat, fds., ...) {
  set_fds(list(assay_data = amatrix, features = fdat, samples = pdat), fds.)
}

#' @noRd
bb.ExpressionSet <- function(amatrix, fdat, pdat, fds., ...) {
  reqpkg("Biobase")
  out <- Biobase::ExpressionSet(amatrix)
  out <- Biobase::`pData<-`(out, pdat)
  out <- Biobase::`fData<-`(out, fdat)
  set_fds(out, fds.)
}

#' @noRd
bb.SummarizedExperiment <- function(amatrix, fdat, pdat, fds., ...) {
  reqpkg("SummarizedExperiment")
  out <- SummarizedExperiment::SummarizedExperiment(
    list(counts = amatrix),
    rowData = fdat,
    colData = pdat)
  set_fds(out, fds.)
}

#' @noRd
#' @importFrom edgeR DGEList
bb.DGEList <- function(amatrix, fdat, pdat, fds., ...) {
  # Issue #2
  # https://github.com/denalitherapeutics/FacileData/issues/2
  # sample.stats <- samples %>%
  #   with_assay_covariates(c("libsize", "normfactor"), assay_name,
  #                         .fds = .fds) %>%
  #   collect(n = Inf) %>%
  #   mutate(samid=paste(dataset, sample_id, sep='__')) %>%
  #   as.data.frame()
  # rownames(sample.stats) <- sample.stats[["samid"]]
  # sample.stats <- sample.stats[colnames(x),,drop=FALSE]
  out <- DGEList(amatrix, genes = fdat, samples = pdat)
  set_fds(out, fds.)
}

#' @noRd
bb.EList <- function(amatrix, fdat, pdat, fds., ...) {
  out <- new("EList", list(E = amatrix, genes = fdat, targets = pdat))
  set_fds(out, fds.)
}
