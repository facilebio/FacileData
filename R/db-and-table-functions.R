#' Query a table to identify its primary key(s)
#'
#' @export
#' @importFrom DBI dbGetQuery
#' @param x a \code{FacileDataSet} or \code{SQLiteConnection}
#' @param table_name the name of the table to query
#' @return a character vector of primary keys
primary_key <- function(x, table_name) {
  if (is.FacileDataSet(x)) x <- x$con
  stopifnot(is(x, 'SQLiteConnection'))
  assert_string(table_name)
  info <- dbGetQuery(x, sprintf("PRAGMA table_info(%s);", table_name))
  filter(info, pk != 0)$name
}

#' Adds rows to a table in a FacileDataSet
#'
#' This function first checks the data in the target table \code{table_name}
#' to ensure that rows in \code{dat} that exist in \code{table_name} (by
#' checking the primary key) are not added.
#'
#' @export
#' @importFrom DBI dbWriteTable
#' @param dat the \code{data.frame} of rows to add to the table, which must
#'   have a superset of columns present in the \code{table_name} that is being
#'   appended to
#' @param x the \code{FacileDataSet}
#' @param table_name the name of the table in \code{x} to add the rows of
#'   \code{dat} to.
#' @return invisibly returns the conformed version of \code{dat}.
append_facile_table <- function(dat, x, table_name, warn_existing = FALSE) {
  stopifnot(is.FacileDataSet(x))
  target <- try(tbl(x$con, table_name), silent=TRUE)
  if (is(target, 'try-error')) stop("Unknown table to append to: ", table_name)
  dat <- conform_data_frame(dat, target)
  # strip facile_frame class if it's there.
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)

  # Ensure that we don't try to add existing rows into the database
  pk <- primary_key(x, table_name)
  if (length(pk)) {
    skip <- target |>
      semi_join(dat, by=pk, copy=TRUE, auto_index=TRUE) |>
      collect(n=Inf) |>
      mutate(added=FALSE)
    if (nrow(skip) && warn_existing) {
      warning(nrow(skip), "/", nrow(dat), " features already in database",
              immediate.=TRUE)
    }
    add.me <- anti_join(dat, skip, by=pk)
    if (nrow(add.me)) {
      dbWriteTable(x$con, table_name, add.me, append=TRUE)
      add.me$added <- TRUE
    }
    dat <- bind_rows(add.me, skip)
  } else {
    dat$added <- TRUE
    dbWriteTable(x$con, table_name, dat, append=TRUE)
  }

  invisible(dat)
}

# Database Table Accessors =====================================================

#' @export
#' @noRd
assay_info_tbl <- function(x) {
  UseMethod("assay_info_tbl", x)
}

#' @export
assay_info_tbl.FacileDataSet <- function(x) {
  out <- tbl(x$con, 'assay_info')
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}

#' @export
#' @noRd
assay_feature_info_tbl <- function(x) {
  UseMethod("assay_feature_info_tbl", x)
}

#' @export
#' @noRd
assay_feature_info_tbl.FacileDataSet <- function(x) {
  out <- tbl(x$con, 'assay_feature_info')
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}

#' @export
#' @noRd
assay_sample_info_tbl <- function(x) {
  UseMethod("assay_sample_info_tbl", x)
}

#' @export
#' @noRd
assay_sample_info_tbl.default <- function(x) {
  # Currently we are accessing directly the assay_sample_info tbl to get
  # assay_sample covariates, which needs to change. See issue #2:
  # https://github.com/facilebio/FacileData/issues/2
  stop("assay_sample_info_tbl not implemented for: ", class(x))
}

#' @export
#' @noRd
assay_sample_info_tbl.FacileDataSet <- function(x) {
  out <- tbl(x$con, 'assay_sample_info')
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}

#' @export
#' @noRd
feature_info_tbl <- function(x, assay_name = NULL) {
  UseMethod("feature_info_tbl", x)
}

#' @export
#' @noRd
feature_info_tbl.FacileDataSet <- function(x, assay_name = NULL) {
  out <- tbl(x$con, 'feature_info')
  if (!is.null(assay_name)) {
    assert_string(assay_name)
    assay.info <- assay_info_tbl(x) |>
      filter(assay == assay_name) |>
      collect()
    if (nrow(assay.info) == 0) {
      stop("Unknown assay: ", assay_name)
    }
    afi <- assay_feature_info_tbl(x) |>
      filter(assay == assay_name)
    out <- semi_join(out, afi, by=c('feature_type', 'feature_id'))
  }
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}

#' @export
#' @noRd
gene_info_tbl <- function(x) {
  UseMethod("gene_info_tbl", x)
}

#' Mimics the old `gene_info` table.
#'
#' @export
gene_info_tbl.FacileDataSet <- function(x) {
  # TODO: This function needs to be removed and the code that relies on
  # gene_info_tbl should be updated.
  ## Columns:
  ## feature_id|feature_type|symbol|n_exons|length|source|hdf5_index
  hdf5.info <- assay_feature_info_tbl(x) |>
    filter(assay == default_assay(x))

  gi <- feature_info_tbl(x) |>
    filter(feature_type == 'entrez') |>
    select(feature_id, feature_type, symbol=name, n_exons=-1,
           # length=effective_length,
           source) |>
    inner_join(hdf5.info, by='feature_id') |>
    as_facile_frame(x, .valid_sample_check = FALSE)
}

#' Mimics old sample_stats table
#'
#' This function needs to be removed and the code that relies on
#' sample_stats_tbl be updated.
#' @export
#' @noRd
sample_stats_tbl <- function(x) {
  UseMethod("sample_stats_tbl")
}

#' @export
#' @noRd
sample_stats_tbl.FacileDataSet <- function(x) {
  assay_sample_info_tbl(x) |>
    select(dataset, sample_id, libsize, normfactor) |>
    as_facile_frame(x, .valid_sample_check = FALSE)
}

#' @export
#' @noRd
sample_covariate_tbl <- function(x) {
  UseMethod("sample_covariate_tbl", x)
}

#' @export
#' @noRd
sample_covariate_tbl.FacileDataSet <- function(x) {
  out <- tbl(x$con, 'sample_covariate')
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}

#' @export
#' @noRd
sample_info_tbl <- function(x) {
  UseMethod("sample_info_tbl", x)
}

#' @export
#' @noRd
sample_info_tbl.FacileDataSet <- function(x) {
  out <- tbl(x$con, 'sample_info')
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}

## Unexported utility functions ================================================

#' Validates the bits required in a legit FacileDataSet directory.
#' @noRd
validate.facile.dirs <- function(path, data.fn, sqlite.fn, hdf5.fn, meta.fn,
                                 anno.dir) {
  if (!dir.exists(path)) {
    stop("Top level FacileData directory does not exist: ", path)
  }
  path <- normalizePath(path)
  if (!file.exists(data.fn)) {
    stop("Data file (database) does not exists", data.fn)
  } else {
    data.fn <- normalizePath(data.fn)
    if (dirname(data.fn) != path) {
      warning("Data file is not under parent directory", immediate.=TRUE)
    }
  }
  if (!file.exists(sqlite.fn)) {
    stop("Database file does not exists", sqlite.fn)
  } else {
    sqlite.fn <- normalizePath(sqlite.fn)
    if (dirname(sqlite.fn) != path) {
      warning("Database file is not under parent directory", immediate.=TRUE)
    }
  }
  if (!file.exists(hdf5.fn)) {
    warning("HDF5 file does not exists", hdf5.fn, immediate.=TRUE)
  } else {
    hdf5.fn <- normalizePath(hdf5.fn)
    if (dirname(hdf5.fn) != path) {
      warning("HDF5 file is not under parent directory", immediate.=TRUE)
    }
  }
  meta.fn <- assert_valid_meta_file(meta.fn) |> normalizePath()
  if (!dir.exists(anno.dir)) {
    stop("Directory for custom annotations does not exist: ", anno.dir)
  } else {
    anno.dir <- normalizePath(anno.dir)
    if (dirname(anno.dir) != path) {
      warning("Custom annotation directory not under parent directory.",
              immediate.=TRUE)
    }
  }

  list(path=path, data.fn=data.fn, sqlite.fn=sqlite.fn, hdf5.fn=hdf5.fn,
       meta.fn=meta.fn, anno.dir=anno.dir)
}
