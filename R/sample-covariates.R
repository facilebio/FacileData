#' Provides a summary table of sample covariates.
#'
#' Sumamrizes a set of sample covariates (returned from
#' [fetch_sample_covariates()] at different granulaities.
#'
#' @md
#' @export
#'
#' @param object A sample covariate table, the likes returned from
#'   [fetch_sample_covariates()].
#' @param expanded includes details (rows) for each covariate per level
#'   (or quantile), depending on the covariates `"class"` attribute.
#' @return a tibble of summary sample covariate information with the following
#'   columns:
#'   * `variable`: name of the variable
#'   * `class`: class of variable (real, categorical)
#'   * `nsamples`: the number of samples that have this variable defined
#'   * `level`: the level (or quantile) of the covariate
#'     (included only when `expanded == TRUE`)
#'   * `ninlevel`: the number of samples with this covariate value
#'     (included only when `expanded == TRUE`)
#' @examples
#' fds <- exampleFacileDataSet()
#' covs <- fetch_sample_covariates(fds)
#' smry <- summary(covs)
#' details <- summary(covs, expanded = TRUE)
#' catdeetz <- covs |>
#'   filter(class == "categorical") |>
#'   summary(expanded = TRUE)
summary.eav_covariates <- function(object, ..., expanded = FALSE,
                                   droplevels = TRUE) {
  object <- assert_sample_covariates(object)
  .fds <- assert_facile_data_store(fds(object))
  with.source <- is.character(object[["source"]])

  dat <- collect(object, n = Inf)
  if (!with.source) dat <- mutate(dat, source = NA)

  if (expanded) {
    covdef <- covariate_definitions(.fds)
    res <- dat |>
      group_by(variable, class) |>
      do({
        value <- cast_covariate(.$variable[1L], .$value, covdef, .fds)
        clz <- .$class[1L]

        if (clz %in% c("categorical", "logical") && is.atomic(value)) {
          levels <- table(value)
        } else if (clz == "real" && is.atomic(value)) {
          # small vectors, or those with a lot of NAs produce quantiles that
          # are not unique, and throw errors
          qtl <- unique(quantile(value))
          if (length(qtl) == 1L) {
            qtl <- sort(c(qtl  -1, qtl))
          }
          bins <- cut(value, qtl, include.lowest = TRUE)
          levels <- table(bins)
        } else {
          levels <- c(all = length(value))
        }
        out <- tibble(
          nsamples = length(value),
          source = .$source[1L],
          level = names(levels),
          ninlevel = as.integer(levels))
        if (droplevels) {
          out <- filter(out, ninlevel > 0L)
        }
      })
  } else {
    res <- dat |>
      group_by(variable, class) |>
      summarize(
        ndatasets = length(unique(dataset)),
        nsamples = n(),
        nlevels = {
          if (class[1L] %in% c("categorical", "logical"))
            length(unique(value))
          else
            NA_integer_
        },
        IQR = {
          if (class[1L] == "real")
            IQR(as.numeric(value))
          else
            NA_real_
        })
  }

  res <- ungroup(res)
  if (!with.source) {
    res <- mutate(res, source = NULL)
  }

  as_facile_frame(res, .fds, .valid_sample_check = FALSE)
}

#' Summaries of long covariates look like this:
#' expanded = FALSE
#' variable         class       ndatasets nsamples nlevels   IQR
#' <chr>            <chr>           <int>    <int>   <int> <dbl>
#' 1 cell_abbrev      categorical         3      130      11   NA 
#' 2 cell_type        categorical         3      130      11   NA 
#' 3 cond             categorical         3      130       3   NA 
#' 4 condition        categorical         3      130       3   NA 
#' 5 diabetes_history categorical         3      130       2   NA 
#' 6 donor_id         categorical         3      130      14   NA 
#' 7 eGFR             categorical         3      130       4   NA 
#' 8 group            categorical         3      130      33   NA 
#' 9 hypertension     categorical         3      130       2   NA 
#' 10 ncells           real                3      130      NA  360.
#' 11 sex              categorical         3      130       2   NA 
#' 
#' @export
#' @noRd
summary.wide_covariates <- function(object, ..., expanded = FALSE,
                                    droplevels = TRUE) {
  if (FALSE) {
    object <- an_fds() |> 
      samples() |> 
      with_sample_covariates()
  }
  cols <- setdiff(colnames(object), c("dataset", "sample_id"))
  nvars <- length(cols)

  out <- lapply(seq_along(cols), function(idx) {
    vname <- cols[idx]
    vals.all <- object[[vname]]
    notna <- !is.na(vals.all)
    vals <- vals.all[notna]
    eavdef <- eavdef_for_column(vals.all, vname)
    vclass <- eavdef$class
    ndatasets <- length(unique(object[["dataset"]][notna]))
    nsamples <- sum(notna)
    
    if (vclass %in% c("categorical", "logical")) {
      vtally <- table(vals)
      iqr <- NA_real_
      nlevels <- length(vtally)
    } else if (vclass == "real") {
      iqr <- IQR(vals)
      # small vectors, or those with a lot of NAs produce quantiles that
      # are not unique, and throw errors
      qtl <- unique(quantile(vals))
      if (length(qtl) == 1L) {
        qtl <- sort(c(qtl - 1, qtl))
      }
      bins <- cut(vals, qtl, include.lowest = TRUE)
      vtally <- table(bins)
      nlevels <- NA_integer_
    }
  
    if (expanded) {
      res <- dplyr::tibble(
        variable = vname,
        class = vclass,
        nsamples = nsamples,
        level = names(vtally),
        ninlevel = as.numeric(vtally))
    } else {
      res <- dplyr::tibble(
        variable = vname,
        class = vclass,
        ndatasets = ndatasets,
        nsamples = nsamples,
        nlevels = nlevels,
        IQR = iqr)
    }
    res
  })
  out <- dplyr::bind_rows(out)
  as_facile_frame(out, fds(object), .valid_sample_check = FALSE)
}

sample_covariates.facile_frame <- function(x, ...) {
  stop("what is this for")
  .fds <- fds(x)
}

#' Fetch rows from sample_covariate table for specified samples and covariates
#'
#' @export
#' @rdname sample-covariates
#' @family API
#'
#' @param x a \code{FacileDataSet} connection
#' @param samples a samples descriptor \code{tbl_*}
#' @param covariates character vector of covariate names
#' @param custom_key The key to use to fetch more custom annotations over
#'   the given samples
#' @return rows from the \code{sample_covariate} table
#' @examples
#' # Retrieves a long EAV table of covariates
#' long_covs <- an_fds() |> 
#'   fetch_sample_covariates(c("cell_abbrev", "sex", "condition"))
#' lc2 <- some_samples() |> 
#'   fetch_sample_covariates(c("cell_abbrev", "sex", "condition"))
#' wide_covs <- some_samples() |> 
#'   with_sample_covariates(c("cell_abbrev", "sex", "condition"))
fetch_sample_covariates.FacileDataSet <- function(
    x, covariates = NULL, samples = NULL, custom_key = Sys.getenv("USER"),
    with_source = FALSE, ...) {
  if (is.null(samples)) {
    samples <- samples(x)
  } else {
    assert_sample_subset(samples, x)
    samples <- distinct(samples, dataset, sample_id)
  }
  samples <- collect(samples, n = Inf)

  dat <- sample_covariate_tbl(x)
  if (is.character(covariates)) {
    if (length(covariates) == 1L) {
      dat <- filter(dat, variable == covariates)
    } else if (length(covariates) > 1L) {
      dat <- filter(dat, variable %in% covariates)
    }
  }
  dat <- collect(dat, n=Inf)
  dat <- set_fds(dat, x) ## explicitly added here to do `collect` above

  # out <- join_samples(dat, samples, semi=TRUE)
  # out <- collect(out, n = Inf)
  out <- semi_join(dat, samples, by = c("dataset", "sample_id"))
  if (with_source) {
    out <- mutate(out, source = "datastore")
  }

  if (!is.null(custom_key)) {
    custom <- fetch_custom_sample_covariates(x, covariates = covariates,
                                             samples = samples,
                                             custom_key = custom_key,
                                             with_source = with_source)
    out <- bind_rows(collect(out, n=Inf), custom)
  }

  as_facile_frame(out, x, "eav_covariates", .valid_sample_check = FALSE)
}

#' @export
#' @rdname sample-covariates
#' @family API
fetch_sample_covariates.facile_frame <- function(
    x, covariates = NULL, samples = NULL,
    custom_key = Sys.getenv("USER"), with_source = FALSE, ...) {
  if (!is.null(samples)) {
    warning("`samples` ignored when fetching covariates from a facile_frame",
            immediate. = TRUE)
  }

  samples. <- assert_sample_subset(x)
  samples. <- distinct(samples., dataset, sample_id)

  fetch_sample_covariates(fds(x), covariates = covariates, samples = samples.,
                          custom_key = custom_key, with_source = with_source,
                          ...)
}

#' Fetches custom (user) annotations for a given user prefix
#'
#' @export
#' @importFrom jsonlite stream_in
#' @param fds The \code{FacileDataSet}
#' @param samples the facile sample descriptor
#' @param custom_key The key to use for the custom annotation
#' @return covariate tbl
#' @family API
fetch_custom_sample_covariates.FacileDataSet <- function(
    x, covariates = NULL, samples = NULL, custom_key = Sys.getenv("USER"),
    with_source = FALSE, file.prefix = "facile", ...) {
  if (is.null(samples)) {
    samples <- samples(x)
  } else {
    assert_sample_subset(samples, x)
    samples <- distinct(samples, dataset, sample_id)
  }
  samples <- collect(samples, n = Inf)

  out.cols <- c("dataset", "sample_id", "variable", "value", "class", "type",
                "date_entered")
  fpat <- paste0('^', file.prefix, '_', custom_key, "_.*json")
  annot.files <- list.files(path=x$anno.dir, pattern=fpat, full.names=TRUE)

  if (length(annot.files)) {
    annos <- lapply(annot.files, function(fn) stream_in(file(fn), verbose=FALSE))
    out <- bind_rows(annos) |>
      # select_(.dots=out.cols) |>
      select(!!out.cols) |>
      set_fds(x) |>
      join_samples(samples, semi=TRUE)
    ## We weren't saving the type == 'categorical' column earlier. So if this
    ## column is.na, then we force it to 'categorical', because that's all it
    ## realy could have been
    if (nrow(out) && mean(is.na(out$type)) > 0.5) {
      out <- mutate(out, class='categorical', type='user_annotation')
    }
  } else {
    ## Make a dummy, 0 row tibble to send back
    out <- sapply(out.cols, function(x) character(), simplify=FALSE)
    out$date_entered <- integer()
    out <- as.data.frame(out, stringsAsFactors=FALSE) |> as_tibble()
  }

  if (!is.null(covariates)) {
    out <- filter(out, variable %in% covariates)
  }

  if (with_source) {
    out <- mutate(out, source = "userstore")
  }

  as_facile_frame(out, x, "eav_covariates", .valid_sample_check = FALSE)
}

#' Saves custom sample covariates to a FacileDataSet
#'
#' @export
#' @importFrom jsonlite stream_out
#'
#' @param x the \code{FacileDataSet}
#' @param annotation the annotation table of covariate values to a
#'   sample-descriptor-like table
#' @param name the variable name of the covariate
#' @param custom_key the custom key (likely userid) for the annotation
#' @param file.prefix Vincent uses this
#' @param sample_filter_criteria optional list of filtering criteria that were
#'   used to drill down into the samples we have the \code{annotatino}
#'   data.frame for
save_custom_sample_covariates <- function(x, annotation, name = NULL,
                                          class = "categorical",
                                          custom_key = Sys.getenv("USER"),
                                          file.prefix = "facile",
                                          sample_filter_critera = NULL) {
  #' TODO: Figure out how to encode sample_filter_criteria into serialized
  #' (JSON) annotation file
  stopifnot(is.FacileDataSet(x))
  annotation <- collect(annotation, n=Inf)
  assert_columns(annotation, c('dataset', 'sample_id', 'value'))
  if (is.null(name)) name <- annotation$name
  if (!test_string(name)) stop("No name given/inferred for custom annotation")

  if (is.null(custom_key)) custom_key <- 'anonymous'
  custom_key <- assert_string(custom_key) |> make.names()

  annotation[['variable']] <- make.names(name)
  annotation <- annotation[, c('dataset', 'sample_id', 'variable', 'value')]
  annotation[['class']] <- class
  annotation[['type']] <- 'user_annotation'
  annotation[['date_entered']] <- as.integer(Sys.time())

  fn <- paste0(file.prefix, '_', custom_key, '_', name, '_', Sys.Date(),'.json')
  fn <- file.path(x$anno.dir, fn)
  stream_out(x=annotation, con=file(fn))
  invisible(set_fds(annotation, x))
}

#' @export
#' @noRd
with_sample_covariates.facile_frame <- function(x, covariates = NULL,
                                                na.rm = FALSE,
                                                custom_key = Sys.getenv("USER"),
                                                .fds = fds(x), ...) {
  x <- collect(x, n = Inf)
  .fds <- assert_facile_data_store(.fds)
  NextMethod(x, .fds = .fds)
}

#' @export
#' @noRd
with_sample_covariates.tbl <- function(x, covariates = NULL,
                                       na.rm = FALSE,
                                       custom_key = Sys.getenv("USER"),
                                       .fds = NULL, ...) {
  with_sample_covariates.data.frame(collect(x, n = Inf),
                                    covariates = covariates,
                                    na.rm = na.rm, custom_key = custom_key,
                                    .fds = .fds, ...)
}

#' @export
#' @noRd
#' @method with_sample_covariates data.frame
#' @examples
#' efds <- exampleFacileDataSet()
#' s <- filter_samples(efds, indication == "CRC")
#' with_sample_covariates(s, c("sample_type", "stage"))
#' with_sample_covariates(s, c("sample_type", tumor_stage = "stage"))
with_sample_covariates.data.frame <- function(x, covariates = NULL,
                                              na.rm = FALSE,
                                              custom_key = Sys.getenv("USER"),
                                              .fds = NULL, ...) {
  assert_facile_data_store(.fds)
  x <- assert_sample_subset(x) |> collect(n=Inf)
  stopifnot(is.character(covariates) || is.null(covariates))
  if (is.character(covariates)) {
    if (length(covariates) == 0L) return(x)
    covariates <- covariates[!duplicated(covariates)]
    covariates <- nameit(covariates)
  }

  samples <- x |>
    select(dataset, sample_id) |>
    distinct(.keep_all=TRUE)

  covs <- fetch_sample_covariates(.fds, covariates = covariates,
                                  samples = samples,
                                  custom_key = custom_key, ...)
  covs <- spread_covariates(covs, .fds, ...)
  if (!is.null(covariates) && !is.null(names(covariates))) {
    covs <- rename(covs, !!covariates)
  }

  out <- left_join(x, covs, by=c("dataset", "sample_id"),
                   suffix = c("", ".infds"))

  if (na.rm && length(covariates)) {
    keep <- complete.cases(select_(out, .dots=covariates))
    out <- out[keep,,drop=FALSE]
  }

  as_facile_frame(out, .fds, "wide_covariates", .valid_sample_check = FALSE)
}

#' Spreads the covariates returned from database into wide data.frame
#'
#' Samples that did not have a value for a specific covariate are assigned to
#' have NA.
#'
#' @export
#' @importFrom data.table dcast setDT setDF
#' @param x output from \code{fetch_sample_covariates}
#' @param .fds A \code{FacileDataSet} object
#' @return a wide \code{tbl_df}-like object
spread_covariates <- function(x, .fds = fds(x), cov.def = NULL, ...) {
  assert_facile_data_store(.fds)
  x <- assert_sample_covariates(x) |>
    collect(n=Inf)

  ## Ensures we get a row for every sample in x, even if it is missing a value
  ## for the covariate
  dummy <- select(x, dataset, sample_id) |>
    collect(n=Inf) |>
    distinct(.keep_all=TRUE) |>
    mutate(variable='.dummy.', value=NA)

  out <- bind_rows(x, dummy) |>
    setDT() |>
    dcast(dataset + sample_id ~ variable, value.var = "value") |>
    setDF() |>
    mutate(.dummy. = NULL) |>
    as_tibble()

  # sample-covariates -------------------------------------------------------
  if (is.null(cov.def)) {
    cov.def <- covariate_definitions(.fds)
  }

  if (!is.null(cov.def)) {
    do.cast <- setdiff(colnames(out), c('dataset', 'sample_id'))
    ## Don't decode categorical variables of type 'user_annotation'
    user.anno <- filter(x, type == 'user_annotation' & class == 'categorical')
    if (nrow(user.anno)) {
      do.cast <- setdiff(do.cast, unique(user.anno$variable))
    }

    for (cname in do.cast) {
      casted <- cast_covariate(cname, out[[cname]], cov.def)
      if (is.data.frame(casted)) {
        ## casting a survival covariate will return a two column thing with time
        ## and censoring information, so we need to account for that.
        out[[cname]] <- NULL
        out <- bind_cols(out, casted)
      } else {
        out[[cname]] <- casted
      }
    }
  }

  as_facile_frame(out, .fds, "wide_covariates", .valid_sample_check = FALSE)
}
