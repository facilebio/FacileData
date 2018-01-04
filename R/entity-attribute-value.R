# What we know as `pData` in the R/bioc/ExpressionSet world are stored in an
# entity-attribute-value `sample_covariate` table in the FacileDataSet.
# The functions defined here support `pData` <-> eav conversion.

# Extracting data from EAV table into R ========================================

#' Casts the character values of the covariates to their defined types.
#'
#' For most things, a single value will be returned from each cast, but in the
#' case of "time_to_event" data, the value is expended to a two column
#' data.frame with a `tte_<covariate>` column for time to event, and an
#' `event_<covariate>` column to indicate event (1) or right censored (2).
#'
#' @md
#' @export
#'
#' @param covariate the name of the covariate
#' @param values the covariate values (which is a `character`) as it is
#'   pulled from the database.
#' @param cov.def the un-yamled covariate definitions, if missing we rely on
#'   pulling this out from the `FacileDataSet` object `.fds`
#' @param .fds If `missing(cov.def)`, this is the `FacileDataSet` to
#'   get the covariate definitions from.
#' @return values cast to appropriate type if a valid definition was found for
#'   `covariate`, otherwise values is returned "as is". Most of the time
#'   this is a single vector, but others it can be a data.frame (for
#'   `right_censored` data, for instance)
cast_covariate <- function(covariate, values, cov.def, .fds) {
  if (missing(cov.def)) {
    stopifnot(is.FacileDataSet(.fds))
    cov.def <- covariate_definitions(.fds)
  }
  stopifnot(is(cov.def, 'CovariateDefinitions'))
  stopifnot(is.character(values))
  stopifnot(is.character(covariate) && length(covariate) == 1L)

  def <- cov.def[[covariate]]
  if (is.list(def)) {
    if (def$class == 'real') {
      values <- as.numeric(values)
    }
    if (def$class == 'right_censored') {
      values <- decode_right_censored(values, suffix=covariate)
    }
    if (is.character(def$levels)) {
      ## protect against NAing a long list of values in the event that the
      ## levels provided in the meta.yaml file don't include all levels observed
      ## here
      lvls <- c(def$levels, setdiff(values, def$levels))
      values <- factor(values, lvls)
    }
  } else {
    if (covariate != '.dummy.') {
      warning("No covariate definition found for: ", covariate, immediate.=TRUE)
    }
  }

  values
}

#' Retrieve the meta information about a covariate
#'
#' @md
#' @export
#'
#' @param covariate the name of the covariate
#' @param .fds the `FacileDataSet`
#' @param covdefs The `covariate_definitions(.fds)` list
#' @return a list of covariate information with the following elements:
#'   `$name`, `$type`, `$class`, `$description`,
#'   `$label`, `$is.factor`, (and maybe `$levels`)
covariate_meta_info <- function(covariate, .fds, covdefs=NULL) {
  if (is.null(covdefs)) {
    stopifnot(is.FacileDataSet(.fds))
    covdefs <- covariate_definitions(.fds)
  }
  assert_covariate_definitions(covdefs)
  meta <- covdefs[[covariate]]
  if (!is.list(meta)) {
    stop("Covariate `", covariate, "` not found in covariate_definition file")
  }
  meta$name <- covariate
  meta$is.factor <- meta$class == 'categorical' && is.character(meta$levels)
  meta
}

# Creating attribute-value definitions =========================================

#' Create a facile covariate definition file from a sample `pData` data.frame
#'
#' Sample covariates (aka `pData`) are encoded in an
#' [entity-attribute-value (EAV) table](https://en.wikipedia.org/wiki/Entity-attribute-value_model).
#' Metadata about these covariates are stored in a `meta.yaml` file in the
#' `FacileDataSet` directory which enables the `FacileDataSet` to cast the value
#' stored in the EAV table to its native R type. This function generates the
#' list-of-list structure to represent the `sample_covariates` section of the
#' `meta.yaml` file.
#'
#' For simple `pData` covariates, each column is treated independantly from the
#' rest. There are some types of covariates which require multiple columns for
#' proper encoding, such as encoding survival information (see the **Survival**
#' section below). In these cases, the caller needs to provide an entry in the
#' `covariate_def` list that describes which `pData` columns (`varname`) goes
#' into the single facile covariate value. Please reference the
#' **Defining EAV Encodings** section for more information.
#'
#' keeping track of the "time to event" in one column, and a separae column to
#' indicate whether or not the event was a "death" or a censored. Still, these
#' data are stored in the single "value" column of the FacileDataSet's internal
#' entity-attribute-value (`sample_covariate`) table. In order to
#' encode these types of columns correctly, we need to provide more information
#' via the `covariate_def` parameter of this function.
#'
#' @section Survival:
#' Survival data is encoded by two columns. One column to indicate the
#' "time to event" and a second to indicate whether or not the denoted
#' tte is an "event" (1) or "censored" (0). The pair of columns will be encoded
#' into the `FacileDataSet`'s EAV table as a single (numeric) value. The
#' absolute value of the numeric indicates the "time to event" and the sign of
#' the value indicates its censoring status (If there are such data in `x`, it
#' must be in a (`tte_OS`, `event_OS`) pair of columns for "ordinary survival"
#' or a (`tte_PFS`, `event_PFS`) for progression free survival. Please see the
#' "Defining EAV Encodings" sections for more details.
#'
#' @md
#' @rdname eav-metadata
#' @export
#'
#' @param x a `pData` `data.frame`
#' @param covariate_def a named list of covariate definitions. The names of
#'   this list are the names the covariates will be called in the target
#'   `FacileDataSet`. The values of the list are:
#'   * `varname`: a `character()` of the column name(s) in `x` that this
#'      sample covariate was derived from. If more than one column is to be used
#'      for the facile covariate conversion (eg. if we are encoding survival),
#'      then provide a `length() > 1` character vector with the names of the
#'      columns in `x` that were used for the encoding. If this were encoding
#'      survival this might be `c("time", "event") columns, in that order.
#'   * label: a human readable label to use for this covariate in user facing
#'      scenarios in the facileverse.
#'   * `class`: the "facile class" of the covariate. This can either be
#'     `categorical`, `real`, or `right_censored` (for survival).
#'   * `levels`: (optional) if you want a `categorical` to be treated as a
#'     factor if it isn't already encoded as such in the `pData` itself, or if
#'     you want to rearrange the factor levels.
#'   * `type`: (optinal) this is used a a "grouping" level, particularly in
#'     the FacileExplorer. Not including this won't matter.
#'     TODO: talk about covariate groupings in
#'     `FacileExplorer::categoricalCovariateSelect`
#' @return a list-of-lists that encodes the `sample_covariate` section of the
#'   `meta.yaml` file for a `FacileDataSet`.
#'
#' @examples
#' # covariate_def definition to take tte_OS and tte_event columns and turn
#' # into a facile "OS" right_censored survival covariate
#' cc <- list(
#'   OS=list(
#'     varname=c("tte_OS", "tte_event"),
#'     label="Overall Survival",
#'     class="right_censored",
#'     type="clinical",
#'     description="Overall survival in days"))
create_eav_metadata <- function(x, covariate_def = list()) {
  stopifnot(is.data.frame(x))
  if (is.null(covariate_def)) covariate_def <- list()
  stopifnot(is.list(covariate_def))
  if (length(covariate_def)) {
    validate_covariate_def_list(covariate_def, x)
  }

  if ("dataset" %in% colnames(x)) x[['dataset']] <- NULL
  if ("sample_id" %in% colnames(x)) x[['sample_id']] <- NULL

  # generate generic covariate definitions for all columns
  gcd <- lapply(colnames(x), function(name) eavdef_for_column(x, name))
  names(gcd) <- colnames(x)

  # remove definitions in gcd that are provided in covariate_def
  axe <- lapply(covariate_def, '[[', 'varname')
  axe <- unique(unlist(axe, recursive = TRUE, use.names = FALSE))
  gcd[axe] <- NULL

  out <- c(gcd, covariate_def)
  out
}

#' @rdname eav-metadata
#'
#' @section Defining EAV Encodings:
#'
validate_covariate_def_list <- function(x, pdata) {
  # this is named a list of lists
  stopifnot(is.list(x), is.character(names(x)))
  is.lists <- sapply(x, is.list)
  stopifnot(all(is.lists))
  # the names are unique
  stopifnot(length(unique(names(x))) == length(x))

  # each list item has the following elements
  required.elements <- c('varname', 'class')
  # 'label', 'type' (group), and 'description' are also nice to have, but not
  # strictly required.
  has.elems <- sapply(x, function(el) all(required.elements %in% names(el)))
  stopifnot(all(has.elems))

  # varname entries are valid columns in `pdata`
  for (element in names(x)) {
    cnames <- x[[element]]$varname
    if (!is.character(cnames)) {
      stop(sprintf("%s:varname is not a character", element))
    }
    bad.cnames <- setdiff(cnames, colnames(pdata))
    if (length(bad.cnames)) {
      msg <- "%s:varname element '%s' not found as column in `pdata`"
      msg <- sprintf(msg, element, paste(bad.cnames, collapse="&"))
      stop(msg)
    }
  }
  TRUE
}

#' Generate entity-attribute-value definition for a column in a data.frame
#'
#' Creates the miniimal list-definition for a single column in a `pData`
#' `data.frame`. This function is not exported on purpose.
#'
#' @md
#'
#' @param x a `data.frame`
#' @param column the name of the column in the `x` to parse.
#' @return a generic list-of-list defintion for `x[[column]]`
eavdef_for_column <- function(x, column) {
  vals <- x[[column]]
  if (is.null(vals)) stop("Unknown column in x: ", column)

  out <- list(varname=column, label=column, class='categorical', type="generic",
              description="no description provided")
  if (is.numeric(vals)) {
    out[['class']] <- "real"
  }
  if (is.factor(vals)) {
    out[['levels']] <- levels(vals)
  }
  out
}
