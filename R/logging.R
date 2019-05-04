#' Generates a logging message using glue and crayon, with some bells/whistles.
#'
#' Like other logging approaches, each message created with this function is
#' assigned a `level` (priority). If the current logging level, which is
#' returned from a call to `flog_level` (ostensibly determenied by the value
#' of the `"facile.log.level(.*?)"` option) is less than or equal to level of
#' this message, then the message will be generated and sent to `file`.
#' You can include a `namespace` for the message to provide a namespace-specific
#' level/priority hierarchy.
#'
#' Conveninece wrapper functions are provided for each logging level, ie.
#' call `fwarn("message")` instead of flog("message", level = "warn")`. Also,
#' each facile* package provides its own `flog()` function which sets the
#' namespace `ns` parameter to default to a package-specific namespace so you
#' can control logging at the different package level.
#'
#' @sectio Logging Levels:
#' Logging levels are:
#'
#' ```
#' .flog_levels <- c("all"  = 0, "trace" = 1, "debug" = 2, "info" = 3,
#'                   "warn" = 4, "error" = 5, "fatal" = 6)
#' ```
#'
#' @section crayon:
#' Glue lets you put cayon functions in `{}` to stylize output. For instance,
#' you can make "bold and red" the color red and also bold, like so:
#'
#' ```r
#' flog("This is {red}{bold}bold and red{reset}, right?")
#' ```
#'
#' Nice! It would be *cooler* if we could make it a bit more terse, like so:
#'
#' ```r
#' flog("This is rb`bold and red`, right?")
#' ```
#'
#' Colors:
#'
#' * b: blue
#' * c: cyan
#' * g: green
#' * k: black
#' * m: magenta
#' * r: red
#' * y: yellow
#'
#' Styles:
#'
#' * i: italic
#' * s: strong (bold)
#' * S: striketthrough
#' * u: underline
#'
#' @export
#' @importFrom glue glue
#' @importFrom crayon blue cyan green black magenta red yellow
#' @importFrom crayon bold italic strikethrough underline
#'
#' @param ... the string elements to pass into `glue::glue()`
#' @param level the "firing level" of this message. Defaults to "info"
#' @param ns (namespace) if included, then the message checks the
#'   namespace-specific logging priority
#' @param session,file,sep,fill,labels,append sent to [base::cat()]
#' @param newline If `TRUE`, appends a `\n` to the message. By default, this
#'   is `TRUE` when `file` is not `NULL`.
#' @return invisibly returns the text generated in the logging message.
flog <- function(..., level = "info", ns = NULL, session = NULL,
                 file = stderr(), sep = "", fill = FALSE, labels = NULL,
                 append = FALSE, newline = !is.null(file)) {
  level <- assert_choice(level, names(.flog_levels))
  level <- .flog_levels[level]
  if (level < flog_level(ns)) return(invisible(NULL))

  reset <- crayon::reset
  if (missing(session)) {
    session <- try(get("session", envir = parent.frame()), silent = TRUE)
  }
  if (is(session, "session_proxy")) {
    txt <- glue("{bold}[{smod}]:{reset} ", smod = session$ns(NULL))
  } else {
    txt <- ""
  }
  txt <- glue(txt, do.call(glue, list(...)))
  if (newline) txt <- paste0(txt, "\n")
  if (!is.null(file)) {
    cat(txt, file = file, sep = sep, fill = fill, labels = labels,
        append = append)
  }
  invisible(txt)
}

# Default level is "warn"
.flog_levels <- c("all"  = 0, "trace" = 1, "debug" = 2, "info" = 3,
                  "warn" = 4, "error" = 5, "fatal" = 6)

#' Retrieves the currently set logging level
#'
#' @export
#' @param namespace Package (or whoever) can provide a value here to set the
#'   level they want to listen to. If this is `NULL` (default), the top level
#'   `facile.log.level` value will be used.
#' @return the logging level, as an integer (from `FacileData:::.flog_levels`)
flog_level <- function(namespace = NULL) {
  optkey <- "facile.log.level"
  if (!is.null(namespace)) {
    assert_string(namespace)
    optkey <- paste0(optkey, ".", namespace)
  }

  lvl <- getOption(optkey, "warn")
  priority <- .flog_levels[lvl]
  if (is.na(priority)) priority <- .flog_levels["warn"]
  priority
}

#' @noRd
#' @export
ftrace <- function(..., ns = NULL, session = NULL) {
  if (missing(session)) {
    session <- try(get("session", envir = parent.frame()), silent = TRUE)
  }
  flog(..., level = "trace", ns = ns, session = session)
}

#' @noRd
#' @export
fdebug <- function(..., ns = NULL, session = NULL) {
  if (missing(session)) {
    session <- try(get("session", envir = parent.frame()), silent = TRUE)
  }
  flog(..., level = "debug", ns = ns, session = session)
}

#' @noRd
#' @export
finfo <- function(..., ns = NULL, session = NULL) {
  if (missing(session)) {
    session <- try(get("session", envir = parent.frame()), silent = TRUE)
  }
  flog(..., level = "info", ns = ns, session = session)
}

#' @noRd
#' @export
fwarn <- function(..., ns = NULL, session = NULL) {
  if (missing(session)) {
    session <- try(get("session", envir = parent.frame()), silent = TRUE)
  }
  flog(..., level = "warn", ns = ns, session = session)
}

#' @noRd
#' @export
ferror <- function(..., ns = NULL, session = NULL) {
  if (missing(session)) {
    session <- try(get("session", envir = parent.frame()), silent = TRUE)
  }
  flog(..., level = "error", ns = ns, session = session)
}

#' @noRd
#' @export
ffatal <- function(..., ns = NULL, session = NULL) {
  if (missing(session)) {
    session <- try(get("session", envir = parent.frame()), silent = TRUE)
  }
  flog(..., level = "fatal", ns = ns, session = session)
}
