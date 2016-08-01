##' Convert +/- encoding of time to event to time + event.
##'
##' Time to event (and event labels) are encoded as a single real value in the
##' clinical data. The absoulte value of the value is the time to the event, and
##' a negative value indicates "event" and a positive value is right censored.
##'
##' @export
##' @param x the time to event
##' @return two column data.frame
decode_right_censored <- function(x, suffix=NULL) {
  x <- as.numeric(x)
  out <- data.frame(tte=abs(x), event=as.integer(x < 0))
  if (is.character(suffix)) {
    names(out) <- paste0(names(out), '_', sub('^[^a-zA-Z]', '', suffix))
  }
  out
}

##' Encodes time and event status into single numeric vector
##'
##' @export
##' @param time \code{numeric} time to event
##' @param event 0/1 vector encoded in the "R sense". "1" is an event, "0" is
##'   right censored.
##' @param sas.encoding default \code{FALSE}. If \code{TRUE} 1 means "censored"
##'   and 0 is an event.
encode_right_censored <- function(time, event, sas.encoding=FALSE) {
  event <- as.integer(event)
  isna <- is.na(event)
  event[isna] <- 1L
  if (!all(event %in% c(0L, 1L))) {
    stop("values in 'event' that are not 0 or 1")
  }
  if (sas.encoding) {
    event <- ifelse(event == 0L, 1L, 0L)
  }
  out <- ifelse(event == 1L, -1L * time, time)
  out[isna] <- NA
  out
}
