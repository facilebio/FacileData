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
encode_right_censored <- function(time, event) {
  ifelse(event == 1, -1 * time, time)
}
