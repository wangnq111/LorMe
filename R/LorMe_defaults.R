#' Restore LorMe global options to factory defaults
#'
#' Equivalent to `options(LorMe = .LorMe_opts)` but provides a unified entry.
#'
#' @return Invisible NULL
#' @examples
#' LorMe_defaults()   # one-click factory reset
#' @export
LorMe_defaults <- function() {
  options(LorMe = .LorMe_opts)
  invisible(NULL)
}
