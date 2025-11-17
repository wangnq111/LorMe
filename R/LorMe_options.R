#' @title Change LorMe options
#' @param ... Named parameters or "level = list" pairs.
#' @param level Optional; specifies a whole-layer replacement.
#' @return Modified global options (or invisible(NULL) after factory reset).
#' @importFrom utils modifyList
#' @examples
#' LorMe_options()                                    # view current
#' LorMe_options(beta = list(diagram = "circle"))     # whole layer
#' LorMe_options(level = "beta", list(diagram = "circle"))  # within layer
#' @export
LorMe_options <- function(..., level = NULL) {
  opts  <- list(...)
  if (length(opts) == 0L) return(getOption("LorMe"))

  old   <- getOption("LorMe")
  nms   <- names(old)          #
  ###
  if (length(opts) == 1L && !is.null(names(opts)) && names(opts) %in% nms) {
    lvl <- names(opts)
    old[[lvl]] <- opts[[lvl]]
    return(options(LorMe = old))
  }
  ###
  for (lvl in names(opts)) {
    if (!lvl %in% nms)
      stop(sprintf("Level '%s' does not exist. Valid levels: %s",
                   lvl, paste(nms, collapse = ", ")))
    old[[lvl]] <- modifyList(old[[lvl]], opts[[lvl]])
  }
  options(LorMe = old)
}

