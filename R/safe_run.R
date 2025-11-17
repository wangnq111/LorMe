#' Running function with safe mode
#'
#' @param func Functions to run
#'
#' @return Same as of functions
#' @export
#'
#' @examples
#' safe_run(print("Hello"))
safe_run <- function(func) {
  tryCatch(
    func,
    error = function(e) {
      message(conditionMessage(e))
      NULL
    }
  )
}
