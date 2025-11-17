#' @importFrom methods setOldClass setClassUnion setClass setMethod getMethod "slot" "slotNames"
NULL

#（移除原来检测 baseenv() 中 DollarNames 的那段无效代码）
# 注意：R 的补全机制使用的是 utils::.DollarNames；我们在 .onLoad 注册 S3 方法以确保补全生效。

if (requireNamespace("ape", quietly = TRUE)) {

  setOldClass(c("phylo", "ape"))
  setClassUnion("phyloOrNULL", c("phylo", "NULL"))
} else {
  setClassUnion("phyloOrNULL", "NULL")
}

#' LorMe S4 class
#'
#' A container class for microbial ecology data integration, including metadata,
#' feature tables, configuration parameters, and phylogenetic tree.
#'
#' @slot groupfile data.frame: Sample metadata
#' @slot data list: Named list containing feature and taxonomy tables
#' @slot configuration list: Analysis configuration parameters
#' @slot tree phyloOrNULL: Phylogenetic tree (NULL if not provided, requires ape package)
#'
#' @section Methods:
#' \describe{
#'   \item{\code{show(object)}}{Print summary of the LorMe object}
#'   \item{\code{x$name}}{Access slots or data elements via $}
#'   \item{\code{.DollarNames(x, pattern)}}{Tab completion for $ access}
#' }
#'
#' @name LorMe-class
#' @aliases LorMe-class
#' @exportClass LorMe
setClass("LorMe",
         slots = c(
           groupfile = "data.frame",
           data = "list",
           configuration = "list",
           tree = "phyloOrNULL"
         ),
         prototype = list(tree = NULL))

#' Access LorMe object elements via $
#'
#' S4 method to enable $ operator for LorMe objects, supporting both slots and
#' elements in the 'data' slot.
#'
#' @param x A LorMe object
#' @param name Character: Slot name or key in the 'data' slot
#' @return Value of the specified slot or data element
#' @exportMethod $
#' @export
setMethod("$", "LorMe", function(x, name) {
  if (name %in% methods::slotNames(x)) {
    return(methods::slot(x, name))
  } else if (name %in% names(slot(x, "data"))) {
    return(methods::slot(x, "data")[[name]])
  } else {
    stop(sprintf("\"%s\" is not a valid key for LorMe", name))
  }
})

#' Tab completion for LorMe $ access
#'
#' enable tab completion when using $ with LorMe objects, including
#' both slot names and keys in the 'data' slot.
#'
#' @param x A LorMe object
#' @param pattern Character: Prefix pattern for matching
#' @return Character vector of valid keys matching the pattern
#' @importFrom utils .DollarNames
#' @export
.DollarNames.LorMe <- function(x, pattern = "") {
  all_keys <- c(methods::slotNames(x), names(methods::slot(x, "data")))
  unique(all_keys[startsWith(all_keys, pattern)])
}

# Register the S3 method for utils::.DollarNames when the package is loaded.
# This ensures RStudio / utils tab-completion will call .DollarNames.LorMe.
.onLoad <- function(libname, pkgname) {
  # registerS3method is the recommended runtime registration approach
  # (works even without editing NAMESPACE). utils is a base namespace, so this is safe.
  if (exists("registerS3method", envir = asNamespace("base"), inherits = FALSE)) {
    registerS3method(".DollarNames", "LorMe", .DollarNames.LorMe, envir = asNamespace("utils"))
  } else {
    # fallback: try direct registration (very unlikely to be needed)
    assign(".DollarNames.LorMe", .DollarNames.LorMe, envir = asNamespace("utils"))
  }
  invisible()
}

#' Show method for LorMe objects
#'
#' Custom summary display when LorMe object is printed.
#'
#' @param object A LorMe object
#' @aliases show,LorMe-method
#' @importFrom methods show
#' @exportMethod show
setMethod("show", "LorMe", function(object) {
  cat("LorMe object\n")
  cat(" Groupfile rows:", nrow(object@groupfile), "\n")
  cat(" Data keys:", paste(names(object@data), collapse = ", "), "\n")
  cat(" Configuration keys:", paste(names(object@configuration), collapse = ", "), "\n")
  cat(" Tree:", if (is.null(object@tree)) "NULL" else "phylo object (requires ape)", "\n")
})
