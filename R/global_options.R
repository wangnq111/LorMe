
#' @title Global options for LorMe pipeline analysis
#' @examples
#' getOption("LorMe")
#' @name .LorMe_opts
#' @export
.LorMe_opts <- list(
  global = list(Analysis_level = "Base", compare_list = NULL),
  alpha  = list(prefix = ""),
  beta   = list(ptsize = 2, diagram = "stick", ellipse.level = 0.85, facet_row = NULL),
  comp   = list(taxlevel = "Phylum", n = 10, palette = "Spectral", nrow = NULL, rmprefix = NULL),
  diff_bar = list(comparison = NULL, rel_threshold = 0, anno_row = "taxonomy", aes_col = NULL, limit_num = 20),
  deseq  = list(comparison = NULL, cutoff = 1, control_name = NULL, paired = FALSE, subject = NULL),
  indic  = list(func = "r.g", reads = FALSE),
  manh   = list(taxlevel = "Phylum", controlname = NULL, mode = "most", top_n = 10, palette = NULL, select_tax = NULL, rmprefix = NULL),
  sub_net = list(reads = FALSE, n = NULL, threshold = 0.9, rel_threshold = 0, method = "spearman", display = TRUE),
  all_net = list(reads = FALSE, n = NULL, threshold = 0.9, rel_threshold = 0, method = "spearman", display = TRUE),
  metanet = list(tag_threshold = 5)
)


#' Accessed by  `getOption("LorMe")` or `LorMe_options()`
#'
#' @name .onLoad
NULL
.onLoad <- function(libname, pkgname) options(LorMe = .LorMe_opts)
