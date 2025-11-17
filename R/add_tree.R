#' Attach phylogenetic tree to the LorMe object
#'
#' @param taxobj One LorMe object. Should be generated from \code{\link{tax_summary}}
#' @param tree The phylogenetic tree. Recommanded to load by \code{\link[ape]{read.tree}} or \code{\link[phyloseq]{phy_tree}}
#'
#' @return A LorMe object with phylogenetic tree
#' @export
#'
#'
#' @examples
#' data("Two_group")
add_tree <- function(taxobj, tree) {
  if (!inherits(tree, "phylo")) {
    tree <- try(ape::as.phylo(tree), silent = TRUE)
    if (inherits(tree, "try-error"))
      stop("tree cannot be converted to ape::phylo by as.phylo()")
  }

  tips <- tree$tip.label
  taxa <- methods::slot(taxobj, "data")$Base$ID
  if (!all(tips %in% taxa))
    warning("Tree tip labels not fully matched to feature table")

  methods::new("LorMe",
      groupfile = methods::slot(taxobj, "groupfile"),
      data      = methods::slot(taxobj, "data"),
      configuration    = methods::slot(taxobj, "configuration"),
      tree      = tree) %>% return()
}
