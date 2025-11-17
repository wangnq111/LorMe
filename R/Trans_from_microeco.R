#' Convert a microeco object to a LorMe object
#'
#' @param microtab  A \code{\link[microeco:microtable]{microtable}} object
#' @param into Names of separated taxonomy to create as character vector. Must select from c("Domain","Phylum","Class","Order","Family","Genus","Species").
#'             Shortcut input:1)By default."standard":c("Domain","Phylum","Class","Order","Family","Genus","Species"). Used for standard taxonomy annotation to OTU/ASV table.
#'                            2)"complete":c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"). Used for complete taxonomy annotation to meta genomic table.

#' @param outputtax Default:c("Phylum","Genus").Names of output taxonomy level table.  Shortcut input is available with 'standard' and 'complete' same as above.
#' @param reads Logical.True for reads table and FALSE for percentage table. Default: TRUE
##'
#'
#' @return LorMe object containing taxonomy table data frame,containing reads and percentage table for each specified output.
#' @export
#'
#' @note
#' For taxonomy annotation with both 'Domain' and 'Kingdom' level, please set 'into' parameter as 'complete'!!!
#'
#' @examples
#' ## Not run:
#'if (requireNamespace("microeco", quietly = TRUE)) {
#'  data("dataset", package = "microeco")
#'
#'  ## convert microeco object to LorMe format
#'  dataset_obj <- Trans_from_microeco(dataset, outputtax = "standard")
#'
#'  ## specify experimental design
#'  dataset_obj_plan <- object_config(dataset_obj,
#'                           treat_location = 4)
#'
#'  ## generate community-structure PCoA plot (Genus level)
#'  community_structure <- structure_plot(dataset_obj_plan, taxlevel = "Genus")
#'  community_structure$PCoA_Plot
#'}
## End(Not run)
Trans_from_microeco <- function(microtab,
                                into      = "standard",
                                outputtax = c("Phylum", "Genus"),
                                reads     = TRUE) {

  if (!requireNamespace("microeco", quietly = TRUE))
    stop("microeco required: install.packages('microeco')")
  if (!inherits(microtab, "microtable"))
    stop("microtab must be a microtable object")

  ## extract
  abund <- microtab$otu_table
  tax   <- microtab$tax_table
  meta  <- microtab$sample_table

  ## transformation
  inputtable <- as.data.frame(abund)
  inputtable <- data.frame(ID = rownames(inputtable), inputtable)

  taxonomytable <- as.data.frame(tax)
  taxonomytable[is.na(taxonomytable)]="Unassigned"
  taxonomytable <- data.frame(ID = rownames(taxonomytable), taxonomy=do.call(paste, c(taxonomytable, sep = ";")))

  ## encaputuration
  out <- tax_summary(groupfile      = as.data.frame(meta),
                     inputtable     = inputtable[, -1],
                     reads          = reads,
                     taxonomytable  = taxonomytable,
                     into           = into,
                     sep            = ";",
                     outputtax      = outputtax)

  ## add tree
  if (!is.null(microtab$phylo_tree)) {
    phylo_tree <- methods::as(microtab$phylo_tree, "phylo")
    out <- add_tree(out, phylo_tree)
  }

  out
}
