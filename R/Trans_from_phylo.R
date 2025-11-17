#' Convert phyloseq object to LorMe object
#'
#' @param physeq The phyloseq object to convert. See in \code{\link[phyloseq:phyloseq]{phyloseq}} object
#' @param into Names of separated taxonomy to create as character vector. Must select from c("Domain","Phylum","Class","Order","Family","Genus","Species").
#'             Shortcut input:1)By default."standard":c("Domain","Phylum","Class","Order","Family","Genus","Species"). Used for standard taxonomy annotation to OTU/ASV table.
#'                            2)"complete":c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"). Used for complete taxonomy annotation to meta genomic table.

#' @param outputtax Default:c("Phylum","Genus").Names of output taxonomy level table.  Shortcut input is available with 'standard' and 'complete' same as above.
#' @param reads Logical.True for reads table and FALSE for percentage table. Default: TRUE
#'
#' @return LorMe object containing taxonomy table data frame,containing reads and percentage table for each specified output.
#' @export
#'
#' @note
#' For taxonomy annotation with both 'Domain' and 'Kingdom' level, please set 'into' parameter as 'complete'!!!
#'
#' The taxonomyTable of phyloseq object should contained as least standard annotations for convert!
#' @examples
#' ## Not run:
#'\donttest{
#'if (requireNamespace("phyloseq", quietly = TRUE)) {
#'  data("GlobalPatterns", package = "phyloseq")
#'
#'  ## convert phyloseq object to LorMe format
#'  gp_obj <- Trans_from_phylo(GlobalPatterns, outputtax = "standard")
#'
#'  ## specify experimental design
#'  gp_plan <- object_config(gp_obj,
#'                           treat_location = 6,
#'                           rep_location   = 7)
#'
#'  ## generate community-structure PCoA plot (Genus level)
#'  community_structure <- structure_plot(gp_plan, taxlevel = "Genus")
#'  community_structure$PCoA_Plot
#'}
#'}
## End(Not run)
#'
Trans_from_phylo <- function(physeq,
                             into      = "standard",
                             outputtax = c("Phylum", "Genus"),
                             reads     = TRUE) {

  if (!requireNamespace("phyloseq", quietly = TRUE))
    stop("phyloseq required: install.packages('phyloseq')")
  if (!inherits(physeq, "phyloseq"))
    stop("'physeq' must be a phyloseq object")
  #extract from phyloseq
  otu  <- phyloseq::otu_table(physeq)
  tax  <- phyloseq::tax_table(physeq)
  meta <- phyloseq::sample_data(physeq)

  ## data transformation
  inputtable <- as.data.frame(otu)
  inputtable <- data.frame(ID = rownames(inputtable), inputtable)

  taxonomytable <- as.data.frame(tax)
  taxonomytable[is.na(taxonomytable)]="Unassigned"
  taxonomytable <- data.frame(ID = rownames(taxonomytable), taxonomy=do.call(paste, c(taxonomytable, sep = ";")))

  ## tree data
  phylo_tree <- if (!is.null(methods::slot(physeq, "phy_tree"))) {
    methods::as(phyloseq::phy_tree(physeq), "phylo")
  } else {
    NULL
  }

  ## encapsulation
  out <- tax_summary(groupfile      = as.data.frame(meta),
                     inputtable     = inputtable[ , -1],
                     reads          = reads,
                     taxonomytable  = taxonomytable,
                     into           = into,
                     sep            = ";",
                     outputtax      = outputtax)

  ## add tree
  if (!is.null(phylo_tree)) {
    out <- add_tree(out, phylo_tree)
  }
  return(out)
}
