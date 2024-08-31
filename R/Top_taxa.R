###version1.2.1 ###
###author Wangningqi####
#' Calculate top taxa and others
#' @description Top taxa is widely used in data analysis,here we provide a simple function to calculate which simplify your R script.
#' @param input Reads or relative abundance(recommended) of OTU/Taxa/gene data frame,see details in inputformat
#' @param n   Top n taxa remained according to relative abundance
#' @param inputformat
#'
#' 1:data frame with first column of OTUID and last column of taxonomy
#'
#' 2:data frame with first column of OTUID/taxonomy (recommended!!!)
#'
#' 3:data frame of all numeric,with row names of OTUID/taxonomy
#' @param outformat
#'
#' 1. return outformat the same as inputformat
#' 2. return data frame of all numeric with OTU/gene/taxa ID in row names(not available for inputformat 1).
#' @return Data frame with top n taxa
#' @export
#' @import magrittr
#' @author  Wang Ningqi<2434066068@qq.com>
#' @examples
#' ### Data preparation ####
#' data(testotu)
#' require(tidyr); require(magrittr)  ## Or use pipe command in "dplyr"
#'
#' testotu.pct <- data.frame(
#'   OTU.ID = testotu[, 1],
#'   sweep(testotu[, -c(1, 22)], 2, colSums(testotu[, -c(1, 22)]), "/"),
#'   taxonomy = testotu[, 22]
#' )
#'
#' sep_testotu <- Filter_function(
#'   input = testotu,
#'   threshold = 0.0001,
#'   format = 1
#' ) %>%
#'   separate(
#'     ., col = taxonomy,
#'     into = c("Domain", "Phylum", "Order", "Family", "Class", "Genus", "Species"),
#'     sep = ";"
#'   )
#'
#' phylum <- aggregate(
#'   sep_testotu[, 2:21], by = list(sep_testotu$Phylum), FUN = sum
#' )
#'
#' phylum1 <- data.frame(row.names = phylum[, 1], phylum[, -1])
#'
#' ##### Input format 1, top 100 OTU #####
#' top100otu <- Top_taxa(
#'   input = testotu.pct,
#'   n = 100,
#'   inputformat = 1,
#'   outformat = 1
#' )
#'
#' ##### Input format 2, top 15 phylum #####
#' head(phylum)
#' top15phylum <- Top_taxa(
#'   input = phylum,
#'   n = 15,
#'   inputformat = 2,
#'   outformat = 1
#' )
#'
#' ##### Input format 3, top 15 phylum #####
#' head(phylum1)
#' top15phylum <- Top_taxa(
#'   input = phylum1,
#'   n = 15,
#'   inputformat = 3,
#'   outformat = 1
#' )

Top_taxa <- function(input, n, inputformat, outformat) {
  if (inputformat == 1) {
    input1 <-
      data.frame(taxonomy = input[, ncol(input)], input[, -c(1, ncol(input))])
  } else
    if (inputformat == 2) {
      input1 <- input
    } else
      if (inputformat == 3) {
        input1 <- data.frame(OTUID = rownames(input), input)
      } else
      {
        warning("Please choose correct inputformat(1,2,3)")
      }
  top_frame <-input1[order(rowMeans(input1[, -1]), decreasing = TRUE), ]  ##order by total relative abundance#
  extra_rm = grep("Unassigned|Unclassified|norank", top_frame[, 1])
  if (length(extra_rm) != 0) {
    left_n = (1:nrow(top_frame))[-extra_rm]
    left_n= left_n[1:n]
    top_frame = top_frame[c(left_n, (1:nrow(top_frame))[-left_n]),]
  }
  top_frame[(n + 1), 1] <- "Others"

  top_frame[(n + 1), -1] <-
    apply(top_frame[(n + 1):nrow(top_frame), -1], 2, sum) ##calculate others#
  top_frame = top_frame[1:(n + 1), ]
  if (inputformat == 3) {
    rownames(top_frame) = top_frame[, 1]
  }
  if (outformat == 1) {
    return(top_frame)
  } else
    if (outformat == 2) {
      rownames(top_frame) <- top_frame[, 1]
      top_frame[, -1] %>% return()
    } else{
      warning("Please choose correct outformat(1/2)!")
    }
} ##remain top taxa###
