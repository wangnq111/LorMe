####networkstat.Version1.0.0###
###Author:Wang Ningqi####
####度数量 total degree###
#' Igraph network statistics
#'
#' @param input An igraph object.
#' @param report Logical. If print report to console. Default:TRUE
#'
#' @return network statistics
#' @importFrom igraph E V degree edge_density diameter transitivity no.clusters centralization.betweenness centralization.degree average.path.length
#' @export
#' @author  Wang Ningqi <2434066068@qq.com>
network_stat=function(input,report=TRUE){
  if (!inherits(input, "igraph")) {
    stop("The input must be an igraph object.")
  }
  if (report==TRUE) {
    message("Total degree: ", sum(igraph::degree(input)))
    message("Total edges/links: ", length(E(input)))
    message("Total vertices: ", length(V(input)))
    message("Connectance: ", edge_density(input, loops = FALSE))
    message("Average degree: ", mean(igraph::degree(input)))
    message("Diameter: ", diameter(input, directed = FALSE, unconnected = TRUE, weights = NULL))
    message("Average path length: ", average.path.length(input))
    message("Global clustering coefficient: ", transitivity(input))
    message("Number of clusters: ", no.clusters(input))
    message("Betweenness centralization: ", centralization.betweenness(input)$centralization)
    message("Degree centralization: ", centralization.degree(input)$centralization)
  }
  }
