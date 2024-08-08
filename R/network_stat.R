####networkstat.Version1.0.0###
###Author:Wang Ningqi####
####度数量 total degree###
#' Igraph network statistics
#'
#' @param input Igraph graph object
#'
#' @return network statistics
#' @importFrom igraph E V degree edge_density diameter transitivity no.clusters centralization.betweenness centralization.degree average.path.length
#' @export
#' @author  Wang Ningqi <2434066068@qq.com>
network_stat=function(input){
  cat("Total degree:",sum(igraph::degree(input)),"\n") #totaldegree=sum(degree(igraph))#
  #  The size of the graph (number of edges)
  cat("Total edges/links:",length(E(input)),"\n")##num.edges = length(E(igraph1)) ##
  #  Order (number of vertices) of a graph
  cat("Total vertices:",length(V(input)),"\n")#num.vertices = length(V(igraph1))#
  # connectance
  cat("Connectance:",edge_density(input,loops=FALSE),"\n")
  # (Average degree) degree/vertices
  cat("Average degree:",mean(igraph::degree(input)),"\n")
  # Diameter
  cat("Diameter:",diameter(input, directed = FALSE, unconnected = TRUE, weights = NULL),"\n")
  # Average path length
  cat("Average path length:",average.path.length(input),"\n")
  # Clustering coefficient
  cat("Global clustering coefficient:",transitivity(input),"\n")
  cat("Number of clusters:",no.clusters(input),"\n")
  # Betweenness centralization
  cat("Betweenness centralization:",centralization.betweenness(input)$centralization,"\n")
  # Degree centralization
  cat("Degree centralization:",centralization.degree(input)$centralization,"\n")}
