#' Conduct Network analysis based on tax summary object
#'
#' @param taxobj tax summary objects computed by \code{\link{tax_summary}}.
#' @param taxlevel taxonomy levels used for analysis. Must be one of c("Domain","Phylum","Class","Order","Family","Genus","Species","Base")
#' @param reads Logical,default:FALSE. Taxonomy abundance type used in analysis.FALSE for relative abundance, TRUE for absolute abundance.
#' @param n Numeric. Number of sample size indicating kept asv/otu/gene/taxa appearing. Recommended to set more than half of total sample size.
#' @param threshold Numeric.Threshold of absolute correlation value (r value for pearson method and rho value for spearman method).
#' @param method Character, default: "spearman". A character indicating which correlation coefficient method to be computed. One of "pearson" or "spearman"
#' @param display Logical, default:TRUE. If display a preview plot of network based on igraph. FALSE for the first attempt is recommended in case of too many vertices and edges.
#' @return  One list contains nodes information table, adjacency column table, adjacency matrix and 'igraph' object.
#' @export
#'
#' @details
#' 1. We had optimized the correlation algorithm to achieve a faster running speed. It takes less than 2 minute to calculate dataframe correlation and p value which more than 400 samples and 10000 OTUs for computer with dual Core i5 processor.
#' However, too many vertices(>2000) or links(>10000) may slow the statistical process and visualization,so we recommend that in your first attempt,set `display` paramter as `F` to have a preview.
#' Then you can adjust your n/threshold/method paramter to generate a suitable visualization network
#'
#' 2. We display a preview plot so as to adjusting your network. Generally a global figure (like we show in examples) with less than 1000 vertices and 5000 edges/links
#' is recommended. Further more,we recommend you to output the statistics and adjacency table and use software like cytoscape or gephi for better visualization.
#'
#' @note
#' 1. Replicates should be at least 5,more than 8 is recommend.
#'
#' 2. In case of too many edges/links or not a global network plot, you can stop the process immediately to provent wasting too much time.
#'
#' @import igraph
#' @importFrom Hmisc rcorr
#' @importFrom dplyr left_join
#' @importFrom fdrtool fdrtool
#' @examples
#' {
#'   ### Data preparation ###
#'   data("Two_group")
#'   set.seed(999)
#'   ## Analysis
#'   network_results <- network_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus",
#'     n = 10,
#'     threshold = 0.8
#'   )
#'
#'   # Nodes information table
#'   network_nodes <- network_results$Nodes_info
#'   head(network_nodes)
#'
#'   # Adjacency table
#'   network_adjacency <- network_results$Adjacency_column_table
#'   head(network_adjacency)
#'
#'   # Complete adjacency matrix
#'   network_matrix <- network_results$Adjacency_matrix
#'   print(network_matrix[1:10, 1:10])
#'
#'   # igraph object
#'   igraph_object <- network_results$Igraph_object
#'   network_stat(igraph_object)  # In case you want to see statistics again
#'   # or do other analysis based on igraph.
#' }
network_analysis<-function(taxobj,taxlevel,reads=FALSE,n,threshold,method="spearman",display=TRUE){
  if(!taxlevel %in% c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species","Base")){
    stop("Illegal 'taxlevel',please choose among c('Domain','Kingdom','Phylum','Class','Order','Family','Genus','Species','Base'")
  }
  if(reads==FALSE){
    input0=eval(parse(text=paste0("taxobj","$",taxlevel,"_percent")))
  }else{
    input0=eval(parse(text=paste0("taxobj","$",taxlevel)))
    }
  taxonomy=eval(parse(text=paste0("taxobj","$",taxlevel,"_taxonomy")))
  input1=data.frame(input0[,-1],row.names =taxonomy[,1])
  zero_count=function(input){length(which(input==0)) %>% return()}
  zerocount=apply(input1,1,zero_count)
  input1=input1[which(zerocount<=(ncol(input1)-n)),]
  input1=input1[which(rowSums(input1)>0),]
  ##network analysis
    corr=rcorr(as.matrix(t(input1)),type=method)
    cor.p=corr$P;cor.p[is.na(cor.p)]<- 0
    fdr=fdrtool(as.numeric(cor.p), statistic="pvalue",plot=FALSE,verbose = FALSE) ##Global fdr correlation##
    cor.r=corr$r
    cor.q=matrix(fdr$qval,ncol=ncol(cor.p),nrow=nrow(cor.p));cor.q[is.nan(cor.q)]<- 0
    cor.r[cor.q>0.05|abs(cor.r)<threshold] = 0 ##fliter via threshold##
    cor.r[is.na(cor.r)]<-0 ##NA变0##
    cor.r[which(cor.r>0.9999)]<-0 ###对角线的1变成0##
    cor.r1=cor.r;cor.r1[which(cor.r1!=0)]<-1 ##无方向矩阵用于igraph##
    cutoff=which(rowSums(cor.r1)== 0)
    if(length(cutoff)==0){cor.r=cor.r;cor.r1=cor.r1}else{
      cor.r=cor.r[-c(cutoff),-c(cutoff)];cor.r1=cor.r1[-c(cutoff),-c(cutoff)]}
    cor.r1[which(cor.r>0)]<-1;cor.r1[which(cor.r<0)]<- -1;
    adj_matrix=cor.r1
    cor.r1[lower.tri(cor.r1,diag = TRUE)]<-NA
  ###creat adjacency###
  source<-rownames(cor.r1)
  adjacency<- data.frame(source,cor.r1) %>%
    gather("target","value",-source)
  adjacency<- adjacency[adjacency$value %in% c(1,-1),]
  adjacency$value=as.numeric(adjacency$value)
  igraph1<-graph_from_data_frame(adjacency,directed = FALSE)##全局输出##
  network_stat(igraph1)
  if(length(E(igraph1))>10000){warning("\n too many edges/links!Better STOP the process")}
  if(isTRUE(display)){
    plot(igraph1,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
         vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))}
  V(igraph1)$degree=igraph::degree(igraph1)

  V(igraph1)$modularity <- membership(cluster_fast_greedy(igraph1))%>%as.numeric()
  cor.tempr=cor.r[order(rownames(cor.r),decreasing = FALSE),];cor.tempr=cor.tempr[,order(colnames(cor.tempr),decreasing = FALSE)]
  #calculate zipi
  communities <- cluster_fast_greedy(igraph1)
  comm_membership <- igraph::membership(communities)
  node_degree <- igraph::degree(igraph1)
  module_mean <- tapply(node_degree, comm_membership, mean)
  module_sd <- tapply(node_degree, comm_membership, sd)
  Zi <- ifelse(module_sd[comm_membership] != 0,
               (node_degree - module_mean[comm_membership]) / module_sd[comm_membership],0)
  module_connections <- matrix(0, nrow=length(V(igraph1)), ncol=max(comm_membership))
  for (v in V(igraph1)) {
    neighbors_of_v <- neighbors(igraph1, v)
    modules_of_neighbors <- comm_membership[neighbors_of_v]
    module_connections[v, ] <- table(factor(modules_of_neighbors, levels=1:max(comm_membership)))
  }
  Pi <- 1 - rowSums((module_connections / node_degree)^2)
  zi_pi8 <- data.frame(
    name = V(igraph1)$name,
    Zi = Zi,
    Pi = Pi
  )
  #merging
  output=data.frame(nodes_id = V(igraph1)$name, node_degree = V(igraph1)$degree, node_betw=betweenness(igraph1),
                    node_evcent=evcent(igraph1,scale = FALSE)$vector,Clustering_coefficient=transitivity(igraph1,type="local"),
                    No.module = V(igraph1)$modularity,Zi=zi_pi8$Zi,Pi=zi_pi8$Pi)
  output=left_join(output,data.frame(nodes_id=taxonomy[,1],taxonomy)) %>% suppressMessages()
  output=output[,-c(9)]
  Groupfile=eval(parse(text=paste0("taxobj","$","Groupfile")))
  config=list(input0,taxonomy,Groupfile,taxlevel,n,threshold,taxobj$configuration)
  names(config)=c("input_data","input_taxonomy","Groupfile","taxlevel","n","threshold","taxobj_configuration")
  outlist=c(list(output),list(adjacency),list(adj_matrix),list(igraph1),list(config))
  names(outlist)=c("Nodes_info","Adjacency_column_table","Adjacency_matrix","Igraph_object","config")
  return(outlist)
}
