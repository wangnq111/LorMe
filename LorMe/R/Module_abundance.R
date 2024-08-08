#' Calculate network module abundance for each sample
#'
#' @param network_obj Network analysis results generated from \code{\link{network_analysis}}
#' @param No.module Numeric or numeric vector of No.module
#'
#' @return A list containing module abundance in metafile and column table of corresponding data frame
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#'
#' @examples#data preparation
#' data("Two_group")
#' ##network analysis
#' network_results<- network_analysis(taxobj = Two_group,taxlevel = "Genus",n = 10,threshold = 0.8)
#'
#' #one module
#' moduleframe=Module_abundance(network_obj =network_results,No.module = 2 )
#' moduleframe$rowframe   #combine into metafile
#' moduleframe$columnframe #column table
#'
#' #multiple modules
#' moduleframe=Module_abundance(network_results,c(2,3,6))
#' moduleframe$rowframe
#' moduleframe$columnframe #column table can be used in ggplot visualization
Module_abundance=function(
    network_obj,
    No.module
){
  rowframe=network_obj$config$Groupfile
  for(i in No.module){
    input_table = as.data.frame(network_obj$Nodes_info)
    nodes_list=input_table$nodes_id[input_table$No.module==i]
    inputtax=network_obj$config$input_taxonomy
    select_table=network_obj$config$input_data[(inputtax[,1])%in% nodes_list ,]
    rowframe=data.frame(rowframe,colSums(select_table[,-1]))
    colnames(rowframe)[ncol(rowframe)]=paste0("Module",i)
  }
  columnframe=gather(rowframe,"Module","Rel",-c(colnames(network_obj$config$Groupfile)))
  outlist=c(list(rowframe),list(columnframe))
  names(outlist)=c("rowframe","columnframe")
  return(outlist)
}
