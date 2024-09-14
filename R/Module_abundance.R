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
#' require(ggplot2)
#' #one module
#' moduleframe=Module_abundance(network_obj =network_results,No.module = 3 )
#' moduleframe$rowframe   #combine into metafile
#' moduleframe$columnframe #column table
#' #statistics
#' moduleframe$plotlist$Plotobj_Module3$Statistics
#' #extract plot
#' moduleframe$plotlist$Plotobj_Module3$Barplot
#' moduleframe$plotlist$Plotobj_Module3$Boxplot
#' moduleframe$plotlist$Plotobj_Module3$Violinplot
#'
#'
#' #multiple modules
#' moduleframe=Module_abundance(network_results,c(1,3,6))
#' moduleframe$rowframe
#' moduleframe$columnframe #column table can be used in ggplot visualization
#' #same as above to extract plots and statistics
#' moduleframe$plotlist$Plotobj_Module6$Barplot
Module_abundance=function(
    network_obj,
    No.module
){
  rowframe=network_obj$config$Groupfile
  inputdata=network_obj$config$input_data
  columsum=colSums(inputdata[,-1])
  if(columsum[1]>1){
    inputdata[,-1]=sweep(inputdata[,-1],columsum,"/",MARGIN=2)
  }
  for(i in No.module){
    input_table = as.data.frame(network_obj$Nodes_info)
    nodes_list=input_table$nodes_id[input_table$No.module==i]
    inputtax=network_obj$config$input_taxonomy
    select_table=inputdata[(inputtax[,1])%in% nodes_list ,]
    rowframe=data.frame(rowframe,colSums(select_table[,-1]))
    colnames(rowframe)[ncol(rowframe)]=paste0("Module",i)
  }
  columnframe=gather(rowframe,"Module","Rel",-c(colnames(network_obj$config$Groupfile)))
  outplot=list()
  for(i in unique(columnframe[,"Module"])){
    subdata=columnframe[columnframe[,"Module"]==i,]
    results=compare_plot(inputframe = subdata,treat_location = network_obj$config$taxobj_configuration$treat_location,value_location = ncol(subdata),aes_col = network_obj$config$taxobj_configuration$treat_col,point = TRUE,ylab_text = i)
    outplot=c(outplot,list(results))
    names(outplot)[length(outplot)]=paste0("Plotobj_",i)
  }
  outlist=c(list(rowframe),list(columnframe),list(outplot))
  names(outlist)=c("rowframe","columnframe","plotlist")
  return(outlist)
}
