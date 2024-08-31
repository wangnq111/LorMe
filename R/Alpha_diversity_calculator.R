#' Calculate alpha diversity based on tax summary object
#' @description  Calculate alpha diversity for each sample
#'
#' @param taxobj tax summary objects computed by \code{\link{tax_summary}}. Default:NULL.
#' @param taxlevel taxonomy levels used for visualization.Must be one of c("Domain","Phylum","Class","Order","Family","Genus","Species","Base").Default:NULL.
#' @param prefix A character string as prefix of diversity index. Default:""
#'
#' @return 'Alpha_diversity_calculator' returns alpha-diversity of each sample in format of column table (dataframe) combined with group information in meta file.
#' @export
#'
#' @importFrom vegan specnumber estimateR
#' @examples
#' ###data preparation####
#' data("Two_group")
#' require(ggplot2)
#'
#' ###analysis####
#' Alpha_results<- Alpha_diversity_calculator(taxobj = Two_group,taxlevel = "Base")
#'
#' #Check data frame contained alpha diversity
#' head(Alpha_results$alphaframe,5)
#'
#' #Check contained statistics and plot list
#' names(Alpha_results$plotlist)
#'
#' #Check statistics for Shannon
#' Alpha_results$plotlist$Plotobj_Shannon$Statistics
#'
#' #Extract plot for Shannon
#' Alpha_results$plotlist$Plotobj_Shannon$Barplot
#' Alpha_results$plotlist$Plotobj_Shannon$Boxplot
#' Alpha_results$plotlist$Plotobj_Shannon$Violinplot
#'
#'
Alpha_diversity_calculator<- function(taxobj,taxlevel,prefix=""){
  if(is.null(eval(parse(text=paste0("taxobj","$",taxlevel))))){
    warning("Illegal 'taxlevel'!")
    return(NULL)
  }
  if(is.null(taxobj["configuration"])){warning("taxonomic summary object not configured yet, call '?object_config' for configuration")}
  input=eval(parse(text=paste0("taxobj","$",taxlevel)))
  matrix=as.matrix(t(input[,-1]))
  matrix=round(matrix,0)
  groupframe=taxobj$Groupfile
  shannon<-vegan::diversity(matrix,index='shannon');richness<-specnumber(matrix) ##calculate alpha-diversity##
  evenness<-shannon/log(richness);simpson<-vegan::diversity(matrix,"simpson") ##calculate alpha-diversity##
  matrix=round(matrix,0)
  chao<-estimateR(matrix)[2,];ACE<-estimateR(matrix)[4,]    ##calculate alpha-diversity##
  alpha.frame<-data.frame(shannon,richness,evenness,simpson,chao, ACE) %T>%
    {colnames(.)<-c(paste0(prefix,"Shannon"),paste0(prefix,"Species number"),paste0(prefix,"Simpson"),paste0(prefix,"Evenness"),paste0(prefix,"Chao"),paste0(prefix,"ACE"))}
  alpha.frame<-combine_and_translate(inputframe = alpha.frame,groupframe = groupframe,itemname = "Indexname",indexname = "Indexvalue",inputtype = FALSE)
  outplot=list()
  for(i in unique(alpha.frame[,"Indexname"])){
    subdata=alpha.frame[alpha.frame[,"Indexname"]==i,]
    results=compare_plot(inputframe = subdata,treat_location = taxobj$configuration$treat_location,value_location = ncol(subdata),aes_col = taxobj$configuration$treat_col,point = TRUE,ylab_text = i)
    outplot=c(outplot,list(results))
    names(outplot)[length(outplot)]=paste0("Plotobj_",i)
  }
  outlist=c(list(alpha.frame),list(outplot))
  names(outlist)=c("alphaframe","plotlist")
  return(outlist)
}
