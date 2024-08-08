#' Calculate alpha diversity based on tax summary object
#' @description  Calculate alpha diversity for each sample
#'
#' @param taxobj tax summary objects computed by \code{\link{tax_summary}}. Default:NULL.
#' @param taxlevel taxonomy levels used for visualization.Must be one of c("Domain","Phylum","Class","Order","Family","Genus","Species","Base").Default:NULL.
#' @param prefix A character string as prefix of diversity index. Default:""
#'
#' @return 'Alpha_diversity_calculator' returns α-diversity of each sample in format of column table (dataframe) combined with group information in meta file.
#' @export
#'
#' @importFrom vegan specnumber estimateR
#' @examples
#' ###data prepration####
#' data("Two_group")
#'
#' ###analysis####
#'Alpha<- Alpha_diversity_calculator(taxobj = Two_group,taxlevel = "Base")
#'head(Alpha)
#'
#'Alpha<- Alpha_diversity_calculator(taxobj = Two_group,taxlevel = "Genus")
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
  shannon<-vegan::diversity(matrix,index='shannon');richness<-specnumber(matrix) ##calculate α-diversity##
  evenness<-shannon/log(richness);simpson<-vegan::diversity(matrix,"simpson") ##calculate α-diversity##
  matrix=round(matrix,0)
  chao<-estimateR(matrix)[2,];ACE<-estimateR(matrix)[4,]    ##calculate α-diversity##
  alpha.frame<-data.frame(shannon,richness,evenness,simpson,chao, ACE) %T>%
    {colnames(.)<-c(paste0(prefix,"Shannon"),paste0(prefix,"Species number"),paste0(prefix,"Simpson"),paste0(prefix,"Evenness"),paste0(prefix,"Chao"),paste0(prefix,"ACE"))}
  combine_and_translate(inputframe = alpha.frame,groupframe = groupframe,itemname = "Indexname",indexname = "Indexvalue",inputtype = F)%>% return()
}
