###version1.1.1 ###
###author Wangningqi####
#' Calculate alpha diversity based on tax summary object or dataframe table
#' @description  Calculate alpha diversity of each sample
#'
#' @param taxobj tax summary objects computed by \code{\link{tax_summary}}. Default:NULL.
#' @param taxlevel taxonomy levels used for visualization.Must be one of c("Domain","Phylum","Class","Order","Family","Genus","Species","Base").Default:NULL.
#' @param prefix A character string as prefix of diversity index. Default:""
#' @param input Reads or relative abundance of OTU/Taxa/gene data frame,see details in inputformat. (Useless when taxobj is set).
#' @param inputformat (Useless when taxobj is set)
#' 1:data frame with first column of OTUID and last column of taxonomy
#'
#' 2:data frame with first column of OTUID/taxonomy
#'
#' 3:data frame of all numeric
#' @param reads If the input data frame were from reads table or not(relative abundance table).(Useless when taxobj is set).
#'
#' @return when tax taxobj is set, returns column table with group information combined with for alpha-diversity of each sample,else returns data frame for alpha-diversity of each sample
#' @export
#' @note
#' 1.When input data frame is in relative abundance table,Chao and ACE are not available
#'
#' @importFrom vegan specnumber estimateR
#' @author  Wang Ningqi <2434066068@qq.com>
#' @examples
#' ### Data preparation ####
#' data(testotu)
#' groupinformation <- data.frame(
#'   group = c(rep("a", 10), rep("b", 10)),
#'   factor1 = rnorm(10),
#'   factor2 = rnorm(mean = 100, 10),
#'   subject = factor(c(1:10, 1:10))
#' )
#'
#' # Summary OTU table into genus table and phylum table
#' testtax_summary <- tax_summary(
#'   groupfile = groupinformation,
#'   inputtable = testotu[, 2:21],
#'   reads = TRUE,
#'   taxonomytable = testotu[, c(1, 22)]
#' )
#'
#' ### Use taxsummary object as input ###
#' Alpha <- Alpha_diversity_calculator2(
#'   taxobj = testtax_summary,
#'   taxlevel = "Base"
#' )
#' head(Alpha)
#'
#' # In genus level
#' Alpha <- Alpha_diversity_calculator2(
#'   taxobj = testtax_summary,
#'   taxlevel = "Genus",
#'   prefix = "Genus"
#' )
#' head(Alpha)
#'
#' ### Input dataframe from reads table ###
#' Alpha <- Alpha_diversity_calculator2(
#'   input = testotu,
#'   prefix = "Bacterial",
#'   inputformat = 1,
#'   reads = TRUE
#' )
#'
#' ### Input dataframe from relative abundance table ###
#' if (!require(magrittr)) install.packages("magrittr")
#' library(magrittr)
#' Alpha <- Filter_function(
#'   input = testotu,
#'   threshold = 0,
#'   format = 1
#' ) %>%
#'   Alpha_diversity_calculator2(
#'     input = .,
#'     prefix = "Bacterial",
#'     inputformat = 1,
#'     reads = FALSE
#'   )
#' head(Alpha)
Alpha_diversity_calculator2<- function(taxobj=NULL,taxlevel = NULL,prefix="",input,inputformat,reads){
  if(is.null(taxobj)==TRUE){
    if(inputformat==1){matrix <-as.matrix(t(input[,-c(1,ncol(input))]))}else ##delete annotaion##
      if(inputformat==2){matrix <-as.matrix(t(input[,-1]))}else
        if(inputformat==3){matrix <-as.matrix(t(input))}else ##delete annotaion##
        {warning("ERROR!!!!PLEASE CHOOSE correct inputformat(1,2,3)")}
  }else{
    input=eval(parse(text=paste0("taxobj","$",taxlevel)))
    matrix=as.matrix(t(input[,-1]))
    matrix=round(matrix,0)
    groupframe=taxobj$Groupfile
    reads=TRUE
  }
  shannon<-vegan::diversity(matrix,index='shannon');richness<-specnumber(matrix) ##calculate alpha-diversity##
  evenness<-shannon/log(richness);simpson<-vegan::diversity(matrix,"simpson") ##calculate alpha-diversity##
  if(reads==FALSE){
    alpha.frame<-data.frame(shannon,richness,evenness,simpson) %T>%
      {colnames(.)<-c(paste0(prefix,"Shannon"),paste0(prefix,"Species number"),paste0(prefix,"Simpson"),paste0(prefix,"Evenness"))}}else
        if(reads==TRUE){
          matrix=round(matrix,0)
          chao<-estimateR(matrix)[2,];ACE<-estimateR(matrix)[4,]    ##calculate alpha-diversity##
          alpha.frame<-data.frame(shannon,richness,evenness,simpson,chao, ACE) %T>%
            {colnames(.)<-c(paste0(prefix,"Shannon"),paste0(prefix,"Species number"),paste0(prefix,"Simpson"),paste0(prefix,"Evenness"),paste0(prefix,"Chao"),paste0(prefix,"ACE"))}}
  else {warning("PLEASE CHOOSE reads parameter!!!")}
  if(is.null(taxobj)==FALSE){
    combine_and_translate(inputframe = alpha.frame,groupframe = groupframe,itemname = "Indexname",indexname = "Indexvalue",inputtype = FALSE)%>% return()}else{
      return(alpha.frame)
    }
}
