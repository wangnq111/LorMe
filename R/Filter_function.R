###version1.2.1 ###
###author Wangningqi####
#'  Filter OTU/ASV/metagenomic profile/gene profile by threshold
#' @description Sequenced data of taxonomy&gene still remains some sequencing error which we needed to be wiped off before analyzing.
#'     Here we provide function including four formats to wipe them clean.
#' @param input Data frame of absolute abundance of standard OTU table,with the first column of OTUID and the final column of taxonomy annotation.
#'              If your data frame is gene table or not a standard OTU table, please manually transformed into a standard input data frame.
#' @param threshold threshold of filter.Relative abundance for format 1 and 4, reads number for format 2, sample size for format 3
#' @param format
#' 1:filter OTU/gene below overall-sample relative abundance threshold(<)
#'
#' 2:filter OTU/gene below overall-sample reads threshold(<)
#'
#' 3:filter OTU/gene reads 0 over threshold sample size(>)
#'
#' 4:filter OTU/gene below relative abundance threshold in each sample(<)
#' @param report Logical. If print report to console. Default:TRUE
#' @importFrom magrittr %>% %T>%
#' @return Dataframe of OTU/gene in format of absolute abundacne(reads) or relative abundance(%)
#' @export
#' @importFrom magrittr %>%
#' @author  Wang Ningqi
#' @examples
#' ### Data frame with absolute abundance (reads)###
#' ### And first column of OTUID and last column of taxonomy ###
#' data(testotu)
#'
#' #### If your data frame does not contain the OTUID column or taxonomy column,
#' #### you can add a simulated column to fit the input format like testotu ##
#'
#' ### 1. Filter OTU with total relative abundance below 0.0001, return absolute abundance ###
#' filtered_otu <- Filter_function(
#'   input = testotu,
#'   threshold = 0.0001,
#'   format = 1
#' )
#'
#' ### 2. Filter OTU with total relative abundance below 0.0001, return relative abundance ###
#' filtered_otu <- Filter_function(
#'   input = testotu,
#'   threshold = 0.0001,
#'   format = 1
#' )
#'
#' ### 3. Filter OTU with total reads below 20 ###
#' filtered_otu <- Filter_function(
#'   input = testotu,
#'   threshold = 20,
#'   format = 2
#' )
#'
#' ### 4. Filter OTU reads 0 over (>=) 11 samples ###
#' filtered_otu <- Filter_function(
#'   input = testotu,
#'   threshold = 11,
#'   format = 3
#' )
#'
#' ### 5. Filter OTU with relative abundance below 0.0001 in each sample ###
#' filtered_otu <- Filter_function(
#'   input = testotu,
#'   threshold = 0.0001,
#'   format = 4
#' )
Filter_function<-function(input,threshold,format,report=TRUE){
  format_percent <- function(input){sweep(input,2,colSums(input),'/') %>% return()}#Fun:reads2percent#
  input_percentage<-cbind(input[,1],format_percent(input[,-c(1,ncol(input))]),input[,ncol(input)]) %T>%
    {colnames(.)[1]=colnames(input)[1];colnames(.)[ncol(input)]=colnames(input)[ncol(input)]}
  if(format==1){
    inputmeans<-rowMeans(input_percentage[,-c(1,ncol(input_percentage))]) ##form PCT data & relative index###
    output=input[which(inputmeans>threshold),]
    }else
  if (format==2){
    output <- which(rowSums(input[,-c(1,ncol(input))])>=threshold) %>%
   input[.,]}else
  if(format==3){
    zero_count=function(input){length(which(input==0)) %>% return()}  ##Fun:calculate zero count ##
    zerocount=apply(input[,-c(1,ncol(input))],1,zero_count)
    output <- input[which(zerocount<=threshold),]}else
  if (format==4){
    low_count <- function(input){length(which(input<threshold))%>% return()}##Fun:calculate low count
    lowcount=apply(input_percentage[,-c(1,ncol(input_percentage))],1,low_count)
    output <- input[which(lowcount!=(ncol(input)-2)),]}
  else{stop("Errorrrrrrrrrrrrr!!!Wrong format!!!!")}
  if(report==TRUE){
  cat("#Filtering report \n#overview\n")
  report=(colSums(output[,2:(ncol(output)-1)])/colSums(input[,2:(ncol(output)-1)]))
  print(report)
  cat("#Summary\n")
  summary(report)%>%print()
  }
  return(output)
  }
