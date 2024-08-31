###version1.0.1 ###
###author Wangningqi####
#' Combine data for visualization
#' @description Combine group information and index into data frame for visualization(scatter, bar plot, alluvial,box plot etc.).
#' @param inputframe Data frame of index ,sample ID in column,requires all numeric(e.g. result from Alpha_diversity_calculator or Top_taxa function)
#' @param groupframe Data frame of group information(and other abiotic/geographic factors)
#' @param itemname  A character string of your inputframe itemname
#' @param indexname A character string of your inputframe indexname
#' @param inputtype If sample ID were in row and index in column in inputframe.
#'
#' @return key-value pairs data frame
#'
#' @importFrom magrittr %>% %T>%
#' @importFrom tidyr gather
#' @export
#' @author  Wang Ningqi<2434066068@qq.com>
#' @examples
#' {
#'   require(magrittr)
#'   data(testotu)
#'
#'   ## Data preparation ##
#'   Alpha <- Alpha_diversity_calculator2(
#'     input = testotu,
#'     prefix = "Bacterial",
#'     inputformat = 1,
#'     reads = TRUE
#'   )
#'
#'   topotu <- data.frame(
#'     Top_taxa(
#'       input = testotu,
#'       n = 10,
#'       inputformat = 1,
#'       outformat = 1
#'     )[, -1],
#'     row.names = paste0(rep("otu", 11), 1:11)
#'   )
#'
#'   groupinformation1 <- data.frame(
#'     group = c(rep("a", 10), rep("b", 10)),
#'     factor1 = rnorm(10),
#'     factor2 = rnorm(mean = 100, 10)
#'   )
#'
#'   ### Use inputtype = FALSE ###
#'   head(Alpha)
#'   combine_and_translate(
#'     Alpha, groupinformation1,
#'     itemname = "Alpha", indexname = "index",
#'     inputtype = FALSE
#'   )
#'
#'   ### Use inputtype = TRUE ###
#'   head(topotu)
#'   combine_and_translate(
#'     topotu, groupinformation1,
#'     itemname = "OTU", indexname = "reads",
#'     inputtype = TRUE
#'   )
#' }
combine_and_translate<-function(inputframe,groupframe,itemname,indexname,inputtype){
  category<-colnames(groupframe)
  if(inputtype==TRUE){input<-cbind(groupframe,t(inputframe))}else
    if(inputtype==FALSE){input<-cbind(groupframe,inputframe)}
    gather(input,"item","index",-c(category)) %T>%{
      colnames(.)[which(colnames(.)=="item")]=itemname
      colnames(.)[which(colnames(.)=="index")]=indexname} %>%return()}
