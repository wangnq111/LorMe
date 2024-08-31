#' Indicator Analysis
#' @description Performs the indicator analysis based on taxonomic summary object
#' @param taxobj tax summary objects computed by \code{\link{tax_summary}}. Default:NULL.
#' @param taxlevel taxonomy levels used for visualization.Must be one of c("Domain","Phylum","Class","Order","Family","Genus","Species","Base").Default:NULL.
#' @param func Default: "r.g".The function to use for the indicator analysis, see in \code{\link[indicspecies]{multipatt}}
#' @param reads A logical value indicating whether the input data is in terms of raw reads (TRUE) or relative abundance (FALSE)
#'
#' @return A data frame with the results of the indicator analysis, including adjusted p-values, tags and taxonomic information.
#' @export
#'
#' @importFrom dplyr left_join
#' @note
#' This function depends on the following packages: indicspecies, permute.
#' These packages are not automatically loaded and should be installed before using this function.
#' @seealso \code{\link[indicspecies]{multipatt}}, \code{\link[permute]{how}}
#' @examples
#' data("Two_group")
#' if (requireNamespace("indicspecies", quietly = TRUE) &&
#'     requireNamespace("permute", quietly = TRUE)) {
#'   set.seed(999)
#'   indicator_results <- indicator_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus"
#'   )
#'   head(indicator_results)
#' }
indicator_analysis=function(taxobj,taxlevel,func="r.g",reads=FALSE){
  if (!requireNamespace("indicspecies", quietly = TRUE)) {
    stop("The 'indicspecies' package is required for this function. Please install it using install.packages('indicspecies').")
  }
  if (!requireNamespace("permute", quietly = TRUE)) {
    stop("The 'permute' package is required for this function. Please install it using install.packages('permute').")
  }
  if(is.null(taxobj$configuration)){
    stop("taxonomic summary object not configured yet, call '?object_config' for configuration!")
    return(NULL)
  }
  if(is.null(eval(parse(text=paste0("taxobj","$",taxlevel))))){
    stop("Illegal 'taxlevel'!")
    return(NULL)
  }
  groupfile= eval(parse(text=paste0("taxobj$Groupfile")))
  condition= groupfile[eval(parse(text=paste0("taxobj$configuration$treat_location"))) ]
  condition= condition[,1]
  if(length(unique(condition))==1){
    stop("Insufficient group for comparsion!")
    return(NULL)
  }
  taxonomy=eval(parse(text=paste0("taxobj","$",taxlevel,"_taxonomy")))
  if(reads==FALSE){
    input0=eval(parse(text=paste0("taxobj","$",taxlevel,"_percent")))
  }else{
    input0=eval(parse(text=paste0("taxobj","$",taxlevel)))
  }
  inputframe=data.frame(input0[,-1],row.names =taxonomy[,1])
  message("Calculating.......\n")
  indicator_results=indicspecies::multipatt(as.data.frame(t(inputframe)),condition,func = func,control=permute::how(nperm=1000)) %$%
    as.data.frame(sign)
  indicator_results$ID=rownames(indicator_results)
  colnames(indicator_results)[ncol(indicator_results)]=colnames(taxonomy)[1]
  indicator_results_anno=left_join(indicator_results,taxonomy) %>% suppressMessages()
  indicator_results_anno$tag="None"
  if(length(unique(condition))==2){
    indicator_results_anno$tag[indicator_results_anno[,paste0("s.",unique(condition)[1])]==1&indicator_results_anno$p.value<0.05]=unique(condition)[1]
    indicator_results_anno$tag[indicator_results_anno[,paste0("s.",unique(condition)[2])]==1&indicator_results_anno$p.value<0.05]=unique(condition)[2]
  }else if(length(unique(condition))==3){
    if(is.null(taxobj$configuration$treat_order)){
      treat_order=unique(condition) %>% sort()
    }else{
      treat_order=taxobj$configuration$treat_order
    }
    indicator_results_anno$tag[indicator_results_anno[,paste0("s.",treat_order[1])]==1&indicator_results_anno[,paste0("s.",treat_order[2])]==0&indicator_results_anno[,paste0("s.",treat_order[3])]==0&indicator_results_anno$p.value<0.05]=treat_order[1]
    indicator_results_anno$tag[indicator_results_anno[,paste0("s.",treat_order[1])]==0&indicator_results_anno[,paste0("s.",treat_order[2])]==1&indicator_results_anno[,paste0("s.",treat_order[3])]==0&indicator_results_anno$p.value<0.05]=treat_order[2]
    indicator_results_anno$tag[indicator_results_anno[,paste0("s.",treat_order[1])]==0&indicator_results_anno[,paste0("s.",treat_order[2])]==0&indicator_results_anno[,paste0("s.",treat_order[3])]==1&indicator_results_anno$p.value<0.05]=treat_order[3]
    indicator_results_anno$tag[indicator_results_anno[,paste0("s.",treat_order[1])]==1&indicator_results_anno[,paste0("s.",treat_order[2])]==1&indicator_results_anno[,paste0("s.",treat_order[3])]==0&indicator_results_anno$p.value<0.05]=paste0(treat_order[1],"_",treat_order[2])
    indicator_results_anno$tag[indicator_results_anno[,paste0("s.",treat_order[1])]==1&indicator_results_anno[,paste0("s.",treat_order[2])]==0&indicator_results_anno[,paste0("s.",treat_order[3])]==1&indicator_results_anno$p.value<0.05]=paste0(treat_order[1],"_",treat_order[3])
    indicator_results_anno$tag[indicator_results_anno[,paste0("s.",treat_order[1])]==0&indicator_results_anno[,paste0("s.",treat_order[2])]==1&indicator_results_anno[,paste0("s.",treat_order[3])]==1&indicator_results_anno$p.value<0.05]=paste0(treat_order[2],"_",treat_order[3])
    }else{
      for(i in unique(condition)){
        indicator_results_anno$tag[indicator_results_anno[,paste0("s.",i)]==1&rowSums(indicator_results_anno[,1:length(unique(condition))]==1)&indicator_results_anno$p.value<0.05]=i
      }
    }
  colnames(indicator_results_anno)[which(colnames(indicator_results_anno)=="p.value")]="padj"
  indicator_results_anno$baseMean=rowMeans(inputframe)
  return(indicator_results_anno)
}
