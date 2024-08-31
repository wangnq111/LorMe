###version2.0.0 ###
###author Wangningqi####
#' Deseq Analysis Function
#'
#' This function performs a differential expression analysis using the DESeq2 package.
#' It is designed to work with microbiome data and can handle paired or non-paired samples.
#'
#' @param taxobj Taxonomic summary objects containing the taxonomic information, relative abundance, metafile and configuration.
#' @param taxlevel The taxonomic level for the analysis.Must be one of c("Domain","Phylum","Class","Order","Family","Genus","Species","Base")
#' @param comparison A vector of conditions to compare. Default: NULL, all unique conditions are compared (only for Two groups).
#' @param cutoff The log2 fold change cutoff for considering a taxon as differentially expressed.
#' @param control_name Character. The name of the control group for the comparison.
#' @param paired Logical. Should the samples be treated as paired? Default: False
#' @param subject Optional. The subject identifier for paired samples. Default: Null
#'
#' @return A data frame with the results of the differential expression analysis.
#' @export
#' @importFrom stats relevel
#'
#' @note
#' 1. Regulation is judged by cutoff of q-value(adjust p value).Detail see in \code{\link[DESeq2]{DESeq}}
#'
#' 2. For more than two groups in taxobj, the 'comparison' must be assigned.
#'
#' 3. The function requires the 'DESeq2', 'S4Vectors', and 'tibble' packages.
#'
#' @seealso \code{\link[DESeq2]{DESeqDataSetFromMatrix}}, \code{\link[DESeq2]{DESeq}}, \code{\link[S4Vectors]{DataFrame}}, \code{\link[tibble]{as_tibble}}
#' @author  Wang Ningqi <2434066068@qq.com>
#' @examples
#' \donttest{
#' if (requireNamespace("DESeq2", quietly = TRUE) &&
#'     requireNamespace("S4Vectors", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'
#'     ### Data preparation ###
#'     data("Two_group")
#'
#'     ### Deseq analysis ###
#'     deseq_results <- Deseq_analysis(
#'       taxobj = Two_group,
#'       taxlevel = "Genus",
#'       cutoff = 1,
#'       control_name = "Control"
#'     )
#'
#'     # Visualization of volcano plot ##
#'     volcano_plot <- volcano_plot(
#'       inputframe = deseq_results,
#'       cutoff = 1,
#'       aes_col = Two_group$configuration$treat_col
#'     )
#'     volcano_plot$FC_FDR
#'     volcano_plot$Mean_FC
#'
#'     # Visualization of Manhattan plot ##
#'     manhattan_object <- manhattan(
#'       inputframe = deseq_results,
#'       taxlevel = "Phylum",
#'       control_name = "Control",
#'       mode = "most",
#'       top_n = 10,
#'       rmprefix = "p__"
#'     )
#'     manhattan_object$manhattan  # Tradition manhattan plot
#'     manhattan_object$manhattan_circle  # Circular manhattan plot
#'
#'     # For object with more than two groups
#'     ### Data preparation ###
#'     data("Three_group")
#'
#'     # Specific comparison
#'     deseq_results_BFCF <- Deseq_analysis(
#'       taxobj = Three_group,
#'       taxlevel = "Genus",
#'       comparison = c("BF", "CF"),
#'       cutoff = 1,
#'       control_name = "CF"
#'     )
#'     volcano_plot <- volcano_plot(
#'       inputframe = deseq_results_BFCF,
#'       cutoff = 1,
#'       aes_col = Three_group$configuration$treat_col
#'     )
#'     volcano_plot$FC_FDR
#'   } else {
#'     message(
#'       "The 'DESeq2', 'S4Vectors', and/or 'tibble' package(s) are not installed. ",
#'       "Please install them to use all features of Deseq_analysis."
#'     )
#'   }
#'}
Deseq_analysis<-function(taxobj,taxlevel,comparison=NULL,cutoff,control_name,paired=FALSE,subject=NULL){
  if (!requireNamespace("DESeq2", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("The 'DESeq2', 'S4Vectors', and 'tibble' package(s) are required but not installed. Please install them to use this function.")
  }
  if(is.null(taxobj$configuration)){
    stop("taxonomic summary object not configured yet, call '?object_config' for configuration!")
    return(NULL)
  }
  if(is.null(eval(parse(text=paste0("taxobj","$",taxlevel))))){
    warning("Illegal 'taxlevel'!")
    return(NULL)
  }
  condition=eval(parse(text=paste0("taxobj$Groupfile")))
  condition=condition[eval(parse(text=paste0("taxobj$configuration$treat_location"))) ]
  condition=condition[,1]
  if(is.null(comparison)){
    comparison= unique(condition)
  }else{
    if((which(condition %in% comparison) %>% length())==0){
      stop("comparision does not match, Please check 'comparision'")
    }
    taxobj=sub_tax_summary(taxobj=taxobj,specificnum = which(condition %in% comparison))
    condition=eval(parse(text=paste0("taxobj$Groupfile")))
    condition=condition[eval(parse(text=paste0("taxobj$configuration$treat_location"))) ]
    condition=condition[,1]
  }
  taxonomy=eval(parse(text=paste0("taxobj","$",taxlevel,"_taxonomy")))
  input0=eval(parse(text=paste0("taxobj","$",taxlevel)))
  inputframe=data.frame(input0[,-1],row.names =taxonomy[,1]) %>% round(0)

  group_length=length(comparison)
  if(group_length!=2){
    stop("Comparsion group not assigned!",call. = FALSE)
    return(NULL)
  }
  if(paired==TRUE){
    ofdds = DESeq2::DESeqDataSetFromMatrix(countData = inputframe, S4Vectors::DataFrame(condition,subject),
                                   ~subject+condition)
  }else{
    ofdds=DESeq2::DESeqDataSetFromMatrix(countData=inputframe, S4Vectors::DataFrame(condition), ~ condition)
  }
  ofdds$condition<- relevel(ofdds$condition, ref = control_name)
  resultdds=as.data.frame(tibble::as_tibble(DESeq2::results(DESeq2::DESeq(ofdds))));rownames(resultdds)=rownames(inputframe)
  resultdds$tag[resultdds$padj<0.05 & resultdds$log2FoldChange>cutoff]= comparison[comparison!=control_name]
  resultdds$tag[resultdds$padj<0.05 & resultdds$log2FoldChange<(-cutoff)]= control_name
  resultdds$tag[which(is.na(resultdds$tag))]="None"
  resultdds<-resultdds[which(resultdds$baseMean>0),]
  resultdds$padj[which(is.na(resultdds$padj))]<- 0.999
  resultdds$ID=rownames(resultdds)
  colnames(resultdds)[ncol(resultdds)]=colnames(taxonomy)[1]
  resultdds_anno=left_join(resultdds,taxonomy) %>% suppressMessages()
  return(resultdds_anno)
}
