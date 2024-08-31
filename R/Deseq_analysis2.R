###version1.1.0 ###
###author Wangningqi####
#'  Deseq analysis
#' @description Deseq analysis
#' @param inputframe Otu/gene/taxa table with all integer numeric variables.Rownames must be  Otu/gene/taxa names,colnames must be sample names with control in front and treatment behind. Reads table is recommended.
#' @param condition A character string which indicates group of samples
#' @param cutoff threshold of log2(Foldchange).Detail see in \code{\link[DESeq2]{DESeq}}
#' @param control_name A character indicating the control group name
#' @param paired Logical to determine if paired comparision would be used. TRUE or FALSE.
#' @param subject A character string which indicates paired design of samples
#'
#' @return Statistics dataframe of all otu/gene/taxa
#' @export
#' @note
#' 1. Inputframe must be all integer numeric variables without NA/NAN/inf! In case your data is not an integer one,a practical method is to multiply them in equal proportion(eg. x 1e6) then round them into integer
#'
#' 2. Regulation is judged by cutoff of qvalue(adjust p value).Detail see in \code{\link[DESeq2]{DESeq}}
#'
#' 3. Set cutoff as 1 is recommened.In case of too few taxa(eg. Phylum level deseq),cutoff can be set to 0.
#'
#' 4. if control_name is not given, the control group will be set according to ASCII
#'
#' 5. The function requires the 'DESeq2', 'S4Vectors', and 'tibble' packages.
#' @author  Wang Ningqi <2434066068@qq.com>
#' @importFrom stats relevel
#' @seealso \code{\link[DESeq2]{DESeqDataSetFromMatrix}}, \code{\link[DESeq2]{DESeq}}, \code{\link[S4Vectors]{DataFrame}}, \code{\link[tibble]{as_tibble}}
#' @examples
#' \donttest{
#' {
#'   ### Data preparation ###
#'   data(testotu)
#'   rownames(testotu) <- testotu[, 1]
#'   inputotu <- testotu[, -c(1, ncol(testotu))]
#'   head(inputotu)
#'   group <- c(rep("a", 10), rep("b", 10))
#'
#'   ### DESeq analysis ###
#'   if (requireNamespace("DESeq2", quietly = TRUE) &&
#'     requireNamespace("S4Vectors", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'     Deseqresult <- Deseq_analysis2(
#'       inputframe = inputotu,
#'       condition = group,
#'       cutoff = 1,
#'       control_name = "b"
#'     )
#'
#'     ### Paired DESeq analysis ###
#'     subject <- factor(c(1:10, 1:10))
#'     Deseqresult <- Deseq_analysis2(
#'       inputframe = inputotu,
#'       condition = group,
#'       cutoff = 1,
#'       control_name = "b",
#'       paired = TRUE,
#'       subject = subject
#'     )
#'   }
#'}
#'}
Deseq_analysis2<-function(inputframe,condition,cutoff,control_name,paired,subject){
  if (!requireNamespace("DESeq2", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("The 'DESeq2', 'S4Vectors', and 'tibble' package(s) are required but not installed. Please install them to use this function.")
  }
  if(missing(paired)){paired=FALSE}
  if(missing(subject)){subject=NULL}
  if(paired==TRUE){
    ofdds = DESeq2::DESeqDataSetFromMatrix(countData = inputframe, S4Vectors::DataFrame(condition,subject),
                                    ~subject+condition)
  }else{
    ofdds=DESeq2::DESeqDataSetFromMatrix(countData=inputframe, S4Vectors::DataFrame(condition), ~ condition)}
  ofdds$condition<- relevel(ofdds$condition, ref = control_name)
  resultdds=as.data.frame(tibble::as_tibble(DESeq2::results(DESeq2::DESeq(ofdds))));rownames(resultdds)=rownames(inputframe)
  resultdds$threshold[resultdds$padj<0.05 & resultdds$log2FoldChange>cutoff]= "up"
  resultdds$threshold[resultdds$padj<0.05 & resultdds$log2FoldChange<(-cutoff)]= "down"
  resultdds$threshold[which(is.na(resultdds$threshold))]="none"
  resultdds<-resultdds[which(resultdds$baseMean>0),]
  resultdds$padj[which(is.na(resultdds$padj))]<- 0.999
  return(resultdds)}
