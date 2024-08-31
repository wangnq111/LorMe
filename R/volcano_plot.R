###version2.0.0 ###
###author Wangningqi####
#' Generate Volcano plot base on Deseq_analysis or indicator_analysis results
#'
#' @param inputframe A data frame containing the results based on \code{\link{Deseq_analysis}} or \code{\link{indicator_analysis}} (only two group indicators)
#' @param cutoff  A numeric value specifying the fold change cutoff,should be the same as in \code{\link{Deseq_analysis}}
#' @param aes_col A named vector of colors to be used in the plots
#'
#' @return A list of two ggplot objects, one for the fold change versus adjusted p-value plot and
#'         another for the mean abundance versus fold change or enrichment factor plot.
#' @export
#' @author  Wang Ningqi <2434066068@qq.com>
#' @import ggplot2
#' @examples
#' ###data prepration###
#'\donttest{
#' {
#'   # Load data
#'   data("Two_group")
#'
#'   # Define color based on treatment column
#'   mycolor <- Two_group$configuration$treat_col
#'
#'   ### DESeq analysis ###
#'   deseq_results <- Deseq_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus",
#'     cutoff = 1,
#'     control_name = "Control"
#'   )
#'
#'   ### Or indicator analysis ###
#'   indicator_results <- indicator_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus"
#'   )
#'
#'   # Create volcano plot for DESeq results
#'   volcano_plot <- volcano_plot(
#'     inputframe = deseq_results,
#'     cutoff = 1,
#'     aes_col = mycolor
#'   )
#'   print(volcano_plot$FC_FDR)  # Fold Change and FDR values
#'   print(volcano_plot$Mean_FC)  # Mean Fold Change values
#'
#'   # Create volcano plot for indicator results
#'   volcano_plot <- volcano_plot(
#'     inputframe = indicator_results,
#'     cutoff = 1,
#'     aes_col = mycolor
#'   )
#'   print(volcano_plot$FC_FDR)  # Fold Change and FDR values
#'   print(volcano_plot$Mean_FC)  # Mean Fold Change values
#' }
#' }
volcano_plot=function(inputframe,cutoff=NULL,aes_col=c("#FE5C5C","#75ABDE")){
  inputframe=inputframe[which(!is.na(inputframe$padj)),]
  tag=unique(inputframe$tag)[unique(inputframe$tag)!="None"]
  if(is.null(names(aes_col)[1])){
    names(aes_col)[1]=tag[1]
    names(aes_col)[2]=tag[2]
    warning("Color names not assigned,use automatic assignment")}
  if(length(aes_col)==2){
    aes_col=c(aes_col,"None"="gray")
  }else if (length(aes_col)==3){
    aes_col=aes_col
  }else {
    warning("Invalid 'aes_col',please check!")
  }

  if("log2FoldChange" %in% colnames(inputframe)){
    if(is.null(cutoff)){
      warning("'cutoff' not set!")
      return(NULL)
      }
    inputframe=inputframe[inputframe$padj<0.999,]
    FC_FDR=ggplot(inputframe,aes(x=inputframe[,'log2FoldChange'],y=-log10(inputframe[,'padj'])))+
      geom_point(aes(fill=inputframe[,'tag']),alpha=0.6,color="black",pch=21)+
      labs(x="Log2(FoldChange)",y="-Log10(q-value)",fill="")+
      geom_hline(yintercept =-log10(0.05),lty=2,color="grey")+
      geom_vline(xintercept = c(cutoff,-cutoff),lty=2,color="grey")+
      scale_fill_manual(values = aes_col)+theme_zg()
    Mean_FC=ggplot(inputframe,aes(x=inputframe[,'baseMean'],y=inputframe[,'log2FoldChange']))+
      geom_point(aes(fill=inputframe[,'tag']),alpha=0.6,color="black",pch=21)+
      labs(x="Basemean",y="Log2(FoldChange)",fill="")+
      geom_hline(yintercept =0,lty=2,color="grey")+
      scale_fill_manual(values = aes_col)+theme_zg()
  }else if("stat" %in% colnames(inputframe)){
    inputframe$stat[inputframe[,1]==1]=-inputframe$stat[inputframe[,1]==1]
    FC_FDR=ggplot(inputframe,aes(x=inputframe[,'stat'],y=-log10(inputframe[,'padj'])))+
      geom_point(aes(fill=inputframe[,'tag']),alpha=0.6,color="black",pch=21)+
      labs(x="Enrichement factor",y="-Log10(q-value)",fill="")+
      geom_hline(yintercept =-log10(0.05),lty=2,color="grey")+
      geom_vline(xintercept = c(0),lty=2,color="grey")+
      scale_fill_manual(values = aes_col)+theme_zg()
    Mean_FC=ggplot(inputframe,aes(x=inputframe[,'baseMean'],y=inputframe[,'stat']))+
      geom_point(aes(fill=inputframe[,'tag']),alpha=0.6,color="black",pch=21)+
      labs(x="Relative abundance",y="Enrichement factor",fill="")+
      geom_hline(yintercept =0,lty=2,color="grey")+
      scale_fill_manual(values = aes_col)+theme_zg()
  }
  outlist=c(list(FC_FDR),list(Mean_FC))
  names(outlist)=c("FC_FDR","Mean_FC")
  return(outlist)
}

