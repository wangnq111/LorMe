#' Manhattan Plot Generator
#' @description Generate Manhattan Plot  base on Deseq_analysis or indicator_analysis results
#'
#' @param inputframe A data frame generated from \code{\link{Deseq_analysis}} or \code{\link{indicator_analysis}}
#' @param taxlevel Taxonomy levels used for visualization.Must be one of c("Domain","Phylum","Class","Order","Family","Genus","Species","Base").Default:NULL.
#' @param control_name Character. The name of the control group for the comparison.
#' @param mode The mode for selecting which taxa to plot: "all" for all taxa, "most" for the top N taxa, and "select" for specific taxa selection
#' @param top_n The number of top taxa to plot when mode is set to "most"
#' @param palette Character. Palette for visualization,default:"Set1".Optional palette same as 'RColorBrewer'. "Plan1" to "Plan10" were also optional,see in \code{\link{color_scheme}}
#' @param select_tax A vector of taxa to be selected for plotting when mode is "select".
#' @param rmprefix A string prefix to be removed from the taxonomic annotation
#'
#'
#' 1) \code{\link{ggplot2::scale_colour_brewer}}
#'
#' 2) \code{\link{color_scheme}}. Can be "Plan1" to "Plan10"
#'
#' 3) Character string with specific color.
#' @param select_tax A vector of taxa to be selected for plotting when mode is "select".
#' @param rmprefix A string prefix to be removed from the taxonomic annotation.Default:NULL.
#'
#' @return a list containing the Manhattan plot, circular Manhattan plot, source data, and color assignments
#' @export
#'
#' @importFrom ggplot2 position_jitter scale_shape_manual geom_text scale_color_manual
#' @importFrom utils head
#' @importFrom stringr str_replace
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @examples
#' \donttest{
#' {
#'   # Data preparation
#'   data("Two_group")
#'
#'   # DESeq analysis
#'   deseq_results <- Deseq_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Base",
#'     cutoff = 1,
#'     control_name = "Control"
#'   )
#'
#'   # Indicator analysis
#'   indicator_results <- indicator_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus"
#'   )
#'
#'   # Show all with Manhattan plot
#'     manhattan_object <- manhattan(
#'       inputframe = deseq_results,
#'       taxlevel = "Phylum",
#'       control_name = "Control"
#'     )
#'     print(manhattan_object$manhattan)  # Tradition Manhattan plot
#'     print(manhattan_object$manhattan_circle)  # Circular Manhattan plot
#'     print(manhattan_object$sourcedata)  # Source data for plot
#'     print(manhattan_object$aes_color)  # Aesthetic color for plot
#'
#'   # Top 8 Phyla with most taxon
#'     manhattan_object <- manhattan(
#'       inputframe = indicator_results,
#'       taxlevel = "Phylum",
#'       control_name = "Control",
#'       mode = "most",
#'       top_n = 8,
#'       palette = "Set1"
#'     )
#'     print(manhattan_object$manhattan)
#'
#'   # Specific phyla
#'   # Top nine dominant phyla
#'     community <- community_plot(
#'       taxobj = Two_group,
#'       taxlevel = "Phylum",
#'       n = 9,
#'       palette = "Paired",
#'       rmprefix = "p__"
#'     )
#'
#'     manhattan_object <- manhattan(
#'       inputframe = indicator_results,
#'       taxlevel = "Phylum",
#'       control_name = "Control",
#'       mode = "select",
#'       palette = community$filled_color,
#'       select_tax = names(community$filled_color),
#'       rmprefix = "p__"
#'     )
#'     print(manhattan_object$manhattan)
#'     print(manhattan_object$manhattan_circle)
#' }
#'}
manhattan=function(inputframe,taxlevel="Phylum",control_name,mode="all",top_n=NULL,palette="Set1",select_tax=NULL,rmprefix=NULL){
  if(!taxlevel %in% colnames(inputframe)){
    warning("'taxlevel' not included in inputframe!")
    return(NULL)
  }
  if(!control_name %in% inputframe$tag){
    warning("'control_name' not included in treatment!")
    return(NULL)
  }

  inputframe=inputframe[!is.na(inputframe[,'padj']),]
  tax_stat=table(inputframe[,taxlevel]) %>% sort()
  rm_tax=tax_stat[which(tax_stat<=2)] %>% names()
  inputframe=inputframe[!inputframe[,taxlevel] %in% rm_tax,]
  if(is.null(rmprefix)==FALSE){
    inputframe[,taxlevel]=str_replace(inputframe[,taxlevel],rmprefix,"")
  }
  tax_order=table(inputframe[,taxlevel]) %>% sort() %>% names() %>% rev()
  if(mode=="all"){
    inputframe[,taxlevel]=factor(inputframe[,taxlevel],levels =tax_order)
    inputframe_select=inputframe
  }else if(mode=="most"){
    if(is.null(top_n)){
      warning("'top_n' not assigned!")
      return(NULL)}
    select_tax_order=tax_order %>% head(top_n)
    inputframe_select=inputframe[which(inputframe[,taxlevel] %in% select_tax_order),]
    inputframe_select[,taxlevel]=factor(inputframe_select[,taxlevel],levels =select_tax_order)
  }else if(mode=="select"){
    if(is.null(select_tax)){
      warning("'select_tax' not assigned!")
      return(NULL)}
    inputframe_select=inputframe[inputframe[,taxlevel] %in% select_tax,]
    inputframe_select[,taxlevel]=factor(inputframe_select[,taxlevel],levels = tax_order[tax_order %in% select_tax] )
  }else{
    warning("Invalid 'mode', please choose among 'all','most' and 'select'")
    return(NULL)
  }
  aes_shape=c(25,20,17)
  treatment_unique=unique(inputframe$tag)
  treatment_name= treatment_unique[!treatment_unique%in% c("None",control_name)]
  names(aes_shape)=c(control_name,"None",treatment_name)

  #visualization
  p_width=ggplot(inputframe_select, aes(x=inputframe_select[,taxlevel], y=-log10(inputframe_select[,'padj']))) +
    geom_hline(yintercept=-log10(0.05), linetype=2, color="lightgrey") +
    geom_point(aes(color=inputframe_select[,taxlevel],shape=inputframe_select[,'tag']),alpha=.8,position=position_jitter(0.5)) +
    scale_shape_manual(values=aes_shape)+
    #scale_size_manual(values=c("IIII"=2,"III"=1.5,"II"=1,"I"=0.5))+
    labs(x=NULL, y="-log10(q value)",shape="")+
    scale_y_continuous(expand = c(0,0))+
    theme_zg()+
    guides(color="none")+
    theme(legend.position="top", #legend
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.1, vjust = 0.1))
  #
  max_y=min(log10(inputframe_select$padj))
  circile_text=data.frame(x=0,y=-seq(1,-max_y,2),label=seq(1,-max_y,2))
  p_circle=ggplot(inputframe_select, aes(x=inputframe_select[,taxlevel], y=log10(inputframe_select[,'padj']))) +  #
    geom_hline(yintercept=0, linetype=2, color="lightgrey") +
    geom_hline(yintercept=log10(0.05), linetype=2, color="lightgrey") +
    geom_hline(yintercept=log10(0.001), linetype=2, color="lightgrey") +
    geom_hline(yintercept=log10(1e-5), linetype=2, color="lightgrey") +
    geom_point(aes(color=inputframe_select[,taxlevel],shape=inputframe_select[,'tag']),alpha=.8,position=position_jitter(0.5)) +
    scale_shape_manual(values=aes_shape)+
    #scale_size_manual(values=c("IIII"=2,"III"=1.5,"II"=1,"I"=0.5))+
    scale_y_continuous(limits = c(1.5*max_y,0.5))+
    labs(x=NULL, y="",shape="")+
    theme_zg()+
    guides(color="none")+coord_polar(theta = "x")+
    theme(#legend.position="top",
      panel.background = element_blank(),
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank(),
      #panel.grid = element_blank(),
      axis.text.x = element_text())
  p_circle=p_circle+geom_text(data=circile_text,aes(x=circile_text[,'x'],y=circile_text[,'y'],label=circile_text[,'label']),size=4)
  #color assign
  tax_length=inputframe_select[,taxlevel] %>% unique() %>% length()
  if(length(palette)>1){
    tax_col=palette
  }else if(palette %in% paste0("Plan",1:10)){
    tax_col=color_scheme(Plan = palette,
                         expand = tax_length,show=FALSE)
  }else{
    color_n<-brewer.pal.info[palette,]$maxcolors
    getPalette <-colorRampPalette(brewer.pal(color_n, palette))
    tax_col=getPalette(tax_length)
    names(tax_col)=inputframe_select[,taxlevel] %>% unique()
  }
  p_width=p_width+
    scale_color_manual(values=tax_col)
  p_circle=p_circle+
    scale_color_manual(values=tax_col)
  output=c(list(p_width),list(p_circle),list(inputframe_select),list(tax_col))
  names(output)=c("manhattan","manhattan_circle","sourcedata","aes_color")
  return(output)
}
