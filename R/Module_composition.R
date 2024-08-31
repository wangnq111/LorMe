#' Pie chart for network module composition
#' @description This function analyzes the composition of modules within a network object,
#' providing a visual and data summary based on taxonomic levels.
#'
#' @param network_obj Network analysis results generated from \code{\link{network_analysis}}
#' @param No.module Numeric or numeric vector of No.module
#' @param taxlevel  taxonomy levels used for visualization.Must be one of c("Domain","Phylum","Class","Order","Family","Genus","Species","Base").Default:NULL.
#' @param mode The mode for selecting which taxa to plot: "all" for all taxa, "most" for the top N taxa, and "select" for specific taxa selection
#' @param top_n The number of top taxa to plot when mode is set to "most"
#' @param palette Character. Palette for visualization,default:"Set1".See optional palette in same as 'RColorBrewer'. And "Plan1" to "Plan10" were also optional,see in \code{\link{color_scheme}}
#' @param select_tax A vector of taxa to be selected for plotting when mode is "select".
#' @param rmprefix A string prefix to be removed from the taxonomic annotation
#'
#' @return The function returns a list containing pie chart of specific module,corresponding source data and color assignments
#' @export
#' @importFrom stringr str_replace
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar labs theme
#' @importFrom stats aggregate
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom grDevices colorRampPalette
#
#' @examples
#'  #Data loading
#'  data("Two_group")
#'
#'  # Network analysis
#'  network_Two_group <- network_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus",
#'     reads = TRUE,
#'     n = 8,
#'     threshold = 0.7
#'   )
#'
#'   # Show all modules
#'   module_results <- Module_composition(
#'     network_obj = network_Two_group,
#'     No.module = c(2, 5),
#'     taxlevel = "Phylum"
#'   )
#'   print(module_results$Module5$Pie)
#'   print(module_results$Module2$Pie)  # View pie chart
#'   head(module_results$Module2$source_data_Module2)  # View source data for pie chart
#'   print(module_results$aes_color)  # Check aesthetic color
#'
#'   # Show taxa with top five frequency
#'   module_results <- Module_composition(
#'     network_obj = network_Two_group,
#'     No.module = c(2, 5),
#'     taxlevel = "Phylum",
#'     mode = "most",
#'     top_n = 5
#'   )
#'   print(module_results$Module2$Pie_plot_Module2)
#'
#'   # Show specific taxa
#'   community <- community_plot(
#'     taxobj = Two_group,
#'     taxlevel = "Phylum",
#'     n = 5,
#'     palette = "Paired"
#'   )  # Get top 5 dominant phyla
#'   top5_phyla <- names(community$filled_color)
#'
#'   module_results <- Module_composition(
#'     network_obj = network_Two_group,
#'     No.module = c(2, 5),
#'     taxlevel = "Phylum",
#'     mode = "select",
#'     palette = community$filled_color,
#'     select_tax = top5_phyla
#'   )
#'   print(module_results$Module2$Pie_plot_Module2)
#'
#'   # Specific taxa with no prefix 'p__'
#'   module_results <- Module_composition(
#'     network_obj = network_Two_group,
#'     No.module = 2,
#'     taxlevel = "Phylum",
#'     mode = "select",
#'     select_tax = c("Proteobacteria", "Actinobacteria")
#'   )
#'   print(module_results$Module2$Pie_plot_Module2)
#'
#'   # Remove 'p__' prefix
#'   module_results <- Module_composition(
#'     network_obj = network_Two_group,
#'     No.module = 2,
#'     taxlevel = "Phylum",
#'     mode = "most",
#'     top_n = 5,
#'     palette = "Set2",
#'     rmprefix = "p__"
#'   )
#'   print(module_results$Module2$Pie_plot_Module2)
Module_composition = function(
    network_obj,
    No.module,
    taxlevel = "Phylum",
    mode = "all",
    top_n=NULL,
    palette="Set1",
    select_tax=NULL,
    rmprefix=NULL
) {
  if(!taxlevel %in%c("Domain","Phylum","Class","Order","Family","Genus","Species")){
    warning("Invalid input in 'taxlevel',must select among  'Domain''Phylum''Class''Order''Family''Genus''Species'")
  }else if(!taxlevel %in% (network_obj$Nodes_info %>% colnames())){
    warning("Invalid taxlevel for network object")
  }else{
    output_list=list()
    input_table = as.data.frame(network_obj$Nodes_info)
    if(is.null(rmprefix)==FALSE){
      input_table[,taxlevel]=str_replace(input_table[,taxlevel],rmprefix,"")
    }
    for(i in No.module){
      tax_freq = input_table[, taxlevel]
      tax_freq=tax_freq[input_table$No.module == i] %>%
        table() %>% sort(decreasing = TRUE) %>% as.data.frame()
      colnames(tax_freq)[1] = "Tax"
      tax_freq["Tax"]=as.character(tax_freq$Tax)
      if (mode == "all") {
        tax_freq_summary = tax_freq
      }else if(mode=="most"){
        if(is.null(top_n)){
          warning("'top_n' not assigned!")
          return(NULL)}
        kept_tax=table(input_table[,taxlevel]) %>% sort() %>% names() %>% tail(top_n)
        tax_freq$Tax[!(tax_freq["Tax"]%>% as.matrix() %>% as.character()) %in% kept_tax] ="Others"
        tax_freq_summary=aggregate(tax_freq$Freq,by=list(Tax=tax_freq$Tax),FUN=sum)
        colnames(tax_freq_summary)[2]="Freq"
        tax_freq_summary=tax_freq_summary[order(tax_freq_summary["Freq"]%>% as.matrix() %>%as.numeric(),decreasing = TRUE),]
      }else if(mode=="select"){
        if(is.null(select_tax)){
          warning("'select_tax' not assigned!")
          return(NULL)}
        tax_freq$Tax[!(tax_freq["Tax"]%>% as.matrix() %>% as.character())%in% select_tax] ="Others"
        tax_freq_summary=aggregate(tax_freq$Freq,by=list(Tax=tax_freq$Tax),FUN=sum)
        colnames(tax_freq_summary)[2]="Freq"
        tax_freq_summary=tax_freq_summary[order(tax_freq_summary["Freq"]%>% as.matrix() %>%as.numeric(),decreasing = TRUE),]
      }else{
        warning("Invalid 'mode', please choose among 'all','most' and 'select'")
        return(NULL)
      }
      tax_freq_summary$Tax=factor(tax_freq_summary$Tax,levels=c(tax_freq_summary$Tax[which(tax_freq_summary$Tax!="Others")],"Others"))
      if(mode=="all"){
        tax_length=unique(input_table$Phylum) %>% length()
      }else if(mode=="most"){
        tax_length=top_n+1
      }else if(mode=="select"){
        tax_length=length(select_tax)+1
      }
      if(length(palette)>1){
        tax_col=palette
      }else if(palette %in% paste0("Plan",1:10)){
        tax_col=color_scheme(Plan = palette,
                             expand = tax_length,show=FALSE)
      }else{
        color_n<-brewer.pal.info[palette,]$maxcolors
        getPalette <-colorRampPalette(brewer.pal(color_n, palette))
        tax_col=getPalette(tax_length)
        if(mode=="all"){
          names(tax_col)=input_table[,taxlevel] %>% unique()
        }else if(mode=="select"){
          names(tax_col)=c(select_tax,"Others")
        }else if(mode=="most"){
          names(tax_col)=c(kept_tax,"Others")
        }
      }
      p=ggplot(tax_freq_summary, aes_string(x="0", y='Freq', fill='Tax'))+
        geom_bar(width = 1, stat = "identity",show.legend = TRUE)+
        coord_polar("y", start=0)+
        scale_fill_manual(values = tax_col)+
        labs(x="",y="",title = paste0("Module #",i),fill=taxlevel)+
        theme(panel.grid = element_blank(), panel.background = element_blank(),
              axis.text = element_blank(),axis.ticks=element_blank(),plot.title = element_text(hjust = 0.5,size=8))
      plot=p
      output=list(plot,tax_freq_summary)
      names(output)=c(paste0("Pie_plot_Module",i),paste0("source_data_Module",i))
      output_list=c(output_list,list(output))
      names(output_list)[length(output_list)] =paste0("Module",i)
    }
  }
  output_list=c(output_list,list(tax_col))
  names(output_list)[length(output_list)]="aes_color"
  return(output_list)
}
