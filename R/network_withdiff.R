#' Network Analysis with Differential Species
#'
#' Meta network analysis integrating differential taxon into a network analysis
#' @param network_obj Network analysis results generated from \code{\link{network_analysis}}
#' @param diff_frame Differential analysis results generated from \code{\link{indicator_analysis}} or \code{\link{Deseq_analysis}}.
#' @param aes_col A named vector of colors to be used to highlight differential taxon vertices
#' @param tag_threshold Numeric. A threshold for the minimum number of differential taxon to display.
#'
#' @return A list containing the configured igraph object, vertices coordinates, parameters, and tag statistics.
#' @export
#'
#' @importFrom dplyr left_join
#' @importFrom igraph V cluster_fast_greedy as.undirected membership with_fr
#' @importFrom tidyr spread
#' @importFrom graphics par legend
#'
#' @examples
#' \donttest{
#' {
#'   # Data preparation
#'   data("Two_group")
#'   set.seed(999)
#'   # Analysis
#'   network_results <- network_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus",
#'     n = 10,
#'     threshold = 0.8
#'   )
#'   indicator_results <- indicator_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus"
#'   )
#'   deseq_results <- Deseq_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Genus",
#'     cutoff = 1,
#'     control_name = "Control"
#'   )
#'
#'   # Visualize
#'   network_diff_obj <- network_withdiff(
#'     network_obj = network_results,
#'     diff_frame = indicator_results
#'   )
#'   # Check contained tags for each model
#'   print(network_diff_obj$tag_statistics$sum_of_tags)
#'   # Check contained different tags for each model
#'   print(network_diff_obj$tag_statistics$detailed_tags)
#'
#'   # Re-visualize
#'   network_visual_re(
#'     network_visual_obj = network_diff_obj,
#'     module_paint = TRUE,
#'     module_num = c(1, 4)
#'   )  # Show module with most Treatment indicators
#'
#'   my_module_palette <- color_scheme(
#'     c("#83BA9E", "#F49128"),
#'     5
#'   )
#'   network_visual_re(
#'     network_visual_obj = network_diff_obj,
#'     module_paint = TRUE,
#'     module_num = c(1, 4, 6, 3, 8),
#'     module_palette = my_module_palette
#'   )  # Show module with most Treatment indicators
#'
#'   # Available also for DESeq analysis results
#'   network_diff_obj <- network_withdiff(
#'     network_obj = network_results,
#'     diff_frame = deseq_results
#'   )
#'
#'   # Parameter adjustment
#'   network_diff_obj <- network_withdiff(
#'     network_obj = network_results,
#'     diff_frame = indicator_results,
#'     tag_threshold = 20
#'   )  # The 'tag_threshold' set too high
#'
#'   network_diff_obj <- network_withdiff(
#'     network_obj = network_results,
#'     diff_frame = indicator_results,
#'     tag_threshold = 10
#'   )  # Set lower
#'   # Check contained tags for each model
#'   print(network_diff_obj$tag_statistics$sum_of_tags)
#'   # Check contained different tags for each model
#'   print(network_diff_obj$tag_statistics$detailed_tags)
#'
#'   network_diff_obj <- network_withdiff(
#'     network_obj = network_results,
#'     diff_frame = indicator_results,
#'     tag_threshold = 1
#'   )  # Set too low
#'
#'   # Another example
#'   data("Three_group")
#'   network_results <- network_analysis(
#'     taxobj = Three_group,
#'     taxlevel = "Genus",
#'     n = 15,
#'     threshold = 0.9
#'   )
#'   indicator_results <- indicator_analysis(
#'     taxobj = Three_group,
#'     taxlevel = "Genus"
#'   )
#'
#'   tag_color <- c(
#'     "CF" = "#F8766D",
#'     "CF_OF" = "#FFFF00",
#'     "OF" = "#00BA38",
#'     "OF_BF" = "#800080",
#'     "BF" = "#619CFF",
#'     "CF_BF" = "#00FFFF"
#'   )
#'   network_diff_obj <- network_withdiff(
#'     network_obj = network_results,
#'     diff_frame = indicator_results,
#'     aes_col = tag_color,
#'     tag_threshold = 10
#'   )
#'
#'   # Re-visualize
#'   print(network_diff_obj$tag_statistics$detailed_tags)
#'   network_visual_re(
#'     network_visual_obj = network_diff_obj,
#'     module_paint = TRUE,
#'     module_num = c(8, 10, 11)
#'   )  # Show module with most BF indicators
#'   network_visual_re(
#'     network_visual_obj = network_diff_obj,
#'     module_paint = TRUE,
#'     module_num = c(1, 6, 8)
#'   )  # Show module with most BF and OF_BF indicators
#' }
#' }
network_withdiff=function(network_obj,diff_frame,aes_col=NULL,tag_threshold=5){
  adjfile=network_obj$Adjacency_column_table
  verticefile=network_obj$Nodes_info
  verticefile$ID=verticefile$nodes_id
  if(network_obj$config$taxlevel!="Base"){
    matchID=paste0(network_obj$config$taxlevel,"ID")
    colnames(verticefile)[ncol(verticefile)]=matchID
  }else{
    matchID="ID"
  }
  verticefile=left_join(verticefile,diff_frame[,c(matchID,"tag")]) %>% suppressMessages()
  t_net_its <- graph_from_data_frame(adjfile,directed=FALSE,vertices = verticefile)
  #color
  if(length(aes_col)==1){
    warning("Only one element found in 'aes_col'!Automaticlly assigned yet")
    aes_col=NULL
  }
  if(is.null(aes_col)){
    aes_col=color_scheme("Plan3",length(unique(diff_frame$tag))-1,show = FALSE) %>% suppressMessages()
    names(aes_col)=unique(diff_frame$tag)[unique(diff_frame$tag)!="None"]
  }else{
    if(is.null(names(aes_col))){
      aes_col=color_scheme(aes_col,length(unique(diff_frame$tag))-1,show = FALSE)%>% suppressMessages()
      names(aes_col)=unique(diff_frame$tag)[unique(diff_frame$tag)!="None"]
    }else{
      aes_collength=names(aes_col) %in% unique(diff_frame$tag)%>% length()
      taglength=unique(diff_frame$tag)
      taglength=taglength[which(taglength!="None")]%>% length()
      if(aes_collength!=taglength){
        warning("Names of 'aes_col' do not match tags! Please check!")
        return(NULL)
      }else(aes_col=aes_col)
    }
  }
  V(t_net_its)$color <- V(t_net_its)$tag
  V(t_net_its)$color[!V(t_net_its)$color =="None"]=aes_col[V(t_net_its)$color[!V(t_net_its)$color =="None"]]
  V(t_net_its)$color[V(t_net_its)$color=="None"]="gray"
  #size
  V(t_net_its)$size <- V(t_net_its)$tag
  V(t_net_its)$size[V(t_net_its)$size =="None"] <- 3
  V(t_net_its)$size[V(t_net_its)$size %in% names(aes_col)] <- 6

  vertice_size <- as.numeric(V(t_net_its)$size)

  cfg_t=cluster_fast_greedy(as.undirected(t_net_its))
  t_modules=sort(table(membership(cfg_t)),decreasing = TRUE)
  t_mods_list_cs <- list()
  for (i in names(t_modules)){
    x1 <- names(membership(cfg_t)[membership(cfg_t)==i])
    x2 <- x1[x1 %in% verticefile$nodes_id]
    t_mods_list_cs[[i]] <- as.numeric(V(t_net_its)[x2])
  }

  all_stat=table(verticefile$No.module[verticefile$tag!="None"]) %>% as.data.frame()
  colnames(all_stat)=c("No.module","sum_tag_number")
  all_stat=all_stat[order(all_stat$sum_tag_number,decreasing = TRUE),]
  show_module=all_stat$No.module[all_stat$sum_tag_number>=tag_threshold] %>% as.character()
  if(length(show_module)==0){
    warning("None of modules contained more than ",tag_threshold," tag vertices, please decrease 'tag_threshold'")
    return(NULL)
  }
  detailed_stat=table(verticefile[verticefile$tag!="None",c("No.module","tag")])  %>%
    as.data.frame() %>% spread("tag","Freq")
  coords_t_its <- layout_(t_net_its,with_fr(niter=9999, grid="auto"))
  t_cols <- color_scheme(c("aquamarine3","antiquewhite2","goldenrod2"),length(show_module),show = FALSE)
  show_list=t_mods_list_cs[names(t_mods_list_cs) %in% show_module]
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1,1), mar=c(2,2,2,4),pty = "m")
  plot(t_net_its,vertex.label=NA,vertex.size=vertice_size,layout=coords_t_its,vertex.shape="circle",
       mark.groups=show_list,mark.col=t_cols, mark.border="gray")
  legend_text=paste0("Module #",show_module)
  #plot.new()
  legend("topright", c(names(aes_col),"None"), col= c(aes_col,"gray"),
         text.col = c(aes_col,"gray"),
         #lty = c(2, -1, 1),
         pch = c(16, 16, 16), trace=TRUE,border="black",bty="n")
  legend("bottomright",legend=legend_text,col=t_cols, bty="n",fill=t_cols,border="gray",xjust = 1)
  parameters=c(list(vertice_size),list(t_mods_list_cs),list(aes_col))
  names(parameters)=c("vertice_size","module_list","aes_col")
  tag_freq=c(list(all_stat),list(detailed_stat))
  names(tag_freq)=c("sum_of_tags","detailed_tags")
  outlist=c(list(t_net_its),list(coords_t_its),list(parameters),list(tag_freq))
  names(outlist)=c("configured_igraph_object","vertices_coordinates","parameters","tag_statistics")
  return(outlist)
}
