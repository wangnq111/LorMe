#' Network Visualization
#'
#' Visualizes a network based on a network object from \code{\link{network_analysis}}.
#'
#' @param network_obj A network analysis results object generated from \code{\link{network_analysis}}.
#' @param mode The visualization mode, optionally "major_module" or "major_tax".
#' @param major_num The number of major modules to display in the network.
#' @param taxlevel Taxonomy levels used for visualization when mode is "major_tax".
#' @param select_tax A vector of taxa to be selected for displaying in "major_tax" mode.
#' @param palette Character. Palette for visualization.
#' @param vertex.size Numeric. The size of the vertices.
#'
#' @return A list containing the configured igraph object and the coordinates of the vertices, with network visualization displayed in the plots panel.
#' @export
#'
#' @importFrom igraph graph_from_data_frame layout_ V with_fr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics par legend
#'
#' @examples
#' {
#'   # Data preparation
#'   data("Two_group")
#'   set.seed(999)
#'   # Analysis
#'   network_results <- network_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Species",
#'     n = 10,
#'     threshold = 0.6
#'   )
#'
#'   # Default mode
#'   network_visual_obj <- network_visual(network_obj = network_results)
#'
#'   # View again
#'   network_visual_re(network_visual_obj)
#'
#'   # More modules
#'   network_visual_obj <- network_visual(
#'     network_obj = network_results,
#'     major_num = 10
#'   )
#'
#'   # Specific tax
#'   # Generate top 5 phyla for displaying
#'   community <- community_plot(
#'     taxobj = Two_group,
#'     taxlevel = "Phylum",
#'     n = 5,
#'     palette = "Paired"
#'   )
#'   display_phyla <- names(community$filled_color)
#'
#'   network_visual_obj <- network_visual(
#'     network_obj = network_results,
#'     mode = "major_tax",
#'     taxlevel = "Phylum",
#'     select_tax = display_phyla,
#'     palette = community$filled_color
#'   )
#'
#'   # Another sample for specific tax
#'   network_visual_obj <- network_visual(
#'     network_obj = network_results,
#'     mode = "major_tax",
#'     taxlevel = "Phylum",
#'     select_tax = "p__Proteobacteria"
#'   )
#' }
network_visual=function(network_obj,mode="major_module",major_num=5,taxlevel=NULL,select_tax=NULL,palette="Set1",vertex.size=6){
  adjfile=network_obj$Adjacency_column_table
  verticefile=network_obj$Nodes_info

  t_net_its <- graph_from_data_frame(adjfile,directed=FALSE,vertices = verticefile)
  #color
  if(mode=="major_module"){
    module_stat=V(t_net_its)$No.module %>% table() %>% sort(decreasing=TRUE)
    if(major_num>length(module_stat)){
      major_num=length(module_stat)
      warning("Module number less than 'major_num',all modules were colored")
    }
    if(length(palette)>1){
      tax_col=color_scheme(Plan = palette,
                           expand = major_num,show=FALSE)
    }else if(palette %in% paste0("Plan",1:10)){
      tax_col=color_scheme(Plan = palette,
                           expand = major_num,show=FALSE)
    }else{
      color_n<-brewer.pal.info[palette,]$maxcolors
      getPalette <-colorRampPalette(brewer.pal(color_n, palette))
      tax_col=getPalette(major_num)
    }
    select_module=names(module_stat)[1:major_num]
    names(tax_col)=select_module
    V(t_net_its)$color <- V(t_net_its)$No.module
    V(t_net_its)$color[!V(t_net_its)$color %in%select_module]="gray"
    V(t_net_its)$color[V(t_net_its)$color %in%select_module]=tax_col[V(t_net_its)$color[V(t_net_its)$color %in%select_module]]
  }else if(mode=="major_tax"){
    if(!taxlevel %in% colnames(verticefile)){
      warning("'taxlevel' not included in network_obj!")
      return(NULL)
    }
    if(is.null(select_tax)){
      warning("'select_tax' not set! Invalid input!")
      return(NULL)
    }
    if(length(palette)>1){
      if(is.null(names(palette))){
        tax_col=color_scheme(Plan = palette,
                             expand = length(select_tax),show=FALSE)
        names(tax_col)=select_tax
      }else{
        tax_col=palette
      }
    }else if(palette %in% paste0("Plan",1:10)){
      tax_col=color_scheme(Plan = palette,
                           expand = length(select_tax),show=FALSE)
      names(tax_col)=select_tax
    }else{
      color_n<-brewer.pal.info[palette,]$maxcolors
      getPalette <-colorRampPalette(brewer.pal(color_n, palette))
      tax_col=getPalette(length(select_tax))
      names(tax_col)=select_tax
    }

    V(t_net_its)$color <- eval(parse(text=paste0("V(t_net_its)$",taxlevel)))
    V(t_net_its)$color[!V(t_net_its)$color %in% select_tax]="gray"
    V(t_net_its)$color[V(t_net_its)$color %in% select_tax]=tax_col[V(t_net_its)$color[V(t_net_its)$color %in% names(tax_col)]]
  }else{
    warning("Invalid 'mode'! Please select between 'major_module' and 'major_tax'!")
    return(NULL)
  }
  coords_t_its <- layout_(t_net_its,with_fr(niter=9999, grid="auto"))
  plot(t_net_its,vertex.label=NA,vertex.size=vertex.size,layout=coords_t_its,vertex.shape="circle")
  outlist=c(list(t_net_its),list(coords_t_its))
  names(outlist)=c("configured_igraph_object","vertices_coordinates")
  return(outlist)
}
