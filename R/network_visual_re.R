#' Re-visualize network plot from \code{\link{network_visual}} or \code{\link{network_withdiff}}
#' @param network_visual_obj Network object from \code{\link{network_visual}} or \code{\link{network_withdiff}}
#' @param module_paint Logical. If network module should be painted. Only work for network object from \code{\link{network_withdiff}}.
#' @param module_num Numeric indicating which module to be painted.
#' @param module_palette Character string with at least two elements. Palette for painting modules.
#' @param vertex.size Numeric. The size of the vertices, default:6. Only for network object from \code{\link{network_visual}}
#' @param vertex.shape Character. The shape of vertices, default: "circle"
#'
#' @return NULL but visualization in plot panel.
#' @export
network_visual_re=function(network_visual_obj,module_paint=FALSE,module_num=NULL,module_palette=c("aquamarine3","antiquewhite2","goldenrod2"),vertex.size=6,vertex.shape="circle"){
  if(module_paint==FALSE){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1,1), mar=c(2,2,2,4),pty = "m")
    plot(network_visual_obj$configured_igraph_object,
         vertex.label=NA,
         vertex.size=vertex.size,
         layout=network_visual_obj$vertices_coordinates,
         vertex.shape=vertex.shape)
  }else{
    if(is.null(network_visual_obj$parameters)){
      warning("Invalid network_visual_obj for painting module")
      return(NULL)
    }
    t_mods_list_cs=network_visual_obj$parameters$module_list
    if(length(which(module_num %in% names(t_mods_list_cs)==TRUE))!=length(module_num)){
      warning("Invalid module_num! Please check if 'module_num' contained in 'network_visual_obj'")
      return(NULL)
    }
    show_list=t_mods_list_cs[names(t_mods_list_cs) %in% module_num]
    t_cols <- color_scheme(module_palette,length(module_num),show = FALSE)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1,1), mar=c(2,2,2,4),pty = "m")
    plot(network_visual_obj$configured_igraph_object,vertex.label=NA,vertex.size=network_visual_obj$parameters$vertice_size,layout=network_visual_obj$vertices_coordinates,vertex.shape="circle",
         mark.groups=show_list,mark.col=t_cols, mark.border="gray")
    legend_text=paste0("Module #",module_num)
    #plot.new()
    aes_col=network_visual_obj$parameters$aes_col
    legend("topright", c(names(aes_col),"None"), col= c(aes_col,"gray"),
           text.col = c(aes_col,"gray"),
           #lty = c(2, -1, 1),
           pch = c(16, 16, 16), trace=TRUE,border="black",bty="n")
    legend("bottomright",legend=legend_text,col=t_cols, bty="n",fill=t_cols,border="gray",xjust = 1)
  }
}
