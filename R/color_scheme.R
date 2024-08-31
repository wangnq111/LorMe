#' Get color scheme
#' @description
#' color_scheme() can generate color scheme from nine color scheme database and expand into colorRamp
#'
#' @param Plan Character, 'Plan1' to 'Plan10' are optional.
#' @param expand Numeric, default:NULL. Numeric indicating numbers to expand color scheme into colorRamp
#' @param names Character string. Names to assign for color scheme.
#' @param show Logical. If show assigned color in plot panel. Default:TRUE.
#'
#' @return If parameter 'names' is not given, 'color_scheme' returns character string including color scheme.When 'names' is set,'color_scheme' returns named vector of color scheme.
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom scales show_col
#'
#' @note
#' 1.Parameter 'names' is strongly recommended to assign for fixed color scheme, see details in ggplot::scale_color_manual
#' @examples
#' ### Commonly used example ###
#' my_color <- color_scheme(
#'   Plan = "Plan1",
#'   names = c("Treatment1", "Treatment2")
#' )
#'
#' ### Generate colorRamp still based on 'Plan1'
#' my_color <- color_scheme(
#'   Plan = "Plan1",
#'   expand = 4,
#'   names = c("Treatment1", "Treatment2", "Treatment3", "Treatment4")
#' )
#'
#' ### View color scheme from plan1 to plan10 in 'Plots' interface ###
#' color_scheme(Plan = "Plan1")
#' color_scheme(Plan = "Plan2")
#' color_scheme(Plan = "Plan3")
#' color_scheme(Plan = "Plan4")
#' color_scheme(Plan = "Plan5")
#' color_scheme(Plan = "Plan6")
#' color_scheme(Plan = "Plan7")
#' color_scheme(Plan = "Plan8")
#' color_scheme(Plan = "Plan9")
#' color_scheme(Plan = "Plan10")
color_scheme=function(Plan,expand=NULL,names=NULL,show=TRUE){
  if(inherits(Plan, "character") && length(Plan) == 1){
  if(Plan=="Plan1"){
    color_Plan=c("#E69F00","#56B4E9")
  }else if(Plan=="Plan2"){
    color_Plan=c("#FE5D5D","#71C9DD","#33B39F","#6376A0","#F5AF98")
  }else if(Plan=="Plan3"){
    color_Plan=c("#4070AF", "#8CA5BB", "#D8DDC7","#FCD39B", "#F18159", "#D42C24")
  }else if(Plan=="Plan4"){
    color_Plan=c("#35A585","#EAE48E","#006FB0","#CC78A6","#F2C661","#56B4E9")
  }else if(Plan=="Plan5"){
    color_Plan=c("#f49128","#194a55","#187c65","#f26115","#c29f62","#83ba9e")
  }else if(Plan=="Plan6"){
    color_Plan=c("#c62d17","#023f75","#ea894e","#266b69","#eb4601","#f6c619")
  }else if(Plan=="Plan7"){
    color_Plan=c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")
  }else if(Plan=="Plan8"){
    color_Plan=c("#E64B35FF","#4DBBD5FF","#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF","#B09C85FF")
  }else if(Plan=="Plan9"){
    color_Plan=c("#6a73cf","#edd064","#0eb0c8","#f2ccac","#a1d5b9","#e1abbc","#fa6e01","#2f2f2f","#972b1d","#e6a84b","#4c211b","#ff717f")
  }else if(Plan=="Plan10"){
    color_Plan=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99" ,"#E31A1C" ,"#FDBF6F", "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99","#B15928")
  }else{warning("Please choose correct Plan! ('Plan1' to 'Plan9')")}
  }else{
    color_Plan=Plan
  }
  if(!is.null(expand)){
    assign_col=colorRampPalette(color_Plan)(expand)
  }else{
    assign_col=color_Plan
  }
  if(length(names)>length(assign_col)){
    assign_col=colorRampPalette(color_Plan)(length(names))
  }
  names(assign_col)=names
  if(isTRUE(show)){
    show_col(assign_col)
    message("Color scheme generated, see in your plot interface")
    }
  return(assign_col)
}


