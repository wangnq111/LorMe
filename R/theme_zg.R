#' A classic theme for ggplot
#'
#' @return ggplot theme
#' @export
#' @note Build inside the LorMe package, Please use theme_zg() as a theme directly
#' @importFrom ggplot2 theme element_rect element_blank unit element_line element_text


theme_zg <- function(){
    theme(panel.background = element_rect(fill="transparent",color="black",size=.4),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.2,"lines"),
          axis.ticks = element_line(color='black',size=.2),
          axis.line = element_blank(),
          axis.title=element_text(colour='black', size=10,face = "bold"),
          axis.text=element_text(colour='black',size=8),
          strip.background = element_rect(fill='transparent',size=.4,colour = "transparent"),
          strip.text.x = element_text(size=8),#, face = "bold"
          legend.title = element_text(colour = "black",size=10,face="bold",hjust=.5),
          legend.text =  element_text(colour = "black",size=8,face="bold"),
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.background = element_rect(fill="transparent"),
          plot.title =  element_text(colour = "black",size=10,face="bold"))
}


