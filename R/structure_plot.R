#' Visualize microbial community composition structure based on tax summary object
#'
#' Function for visualization of microbial structure with PCAplot, PCoAplot and NMDSplot
#'
#' @param taxobj tax summary objects computed by \code{\link{tax_summary}}.
#' @param taxlevel taxonomy levels used for visualization.Must be one of
#'   c("Domain","Phylum","Class","Order","Family","Genus","Species","Base").
#' @param ptsize Numeric, default: 2. Size of point in plot. See
#'   \code{\link[ggplot2]{geom_point}} for details.
#' @param diagram Character, default: NULL. A character indicating group diagram,
#'   should be in c("ellipse", "stick", "polygon").
#' @param ellipse.level Numeric, default: 0.85. The level at which to draw an ellipse,
#'   or, if type = "euclid", the radius of the circle to be drawn. See
#'   \code{\link[ggplot2]{stat_ellipse}} for details.
#' @param facet_row Numeric, default: NULL. Number of rows when wrap panels. See
#'   \code{\link[ggplot2]{facet_wrap}} for details.
#'
#' @return Microbial structure analysis object.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes labs facet_wrap geom_hline geom_vline geom_point stat_ellipse geom_segment geom_polygon guides scale_color_manual
#' @importFrom dplyr left_join
#' @importFrom vegan adonis2
#' @importFrom magrittr %>%
#' @note
#' 1. Do not use NMDS when warning: In metaMDS(t(inputframe)) :stress is (nearly) zero: you may have insufficient data
#'
#' 2. Ellipse not available when replicates less than 3, please use 'stick' or 'polygon' instead
#' @examples
#' ###data preparation####
#' data("Two_group")
#'
#' ###analysis####
#' set.seed(999)
#' community_structure<- structure_plot(taxobj = Two_group,taxlevel = "Base")
#'  #check output list in console (not run)
#'  ######Output list##
#'  #####Plot#
#'  ####PCAplot:named as('PCA_Plot')(1/3)
#'  ####PCoAplot:named as('PCoA_Plot')(2/3)
#'  ####NMDSplot:named as('NMDS_Plot')(3/3)
#'  #####Analysis object#
#'  ####PCA object:named as('PCA_object')
#'  ####PCoA object:named as('PCoA_object')
#'  ####NMDS object:named as ('NMDS_object')
#'  #####Coordinates dataframe#
#'  ####PCA Coordinates dataframe:named as('PCA_coordinates')
#'  ####PCoA Coordinates dataframe:named as('PCoA_coordinates')
#'  ####NMDS Coordinates dataframe:named as('NMDS_coordinates')
#'  ######Done##
#'  #check PERMANOVA results
#'  community_structure$PERMANOVA_statistics
#'
#'  #extract plot
#'  community_structure$PCA_Plot
#'  community_structure$PCoA_Plot
#'  community_structure$NMDS_Plot
#'
#'  #extract object
#'  PCA_obj<- community_structure$PCA_object
#'  print(PCA_obj)
#'
#'  #extract coordinates frame
#'  PCA_coord<- community_structure$PCA_coordinates
#'  head(PCA_coord)
#'
#'  #stick plot
#'  set.seed(999)
#'  community_structure<- structure_plot(taxobj = Two_group,taxlevel = "Base",diagram = "stick")
#'  community_structure$PCoA_Plot
#'
#'  #faced form
#'  data("Facet_group")
#'  set.seed(999)
#'  community_structure<- structure_plot(taxobj = Facet_group,taxlevel = "Genus",diagram = "stick")
#'  community_structure$PERMANOVA_statistics
#'  community_structure$PCA_Plot
#'  community_structure$PCoA_Plot
#'  community_structure$NMDS_Plot
structure_plot=function(taxobj,taxlevel,ptsize=2,diagram=NULL,ellipse.level=0.85,facet_row=NULL){
  if(!taxlevel %in% c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species","Base")){
    warning("Illegal 'taxlevel',please choose among c('Domain','Kingdom','Phylum','Class','Order','Family','Genus','Species','Base'")
  }
  if(!"configuration" %in% names(taxobj)){warnings("taxonomy summary object not configured yet, call '?object_config' for configuration")}
  inputframe=eval(parse(text=paste0("taxobj","$",taxlevel,"_percent")))
  input=data.frame(inputframe[,-1],row.names = paste0(taxlevel,1:nrow(inputframe)))
  groupframe=taxobj$Groupfile
  treat_location=taxobj$configuration$treat_location
  facet_location=taxobj$configuration$facet_location
  rep_location=taxobj$configuration$rep_location
  specific.color=taxobj$configuration$treat_col
  output=list()
  #PERMANOVA##
  if(is.null(facet_location)){
    adonis_results=adonis2(t(inputframe[,-1])~groupframe[,treat_location])
    rownames(adonis_results)[1]=colnames(groupframe)[treat_location]
    attr(adonis_results, "heading")[2]='adonis2'
  }else{
    adonis_results=adonis2(t(inputframe[,-1])~groupframe[,treat_location]*groupframe[,facet_location])
    rownames(adonis_results)[1:3]=c(colnames(groupframe)[treat_location],
                                    colnames(groupframe)[facet_location],
                                    paste0(colnames(groupframe)[treat_location],"*",colnames(groupframe)[facet_location]))
    attr(adonis_results, "heading")[2]='adonis2'
  }
  message("###PERMANOVA has done###")
  ##PCA##
  PCA=Dimension_reduction(input,groupframe,1)
  data.pca=PCA$data.pca
  PCA_summary=summary(data.pca)
  PCA_summary_importance=PCA_summary$importance
  PC1_lab=PCA_summary_importance[2,1] %>% round(4)*100
  PC2_lab=PCA_summary_importance[2,2] %>% round(4)*100
  PCAframe<- as.data.frame(PCA[["outframe"]])
  if(is.null(facet_location)){
    cent = aggregate(PCAframe[,1:2],by=list(groupframe[,treat_location]),FUN = mean)
    colnames(cent) = c(colnames(groupframe)[treat_location],'cent1','cent2')
    PCAframe=left_join(PCAframe,cent)%>%suppressMessages()
  }else{
    cent = aggregate(PCAframe[,1:2],by=list(groupframe[,treat_location],groupframe[,facet_location]),FUN = mean)
    colnames(cent) = c(colnames(groupframe)[treat_location],colnames(groupframe)[treat_location],'cent1','cent2')
    PCAframe=left_join(data.frame(PCAframe,joint=paste0(groupframe[,treat_location],groupframe[,facet_location])),data.frame(cent[,3:4],joint=paste0(cent[,1],cent[,2])))%>%suppressMessages()
    PCAframe=PCAframe[,-which(colnames(PCAframe)=="joint")]
  }
  PCAplot=ggplot(data=PCAframe,aes(x =PCAframe[,1],y=PCAframe[,2])) +
    labs(x=paste0("PC1:",PC1_lab,"%"),y=paste0("PC2:",PC2_lab,"%"),title="PCA",fill="Treatment")
  if(is.null(facet_location)==FALSE){
    PCAplot=PCAplot+facet_wrap(~groupframe[,facet_location],nrow=facet_row,scales = "free")
  }else{
    PCAplot=PCAplot+
      geom_hline(yintercept=0,linetype=4,color="grey") +
      geom_vline(xintercept=0,linetype=4,color="grey")
  }
  if(is.null(diagram)){
    PCAplot=PCAplot+
      geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location])))
  }else {if(diagram=="ellipse"){
    PCAplot=PCAplot+
      geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location]))) +
      stat_ellipse(level=ellipse.level,lty=2,aes(color=factor(groupframe[,treat_location])))
  }else if(diagram=="stick"){
    PCAplot=PCAplot+
      geom_segment(aes(xend = PCAframe[,'cent1'], yend = PCAframe[,'cent2'],color=factor(groupframe[,treat_location])), show.legend = FALSE,size=.3,alpha=.8) +
      geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location])))
  }else if(diagram=="polygon"){
    PCAplot=PCAplot+
      geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location]))) +
      geom_polygon(aes(fill=factor(groupframe[,treat_location])),alpha=.8,show.legend = FALSE)
  }else{warning("Illegal input in parameter `diagram`")}}
  PCAplot=PCAplot+theme_zg()+guides(color="none")
  if(is.null(specific.color)){PCAplot=PCAplot}else{PCAplot=PCAplot+scale_color_manual(values = specific.color)+scale_fill_manual(values = specific.color)}
  output=c(output,list(PCAplot,PCA,PCAframe))
  names(output)[1:3]=c("PCA_Plot","PCA_object","PCA_coordinates")
  message("###PCA has done###")
  ##PCoA##
  PCoA=Dimension_reduction(input,groupframe,2)
  PCO1_lab=PCoA$PCOA$values[1,2] %>% round(4)*100
  PCO2_lab=PCoA$PCOA$values[2,2] %>% round(4)*100
  PCoA_data=as.data.frame(PCoA[["outframe"]])
  if(is.null(facet_location)){
    cent = aggregate(PCoA_data[,1:2],by=list(groupframe[,treat_location]),FUN = mean)
    colnames(cent) = c(colnames(groupframe)[treat_location],'cent1','cent2')
    PCoA_data=left_join(PCoA_data,cent)%>%suppressMessages()
  }else{
    cent = aggregate(PCoA_data[,1:2],by=list(groupframe[,treat_location],groupframe[,facet_location]),FUN = mean)
    colnames(cent) = c(colnames(groupframe)[treat_location],colnames(groupframe)[treat_location],'cent1','cent2')
    PCoA_data=left_join(data.frame(PCoA_data,joint=paste0(groupframe[,treat_location],groupframe[,facet_location])),data.frame(cent[,3:4],joint=paste0(cent[,1],cent[,2])))%>%suppressMessages()
    PCoA_data=PCoA_data[,-which(colnames(PCoA_data)=="joint")]
  }
  PCoAplot=ggplot(data=PCoA_data,aes(x =PCoA_data[,1],y=PCoA_data[,2])) +
    labs(x=paste0("Pco1:",PCO1_lab,"%"),y=paste0("Pco2:",PCO2_lab,"%"),title="PCoA",fill="Treatment")
  if(is.null(facet_location)==FALSE){
    PCoAplot=PCoAplot+facet_wrap(~groupframe[,facet_location],nrow=facet_row,scales = "free")
  }else{
    PCoAplot=PCoAplot+
      geom_hline(yintercept=0,linetype=4,color="grey") +
      geom_vline(xintercept=0,linetype=4,color="grey")
  }
  if(is.null(diagram)){
    PCoAplot=PCoAplot+
      geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location])))
  }else{
    if(diagram=="ellipse"){
      PCoAplot=PCoAplot+
        geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location]))) +
        stat_ellipse(level=ellipse.level,lty=2,aes(color=factor(groupframe[,treat_location])))
    }else if(diagram=="stick"){
      PCoAplot=PCoAplot+
        geom_segment(aes(xend = PCoA_data[,'cent1'], yend = PCoA_data[,'cent2'],color=factor(groupframe[,treat_location])), show.legend = FALSE,size=.3,alpha=.8) +
        geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location])))
    }else if(diagram=="polygon"){
      PCoAplot=PCoAplot+
        geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location]))) +
        geom_polygon(aes(fill=factor(groupframe[,treat_location])),alpha=.8,show.legend = FALSE)
    }else{warning("Illegal input in parameter `diagram`")}}
  PCoAplot=PCoAplot+theme_zg()+guides(color="none")
  if(is.null(specific.color)){PCoAplot=PCoAplot}else{PCoAplot=PCoAplot+scale_color_manual(values = specific.color)+scale_fill_manual(values = specific.color)}
  message("###PCoA has done###")
  output=c(output,list(PCoAplot,PCoA,PCoA_data))
  names(output)[4:6]=c("PCoA_Plot","PCoA_object","PCoA_coordinates")
  ##NMDS##
  NMDS=Dimension_reduction(input,groupframe,3)
  NMDSframe=as.data.frame(NMDS[["outframe"]])
  if(is.null(facet_location)){
    cent = aggregate(NMDSframe[,1:2],by=list(groupframe[,treat_location]),FUN = mean)
    colnames(cent) = c(colnames(groupframe)[treat_location],'cent1','cent2')
    NMDSframe=left_join(NMDSframe,cent)%>%suppressMessages()
  }else{
    cent = aggregate(NMDSframe[,1:2],by=list(groupframe[,treat_location],groupframe[,facet_location]),FUN = mean)
    colnames(cent) = c(colnames(groupframe)[treat_location],colnames(groupframe)[treat_location],'cent1','cent2')
    NMDSframe=left_join(data.frame(NMDSframe,joint=paste0(groupframe[,treat_location],groupframe[,facet_location])),data.frame(cent[,3:4],joint=paste0(cent[,1],cent[,2])))%>%suppressMessages()
    NMDSframe=NMDSframe[,-which(colnames(NMDSframe)=="joint")]
  }
  NMDSplot=ggplot(data=NMDSframe,aes(x =NMDSframe[,'MDS1'],y=NMDSframe[,'MDS2'])) +
    labs(x="NMDS1",y="NMDS2",title=paste0("stress",round(NMDS$NMDSstat$stress,2)),fill="Treatment")
  if(is.null(facet_location)==FALSE){
    NMDSplot=NMDSplot+facet_wrap(~groupframe[,facet_location],nrow=facet_row,scales = "free")
  }else{
    NMDSplot=NMDSplot+
      geom_hline(yintercept=0,linetype=4,color="grey") +
      geom_vline(xintercept=0,linetype=4,color="grey")
  }
  if(is.null(diagram)){
    NMDSplot=NMDSplot+
      geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location])))
  }else{
    if(diagram=="ellipse"){
      NMDSplot=NMDSplot+
        geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location]))) +
        stat_ellipse(level=ellipse.level,lty=2)
    }else if(diagram=="stick"){
      NMDSplot=NMDSplot+
        geom_segment(aes(xend = NMDSframe[,'cent1'], yend = NMDSframe[,'cent2'],color=factor(groupframe[,treat_location])), show.legend = FALSE,size=.3,alpha=.8) +
        geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location])))
    }else if(diagram=="polygon"){
      NMDSplot=NMDSplot+
        geom_point(size=ptsize,alpha=.8,pch=21,color="black",aes(fill=factor(groupframe[,treat_location]))) +
        geom_polygon(aes(fill=factor(groupframe[,treat_location])),alpha=.8,show.legend = FALSE)
    }else{warning("Illegal input in parameter `diagram`")}}
  NMDSplot=NMDSplot+theme_zg()+guides(color="none")
  if(is.null(specific.color)){NMDSplot=NMDSplot}else{NMDSplot=NMDSplot+scale_color_manual(values = specific.color)+scale_fill_manual(values = specific.color)}
  output=c(output,list(NMDSplot,NMDS,NMDSframe))
  names(output)[7:9]=c("NMDS_Plot","NMDS_object","NMDS_coordinates")
  message("###NMDS has done###")
  output=c(output,list(adonis_results))
  names(output)[10]="PERMANOVA_statistics"
  ##list##
  ##list##
  message("##Output list##")
  message("#PERMANOVA#")
  message("PERMANOVA:named as('PERMANOVA_statistics')#")
  message("#Plot#")
  message("PCAplot:named as('PCA_Plot')(1/3)")
  message("PCoAplot:named as('PCoA_Plot')(2/3)")
  message("NMDSplot:named as('NMDS_Plot')(3/3)")
  message("#Analysis object#")
  message("PCA object:named as('PCA_object')")
  message("PCoA object:named as('PCoA_object')")
  message("NMDS object:named as ('NMDS_object')")
  message("#Coordinates dataframe#")
  message("PCA Coordinates dataframe:named as('PCA_coordinates')")
  message("PCoA Coordinates dataframe:named as('PCoA_coordinates')")
  message("NMDS Coordinates dataframe:named as('NMDS_coordinates')")
  message("##Done##")
  return(output)
}
