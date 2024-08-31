#' Comparison plot generator
#' This function help generate comparsion plot including bar plot, box plot, and violin plot
#'
#' @param inputframe A data frame contain information for visualization.
#' @param treat_location Numeric. Treatment column number in inputframe.
#' @param value_location Numeric. Value column number in inputframe.
#' @param aes_col Named character string, default:NULL. A set of aesthetic character to map treatment to.
#' @param point Logical. If draw point on bar, box and violin plot. Default:TRUE.
#' @param facet_location  Numeric, default:NULL. Facet column number in inputframe.
#' @param ylab_text Character. Text for y axis.
#'
#' @import ggplot2
#' @import ggpubr
#' @importFrom tidyr separate
#' @return A list contained plot and statistics
#' @export
#'
#' @examples
#' data("iris")
#' results=compare_plot(inputframe=iris,treat_location=5,
#'                      value_location=1,ylab_text = "Sepal Length")
#'
#' #Check statistics
#' results$Statistics
#' #Extract plot
#' results$Barplot
#' results$Boxplot
#' results$Violinplot
#'
#'
#' iris$Treat2=rep(c(rep("A",25),rep("B",25)),3)
#'
#' results=compare_plot(inputframe=iris,treat_location=5,
#'                      value_location=1,facet_location = 6,
#'                      ylab_text = "Sepal Length")
#'
#' #Check statistics
#' results$Statistics
#' #Extract plot
#' results$Barplot
#' results$Boxplot
#' results$Violinplot
#' #Extract combined plot
#' results$All_Barplot
#' results$All_Boxplot
#' results$All_Violinplot
#'
compare_plot=function(inputframe,treat_location,value_location,aes_col=NULL,point=TRUE,facet_location=NULL,ylab_text=NULL){
  inputframe=as.data.frame(inputframe)
  if(is.numeric(inputframe[,value_location])==FALSE){
    warning("'value_location' indicated non-numeric column,please check")
    return(NULL)
  }
  condition=unique(inputframe[,treat_location])
  if(is.null(aes_col)){
    aes_col=color_scheme("Plan7",length(condition),show = FALSE,names = condition)
  }
  if(length(condition)==1){
    warning("Only one Group contained! Invalid 'inputframe'")
    return(NULL)
  }
  #data description
  if(is.null(facet_location)){
    mean_frame<- aggregate(inputframe[,value_location],by=list(inputframe[,treat_location]),FUN=mean)
    Sd<- aggregate(inputframe[,value_location],by=list(inputframe[,treat_location]),FUN=sd)
    Sd<- Sd[,'x']
    Sd[is.na(Sd)]=0
    Treatment_Name<- mean_frame$Group.1
    N<- table(inputframe[,treat_location]) %>% as.numeric()
    Mean<- mean_frame[,'x']
    SEM<- Sd/(N^0.5)
  }else{
    mean_frame<- aggregate(inputframe[,value_location],by=list(inputframe[,treat_location],
                                                               inputframe[,facet_location]),FUN=mean)
    Sd<- aggregate(inputframe[,value_location],by=list(inputframe[,treat_location],
                                                       inputframe[,facet_location]),FUN=sd)
    Sd<- Sd[,'x']
    Sd[is.na(Sd)]=0
    Treatment_Name<- mean_frame$Group.1
    Facet_name=mean_frame$Group.2
    N<- table(inputframe[,c(treat_location,facet_location)]) %>% as.numeric()
    N<-N[N!=0]
    Mean<- mean_frame[,'x']
    SEM<- Sd/(N^0.5)
  }

  ylimitmax=(1.2*Mean+Sd) %>% max()
  ylimitmin=(0.9*Mean-Sd) %>% max()
  ylimitmax2=1.2*max(inputframe[,value_location])
  ylimitmin2=0.9*min(inputframe[,value_location])
  if(is.null(facet_location)){
    input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM)
  }else{
    input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,'Group.2'])
    colnames(input_mean_frame)[6]=colnames(inputframe)[facet_location]
  }
  #statistics
  if(is.null(facet_location)){
    sig_results<-auto_signif_test(data =inputframe,treatment_col =treat_location,value_col =value_location,prior = TRUE)
  }else{
    inputframe1=data.frame(inputframe,combinetag=paste0(inputframe[,treat_location],"_",inputframe[,facet_location]))
    all_sig_results<-auto_signif_test(data =inputframe1,treatment_col =ncol(inputframe1),value_col =value_location,prior = TRUE)
    grouped_sig_results=list()
    for(i in unique(inputframe[,facet_location])){
      subd=inputframe[inputframe[,facet_location]==i,]
      if(nrow(subd)==0){
        sig_results=NA
      }else if(length(unique(subd[,treat_location]))==1){
        sig_results=NA
      }else{
        sig_results<-auto_signif_test(data =subd,treatment_col=treat_location,value_col =value_location,prior = TRUE)
      }
      grouped_sig_results=c(grouped_sig_results,list(sig_results))
      names(grouped_sig_results)[length(grouped_sig_results)]=i
    }
  }

  #visualization
  bar=ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+
    scale_y_continuous(expand = c(0,0),limits = c(0,ylimitmax))+
    labs(x='',fill='',y=ylab_text)+
    geom_bar(stat = 'identity',size=0.2,width=0.5,aes(fill=as.factor(Treatment_Name)),color='#000000',alpha=0.8)+
    geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)
  if(ylimitmin<0){
    bar=bar+
      scale_y_continuous(expand = c(0,0),limits = c(ylimitmin2,ylimitmax2))
  }
  box=ggplot(inputframe,aes(x=as.factor(inputframe[,treat_location]),y=inputframe[,value_location]))+
    scale_y_continuous(expand = c(0,(ylimitmax2-ylimitmin2)/50))+ #,limits = c(ylimitmin,ylimitmax)
    labs(x='',fill='',y=ylab_text)+
    geom_boxplot(aes(fill=factor(inputframe[,treat_location])),alpha=0.8,width=0.5,outlier.color =NA,color='#000000',linewidth=0.2)

  violin=ggplot(inputframe,aes(x=as.factor(inputframe[,treat_location]),y=inputframe[,value_location]))+
    scale_y_continuous(expand = c(0,(ylimitmax2-ylimitmin2)/50))+ #,limits = c(ylimitmin,ylimitmax)
    labs(x='',fill='',y=ylab_text)+
    geom_violin(aes(fill=factor(inputframe[,treat_location])),alpha=0.8,width=0.5,color='#000000',linewidth=0.2)
  #add statistics
  if(length(condition)>2){
    if(is.null(facet_location)){
      letterframe=sig_results$comparison_letters
    }else{
      letterframe=data.frame()
      for(i in names(grouped_sig_results)){
        sig_results=grouped_sig_results[[i]]
        letterframe_temp=data.frame(sig_results$comparison_letters,tag=i)
        colnames(letterframe_temp)[ncol(letterframe_temp)]=colnames(inputframe)[facet_location]
        letterframe=rbind(letterframe,letterframe_temp)
      }
    }
    letterp=(1.1*letterframe$Mean+letterframe$std) %>% max()
    letters<- data.frame(letterframe,letterp)
    bar=bar+
      geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=letters[,'Letters']),size=6)

    box=box+
      geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=letters[,'Letters']),size=6)

    violin=violin+
      geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=letters[,'Letters']),size=6)
  }else{
    bar=bar+
      stat_compare_means(data=inputframe,aes(x=as.factor(inputframe[,treat_location]),y=inputframe[,value_location],label = after_stat(get('p.signif'))),size=4,label.x.npc = 0.5,label.y.npc = .8)

    box=box+
      stat_compare_means(data=inputframe,aes(x=as.factor(inputframe[,treat_location]),y=inputframe[,value_location],label = after_stat(get('p.signif'))),size=4,label.x.npc = 0.5,label.y.npc = .95)

    violin=violin+
      stat_compare_means(data=inputframe,aes(x=as.factor(inputframe[,treat_location]),y=inputframe[,value_location],label = after_stat(get('p.signif'))),size=4,label.x.npc = 0.5,label.y.npc = .95)
  }

  if(point==TRUE){
    bar=bar+
      geom_jitter(data=inputframe,pch=21,
                  aes(x=factor(inputframe[,treat_location]),y=inputframe[,value_location],fill=factor(inputframe[,treat_location])),
                  size=1,alpha=.8,show.legend = FALSE,width = 0.2)
    box=box+
      geom_jitter(pch=21,
                  aes(fill=factor(inputframe[,treat_location])),
                  size=1,alpha=.8,show.legend = FALSE,width = 0.2)
    violin=violin+
      geom_jitter(pch=21,
                  aes(fill=factor(inputframe[,treat_location])),
                  size=1,alpha=.8,show.legend = FALSE,width = 0.2)
  }
  if(!is.null(facet_location)){
    bar=bar+
      facet_wrap(~get(colnames(inputframe)[facet_location]),scales="fixed")
    box=box+
      facet_wrap(~get(colnames(inputframe)[facet_location]),scales="fixed")
    violin=violin+
      facet_wrap(~get(colnames(inputframe)[facet_location]),scales="fixed")
  }
  bar=bar+theme_zg()+
    #scale_color_manual(values=aes_col)+
    scale_fill_manual(values=aes_col)
  box=box+theme_zg()+
    #scale_color_manual(values=aes_col)+
    scale_fill_manual(values=aes_col)
  violin=violin+theme_zg()+
    #scale_color_manual(values=aes_col)+
    scale_fill_manual(values=aes_col)

  if(is.null(facet_location)){
    outlist=c(list(sig_results),list(bar),list(box),list(violin),list(aes_col))
    names(outlist)=c("Statistics","Barplot","Boxplot","Violinplot","aes_color")
  }else{
    letterframe=all_sig_results$comparison_letters
    tempinto=colnames(inputframe)[c(treat_location,facet_location)]
    letterframe=tidyr::separate(data = letterframe,col="compare",into =tempinto ,sep = "_")
    letterp=(1.1*letterframe$Mean+letterframe$std) %>% max()
    letters<- data.frame(letterframe,letterp)

    barall=ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean,fill=as.factor(input_mean_frame[,ncol(input_mean_frame)])))+
      scale_y_continuous(expand = c(0,0),limits = c(0,ylimitmax))+
      labs(x='',fill='',y=ylab_text)+
      geom_bar(stat = 'identity',position = position_dodge(width = 0.5),size=0.2,width=0.5,color='#000000',alpha=0.8)+
      geom_errorbar(position = position_dodge(0.5),aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)+
      geom_text(data=letters,aes(x=as.factor(letters[,1]),y=letterp,label=letters[,'Letters'],group=as.factor(letters[,2])),size=6,position=position_dodge(0.5))+
      theme_zg()
    if(ylimitmin<0){
      barall=barall+
        scale_y_continuous(expand = c(0,0),limits = c(ylimitmin2,ylimitmax2))
    }
    boxall=ggplot(inputframe,aes(x=as.factor(inputframe[,treat_location]),y=inputframe[,value_location]))+
      scale_y_continuous(expand = c(0,(ylimitmax2-ylimitmin2)/50))+ #,limits = c(ylimitmin,ylimitmax)
      labs(x='',fill='',y=ylab_text)+
      geom_boxplot(aes(fill=as.factor(inputframe[,facet_location])),alpha=0.8,width=0.5,outlier.color =NA,color='#000000',linewidth=0.2)+
      geom_text(data=letters,aes(x=as.factor(letters[,1]),y=letterp,label=letters[,'Letters'],group=as.factor(letters[,2])),size=6,position=position_dodge(0.5))+
      theme_zg()

    violinall=ggplot(inputframe,aes(x=as.factor(inputframe[,treat_location]),y=inputframe[,value_location]))+
      scale_y_continuous(expand = c(0,(ylimitmax2-ylimitmin2)/50))+ #,limits = c(ylimitmin,ylimitmax)
      labs(x='',fill='',y=ylab_text)+
      geom_violin(aes(fill=as.factor(inputframe[,facet_location])),alpha=0.8,width=0.5,color='#000000',linewidth=0.2)+
      geom_text(data=letters,aes(x=as.factor(letters[,1]),y=letterp,label=letters[,'Letters'],group=as.factor(letters[,2])),size=6,position=position_dodge(0.5))+
      theme_zg()
    sig_results=c(list(grouped_sig_results),list(all_sig_results))
    outlist=c(list(sig_results),list(bar),list(box),list(violin),list(aes_col),
              list(barall),list(boxall),list(violinall))
    names(outlist)=c("Statistics","Barplot","Boxplot","Violinplot","aes_color",
                     "All_Barplot","All_Boxplot","All_Violinplot")
  }
  return(outlist)
}

