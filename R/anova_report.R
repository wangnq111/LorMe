#' Print Analysis of Variance report
#'
#' @param data Data frame containing the treatment, value and other information.
#' @param treatment_col Numeric indicating where treatment locates (column number) in data.
#' @param value_col Numeric indicating where treatment value (column number) in data.
#' @param prior logical. Whether conducted prior comparisons.
#' @param comparison_method Default would automaticly choose method. Method of multiple comparison,must be one of "SNK", "Tukey", "bonferroni","LSD" or "Scheffe".
#' @param equally_rep Logical. Whether all treatments have same number of replication.
#' @param report Logical. If print report to console. Default:TRUE
#'
#' @return
#' anova_report returns list of:
#'
#' 1)basic data description
#'
#' 2)ANOVA model
#'
#' 3)summary of ANOVA model
#'
#' 4)model of multiple comparison
#'
#' 5)difference of multiple comparison
#'
#' 6)letters of multiple comparison, which could be use for visualization.
#'
#' @export
#'
#' @importFrom agricolae LSD.test scheffe.test SNK.test
#' @importFrom stats aov TukeyHSD sd shapiro.test residuals lm
#' @importFrom multcompView multcompLetters2
#' @importFrom DescTools PostHocTest
#' @importFrom magrittr %$% %>%
#'
#' @examples
#' {
#'   #' Data loading from 'agricolae' package
#'   data("cotton", package = "agricolae")
#'
#'   #' ANOVA report with default settings
#'   anova_results <- anova_report(
#'     data = cotton,
#'     treatment_col = 3,
#'     value_col = 5
#'   )
#'   ## Here returns NULL because no significance among groups
#'
#'   ## To conduct prior comparisons
#'   anova_results <- anova_report(
#'     data = cotton,
#'     treatment_col = 3,
#'     value_col = 5,
#'     prior = TRUE
#'   )
#'
#'   ## Here found no difference among groups, thus change to a more sensitive method
#'   ## (maybe illegal, but only as an example)
#'   anova_results <- anova_report(
#'     data = cotton,
#'     treatment_col = 3,
#'     value_col = 5,
#'     prior = TRUE,
#'     comparison_method = "LSD"
#'   )
#'
#'   #' Data loading 'iris' dataset
#'   data("iris")
#'
#'   #' ANOVA report for 'iris' dataset
#'   anova_results <- anova_report(
#'     data = iris,
#'     treatment_col = 5,
#'     value_col = 2
#'   )
#'
#'   ### Extract return
#'
#'   ### Basic data description
#'   print(anova_results$basicdata)
#'
#'   ### ANOVA model
#'   print(anova_results$anova_model)
#'
#'   ### Summary of ANOVA model
#'   print(anova_results$anova_summary)
#'
#'   ### Model of multiple comparison
#'   print(anova_results$multiple_comparison_model)
#'
#'   ### Difference of multiple comparison
#'   print(anova_results$comparison_results)
#'
#'   ### Letters of multiple comparison, which could be used for visualization
#'   print(anova_results$comparison_letters)
#' }
anova_report=function(data,treatment_col,value_col,prior=FALSE,comparison_method="Auto",equally_rep=TRUE,report=TRUE){
  data[,treatment_col]=factor(data[,treatment_col])
  ntreat=data[,treatment_col]  %>% unique() %>% length()
  ###basic####
  mean_frame=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN=mean)
  Sd=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN=sd)
  Sd= Sd[,2]
  Treatment_Name=mean_frame$Group.1
  N=table(data[,treatment_col]) %>% as.numeric()
  Mean=mean_frame[,2]
  SEM=Sd/(N^0.5)
  basic_data=data.frame(Treatment_Name,N,Mean,Sd,SEM)
  print(basic_data)

  aovmodel_data=data.frame(value=data[,value_col],treatment=data[,treatment_col])
  aovmodel=aov(value~treatment,data=aovmodel_data)
  aov_results=summary(aovmodel)
  aov_frame=aov_results[[1]] %>% as.data.frame()
  rownames(aov_frame)[1]="Between Groups"
  aov_frame[1,2:4]=round(as.numeric(aov_frame[1,2:4]),3)
  aov_frame[2,2:3]=round(as.numeric(aov_frame[2,2:3]),3)
  aov_frame[3,]=aov_frame[1,]+aov_frame[2,]
  rownames(aov_frame)[3]="Total"
  aov_frame[3,3:5]=""
  aov_frame[2,4:5]=""
  print(aov_frame)
  if(report==TRUE){
    cat("\n\n ###Data overview#### \n\n")
    cat("###One Way Analysis of Variance begin#### \n")
    cat("\n\n###Statistics on","###Dependent Variable:",colnames(data)[value_col], "####\n")
    cat("\n\n###Conclusion####\n")
    if(as.numeric(aov_frame[1,5])<0.05){cat("The differences in the mean values among the treatment groups are greater than would be expected by chance; there is a statistically significant difference  (P = ",aov_frame[1,5],") \n")
    }else{cat("The differences in the mean values among the treatment groups are not great enough to exclude the possibility that the difference is due to random sampling variability; there is not a statistically significant difference  (P = ",aov_frame[1,5],") \n")}
    cat("\n\n")
  }
  if(as.numeric(aov_frame[1,5])<0.05|prior==TRUE){
    if(report==TRUE){cat("###Multiple Comparison Procedures#### \n\n")}
    ###method choose####
    if(comparison_method=="Auto"){
      if(equally_rep==FALSE){
        comparison_method="Scheffe"
      }else{
        if(ntreat>3){
          comparison_method="Tukey"
        }else{
          if(as.numeric(aov_frame[1,5])>0.001){
            comparison_method="LSD"
          }else{
            comparison_method="bonferroni"
          }
        }
      }
    }else{comparison_method=comparison_method}
    ###multicomp####
    if(comparison_method %in% c("LSD","bonferroni")){
      ###LSD&bonferroni####
      if(comparison_method=="LSD"){
        LSDresults=LSD.test(aovmodel,"treatment")
        comparison=PostHocTest(aovmodel,"treatment",method="lsd")
        comparison=as.data.frame(comparison$treatment)
        if(report==TRUE){
          cat("All Pairwise Multiple Comparison Procedures (Fisher LSD Method) with alpha=0.05 \n\n")
          cat("###LSD statistics#### \n\n")
        }
      }else{
        LSDresults=LSD.test(aovmodel,"treatment",p.adj = "bonferroni")
        comparison=PostHocTest(aovmodel,"treatment",method="bonferroni")
        comparison=as.data.frame(comparison$treatment)
        if(report==TRUE){
          cat("All Pairwise Multiple Comparison Procedures (Bonferroni Method) with alpha=0.05 \n\n")
          cat("###Bonferroni adjusted LSD statistics#### \n\n")
        }
      }
      LSDstatistics=LSDresults$statistics

      ordered_frame=LSDresults$means
      ordered_frame=ordered_frame[order(ordered_frame$value,decreasing = TRUE),]
      letter=data.frame(compare=rownames(ordered_frame),Letters=LSDresults$groups[rownames(ordered_frame),2],type=colnames(data)[treatment_col],ordered_frame)
      colnames(letter)[c(4,6)]=c("Mean","n")
      if(report==TRUE){
        print(LSDstatistics)
        cat("\n\n###Comparions on ",colnames(data)[treatment_col],"####\n\n")
        print(comparison)
        cat("\n\n###Labels####\n\n")
        print(letter)
      }
      list_return=list(basic_data,aovmodel,aov_frame,LSDresults,comparison,letter)
    }else if(comparison_method =="SNK"){
      ###SNK####
      SNKresults=SNK.test(aovmodel,trt="treatment",console = FALSE)
      comparison=PostHocTest(aovmodel,"treatment",method="newmankeuls")
      comparison= as.data.frame(comparison$treatment)
      SNKstatistics=SNKresults$statistics
      ordered_frame=SNKresults$means
      ordered_frame= ordered_frame[order(ordered_frame$value,decreasing = TRUE),]
      letter=data.frame(compare= rownames(ordered_frame),Letters=SNKresults$groups[,2],type=colnames(data)[treatment_col],ordered_frame)
      colnames(letter)[c(4,6)]=c("Mean","n")
      if(report==TRUE){
        cat("All Pairwise Multiple Comparison Procedures (Student-Newman-Keuls Method)) with alpha=0.05 \n\n")
        cat("###Student-Newman-Keuls statistics#### \n\n")
        print(SNKstatistics)
        cat("\n\n###Comparions on ",colnames(data)[treatment_col],"####\n\n")
        print(comparison)
        cat("\n\n###Labels####\n\n")
        print(letter)
        }
      list_return=list(basic_data,aovmodel,aov_frame,SNKresults,comparison,letter)
    }else if(comparison_method =="Tukey"){
      ###TukeyHSD####
      sub_dat=data.frame(value=data[,value_col],treatment=data[,treatment_col])
      HSDresults=TukeyHSD(aovmodel,"treatment")
      comparison=PostHocTest(aovmodel,"treatment",method="hsd")
      comparison=as.data.frame(comparison$treatment)
      tuk=multcompLetters2(value~treatment,HSDresults$treatment[,"p adj"],data=sub_dat)
      Tukey.labels <- data.frame(tuk['Letters'], stringsAsFactors = FALSE)
      Tukey.labels$compare = rownames(Tukey.labels)
      Tukey.labels$type <- colnames(data)[treatment_col]
      mean_sd <- merge(aggregate(sub_dat[,"value"],by=list(sub_dat[,"treatment"]),FUN=mean),
                       aggregate(sub_dat[,"value"],by=list(sub_dat[,"treatment"]),FUN=sd),by="Group.1"
      )
      names(mean_sd) <- c('compare','Mean','std')
      letter=rbind(merge(Tukey.labels,mean_sd,by='compare'))
      if(report==TRUE){
        cat("All Pairwise Multiple Comparison Procedures (TukeyHSD Method)) with alpha=0.05 \n\n")
        cat("###TukeyHSD statistics#### \n\n")
        print(HSDresults)
        cat("\n\n###Comparions on ",colnames(data)[treatment_col],"####\n\n")
        print(HSDresults$treatment)
        cat("\n\n###Labels####\n\n")
        print(letter)
      }
      list_return=list(basic_data,aovmodel,aov_frame,HSDresults,comparison,letter)
    }else if(comparison_method =="Scheffe"){
      ###Scheffe####
      SCFresults=scheffe.test(aovmodel,trt="treatment",console = FALSE)
      comparison=PostHocTest(aovmodel,"treatment",method="scheffe")
      comparison=as.data.frame(comparison$treatment)
      SCFstatistics=SCFresults$statistics
      ordered_frame=SCFresults$means
      ordered_frame= ordered_frame[order(ordered_frame$value,decreasing = TRUE),]
      letter=data.frame(compare=rownames(ordered_frame),Letters=SCFresults$groups[,2],type=colnames(data)[treatment_col],ordered_frame)
      colnames(letter)[c(4,6)]=c("Mean","n")
      if(report==TRUE){
        cat("All Pairwise Multiple Comparison Procedures (Scheffe Method)) with alpha=0.05 \n\n")
        cat("###Student-Newman-Keuls statistics#### \n\n")
        print(SCFstatistics)
        cat("\n\n###Comparions on ",colnames(data)[treatment_col],"####\n\n")
        print(comparison)
        cat("\n\n###Labels####\n\n")
        print(letter)
      }
      list_return=list(basic_data,aovmodel,aov_frame,SCFresults,comparison,letter)
    }else{warning("Wrong method!!! \n\n")}
    names(list_return)=c("basicdata","anova_model","anova_summary","muiltiple_comparision_model","comparision_results","comparison_letters")
  }else{
    list_return=list(basic_data,aovmodel,aov_frame)
    names(list_return)=c("basicdata","anova_model","anova_summary")
  }
  return(list_return)
}
