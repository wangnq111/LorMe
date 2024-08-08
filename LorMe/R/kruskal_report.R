#' Print Kruskal-Wallis Rank Sum Test report
#'
#' @param data Data frame containing the treatment, value and other information.
#' @param treatment_col Numeric indicating where treatment locates (column number) in data.
#' @param value_col Numeric indicating where treatment value (column number) in data.
#' @param prior logical. Whether conducted prior comparisons.
#' @param comparison_method Default would automaticly choose method. Method of multiple comparison,must be one of "SNK" or "Tukey".
#' @param equally_rep Logical. Whether all treatments have same number of replication.
#'
#' @importFrom stats kruskal.test
#' @importFrom agricolae LSD.test SNK.test
#' @importFrom dplyr left_join
#' @importFrom multcompView multcompLetters2
#'
#' @return
#' kruskal_report returns list with
#'
#' 1)basic data description
#'
#' 2)summary of Kruskal-Wallis Rank Sum Test
#'
#' 3)model of multiple comparison
#'
#' 4)difference of multiple comparison
#'
#' 5)letters of multiple comparison, which could be use for visualization.
#'
#' @export
#'
#' @examples
#'data("cotton",package ="agricolae" )
#'kruskal_results=kruskal_report(data = cotton,treatment_col =3,value_col = 5)
#'##here returns NULL because no significance among groups
#'
#'##to conduct prior comparisons.
#'kruskal_results=kruskal_report(data = cotton,treatment_col =3,value_col = 5,prior = TRUE)
#'
#'
#'data("iris")
#'kruskal_results=kruskal_report(data = iris,treatment_col = 5,value_col = 2)
#'
#'###extract return##
#'
#'###basic data description
#'kruskal_results$basicdata
#'
#'###summry of Kruskal-Wallis Rank Sum Test
#'kruskal_results$Kruskal_Wallis_summary
#'
#'###model of multiple comparision
#'kruskal_results$muiltiple_comparision_model
#'
#'###difference of multiple comparision
#'kruskal_results$comparision_results
#'
#'###letters of multiple comparision, which could be use for visualization.
#'kruskal_results$comparison_letters
#'
kruskal_report=function(data,treatment_col,value_col,prior=F,comparison_method="Auto",equally_rep=T){
  ###basic####
  data[,treatment_col]=factor(data[,treatment_col])
  ntreat=data[,treatment_col]  %>% unique() %>% length()
  mean_frame=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN=mean)
  Sd=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN=sd)
  Sd=Sd[,2]
  Treatment_Name=mean_frame$Group.1
  N=table(data[,treatment_col]) %>% as.numeric()
  Mean=mean_frame[,2]
  SEM=Sd/(N^0.5)
  basic_data=data.frame(Treatment_Name,N,Mean,Sd,SEM)
  cat("\n\n ###Data overview#### \n\n")
  print(basic_data)
  cat("###Kruskal-Wallis One Way Analysis of Variance on Ranks#### \n")
  kru_results=kruskal.test(get(colnames(data)[value_col])~get(colnames(data)[treatment_col]),data=data)
  pvalue=kru_results$p.value
  cat("H = ",kru_results$statistic," with ",as.numeric(kru_results$parameter)," degrees of freedom. P = ",pvalue," \n\n")
  if(pvalue<0.05){
    cat("The differences in the median values among the treatment groups are greater than would be expected by chance; there is a statistically significant difference  (P = ",pvalue,") \n\n")
  }else{
    cat("The differences in the median values among the treatment groups are not great enough to exclude the possibility that the difference is due to random sampling variability; there is not a statistically significant difference    (P = ",pvalue,") \n\n")
  }
  aovmodel=data.frame(value=data[,value_col],treatment=data[,treatment_col])
  aovmodel=aov(value~treatment,data=aovmodel)
  aov_results=summary(aovmodel)
  aov_frame=aov_results[[1]] %>% as.data.frame()
  ###multicomparison####
  if(pvalue<0.05|prior==T){
    cat("###Multiple Comparison Procedures#### \n\n")
    ###method choose####
    if(comparison_method=="Auto"){
      if(equally_rep==F){
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

  aovmodel=data.frame(value=data[,value_col],treatment=data[,treatment_col])
  aovmodel= aov(value~treatment,data=aovmodel)
  if(comparison_method =="SNK"){
    ###SNK####
    cat("All Pairwise Multiple Comparison Procedures (Student-Newman-Keuls Method)) with alpha=0.05 \n\n")
    SNKresults=SNK.test(aovmodel,trt="treatment",console = F)
    cat("###Student-Newman-Keuls statistics#### \n\n")
    comparison=PostHocTest(aovmodel,"treatment",method="newmankeuls")
    comparison=as.data.frame(comparison$treatment)
    SNKstatistics=SNKresults$statistics
    print(SNKstatistics)
    cat("\n\n###Comparions on ",colnames(data)[treatment_col],"####\n\n")
    print(comparison)
    ordered_frame=SNKresults$means
    ordered_frame=ordered_frame[order(ordered_frame$value,decreasing = T),]
    letter=data.frame(compare= rownames(ordered_frame),Letters=SNKresults$groups[,2],type=colnames(data)[treatment_col],ordered_frame)
    colnames(letter)[c(4,6)]=c("Mean","n")
    cat("\n\n###Labels####\n\n")
    print(letter)
    list_return=list(basic_data,kru_results,SNKresults,comparison,letter)
  }else if(comparison_method =="Tukey"){
    ###TukeyHSD####
    cat("All Pairwise Multiple Comparison Procedures (TukeyHSD Method)) with alpha=0.05 \n\n")
    sub_dat=data.frame(value=data[,value_col],treatment=data[,treatment_col])
    comparison=PostHocTest(aovmodel,"treatment",method="hsd")
    comparison=as.data.frame(comparison$treatment)
    cat("###TukeyHSD statistics#### \n\n")
    HSDresults=TukeyHSD(aovmodel,"treatment")
    cat("\n\n###Comparions on ",colnames(data)[treatment_col],"####\n\n")
    print(HSDresults$treatment)
    cat("\n\n###Labels####\n\n")
    tuk=multcompLetters2(value~treatment,HSDresults$treatment[,"p adj"],data=sub_dat)
    Tukey.labels <- data.frame(tuk['Letters'], stringsAsFactors = FALSE)
    Tukey.labels$compare = rownames(Tukey.labels)
    Tukey.labels$type <- colnames(data)[treatment_col]
    mean_sd <- merge(aggregate(sub_dat[,"value"],by=list(sub_dat[,"treatment"]),FUN=mean),
                     aggregate(sub_dat[,"value"],by=list(sub_dat[,"treatment"]),FUN=sd),by="Group.1"
    )
    names(mean_sd) <- c('compare','Mean','std')
    letter=rbind(merge(Tukey.labels,mean_sd,by='compare'))
    print(letter)
    list_return=list(basic_data,kru_results,HSDresults,comparison,letter)
  }else if(comparison_method %in% c("LSD","bonferroni")){
    ###LSD&bonferroni####
    if(comparison_method=="LSD"){
      cat("All Pairwise Multiple Comparison Procedures (Fisher LSD Method) with alpha=0.05 \n\n")
      LSDresults=LSD.test(aovmodel,"treatment")
      cat("###LSD statistics#### \n\n")
      comparison=PostHocTest(aovmodel,"treatment",method="lsd")
      comparison=as.data.frame(comparison$treatment)
    }else{
      cat("All Pairwise Multiple Comparison Procedures (Bonferroni Method) with alpha=0.05 \n\n")
      LSDresults=LSD.test(aovmodel,"treatment",p.adj = "bonferroni")
      cat("###Bonferroni adjusted LSD statistics#### \n\n")
      comparison=PostHocTest(aovmodel,"treatment",method="bonferroni")
      comparison=as.data.frame(comparison$treatment)
    }
    LSDstatistics=LSDresults$statistics
    print(LSDstatistics)
    cat("\n\n###Comparions on ",colnames(data)[treatment_col],"####\n\n")
    print(comparison)
    ordered_frame=LSDresults$means
    ordered_frame=ordered_frame[order(ordered_frame$value,decreasing = T),]
    letter=data.frame(compare=rownames(ordered_frame),Letters=LSDresults$groups[rownames(ordered_frame),2],type=colnames(data)[treatment_col],ordered_frame)
    colnames(letter)[c(4,6)]=c("Mean","n")
    cat("\n\n###Labels####\n\n")
    print(letter)
    list_return=list(basic_data,aovmodel,LSDresults,comparison,letter)

  }
  names(list_return)=c("basicdata","Kruskal_Wallis_summary","muiltiple_comparision_model","comparision_results","comparison_letters")
  return(list_return)
  }else{
    list_return=list(basic_data,kru_results)
    names(list_return)=c("basicdata","Kruskal_Wallis_summary")
    return(list_return)
  }
}
