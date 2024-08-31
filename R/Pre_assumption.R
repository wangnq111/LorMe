#version 1.2.0
#' Automatically significance testing
#' @description Automatically choose significance testing
#'
#' @param data Data frame containing the treatment, value and other information.
#' @param treatment_col Numeric indicating where treatment locates (column number) in data.
#' @param value_col Numeric indicating where treatment value (column number) in data.
#' @param paired  Logical indicating whether you want a paired t-test.
#' @param subject_col Only meaningful when Pair is ture. Numeric indicating where subject of treatment (column number) in data.
#' @param prior logical. Whether conducted prior comparisons.
#' @param comparison_method Character string. Only use for more than 2 treatment. Default would automatically choose method. Method of multiple comparison,must be one of "SNK", "Tukey", "bonferroni","LSD" or "Scheffe".
#' @param equally_rep Logical indicating Whether all treatments have same number of replication.
#' @param output A character string indicating output style. Default: "console", which print the report in console. And "file" is available to output report into text-file.
#' @param output_dir Default:"./". Available only when output="file". The direction of output file.
#' @param filename A character string indicating file name of output file. Only work when output set as 'file'.
#' @param report Logical. If print report to console. Default:TRUE
#'
#' @return
#' auto_signif_test returns results of significant test and print report in console or file. See details in example.
#'
#' See results return in \code{\link{t_test_report}}, \code{\link{wilcox_test_report}}, \code{\link{anova_report}}, \code{\link{kruskal_report}}.
#' @note
#' 1.when choose output="file", once caused error that terminate the program, use 'sink()' to end the written of exist files.
#'
#' 2.Please confirm your data is in format of dataframe, else may cause bug! (e.g. Do not use 'read.xlsx' to load data into tibble format)
#' @export
#'
#' @importFrom HH hov
#'
#' @examples
#' ### Here shows different types of experimental design ###
#' data("cotton", package = "agricolae")
#'
#' ### Two randomly designed groups ###
#' sig_results <- auto_signif_test(
#'   data = cotton,
#'   treatment_col = 1,
#'   value_col = 5
#' )
#'
#' ### Two paired design groups ###
#' sig_results <- auto_signif_test(
#'   data = cotton,
#'   treatment_col = 1,
#'   value_col = 5,
#'   paired = TRUE,
#'   subject_col = 2
#' )
#'
#' ### More than two randomly designed groups ###
#' sig_results <- auto_signif_test(
#'   data = cotton,
#'   treatment_col = 2,
#'   value_col = 5
#' )
#' head(sig_results)  # Check outputs
#'
#' ### Conduct prior comparisons ###
#' sig_results <- auto_signif_test(
#'   data = cotton,
#'   treatment_col = 2,
#'   value_col = 5,
#'   prior = TRUE
#' )
#' head(sig_results)  # Check outputs
#' print(sig_results$basicdata)  # Check statistical summary
#' print(sig_results$anova_model)  # Extract ANOVA model
#' print(sig_results$anova_summary)  # Check ANOVA summary
#' print(sig_results$multiple_comparison_model)  # Extract multiple comparison model
#' print(sig_results$comparison_results)  # Check between-group comparison
#' print(sig_results$comparison_letters)  # Check letters (can be used in visualization)
#'
#' ## Change multiple comparison method (maybe not illegal!!)
#' sig_results <- auto_signif_test(
#'   data = cotton,
#'   treatment_col = 2,
#'   value_col = 5,
#'   prior = TRUE,
#'   comparison_method = "LSD"
#' )
#' head(sig_results)  # Check outputs
#' print(sig_results$comparison_letters)  # Note that letters become different
#'
auto_signif_test=function(data,treatment_col,value_col,paired,subject_col,prior=FALSE,comparison_method=NULL,equally_rep=TRUE,output="console",output_dir="./",filename="auto_signif_test",report=TRUE){
  data=as.data.frame(data)
  if(length(unique(data[,treatment_col]))==1){
    warning("Only one group detected, please check you data or 'treatment_col'")
    return(NULL)
  }
  if(is.null(comparison_method)){comparison_method="Auto"}
  if(output=="file"){
    filename=paste0(output_dir,filename,".txt")
    sink(filename)}
  data[,treatment_col]=factor(data[,treatment_col])
  ##pre_hypo####
  if(report==TRUE){cat("###Distribution hypothesis####\n")}
  if(missing(paired)){paired=FALSE}
  if(missing(subject_col)){subject_col=NULL}
  ntreat=data[,treatment_col]  %>% unique() %>% length()
  if(ntreat==1){warnings("Only one treatment detected!!!")}else
    if(ntreat==2){
      if(paired==FALSE){model=lm(data[,value_col]~data[,treatment_col])
      equal_variance=hov(data[,value_col]~data[,treatment_col])
      normality=shapiro.test(residuals(model))}else{
        #paired differences##
        data_diff=data[order(data[,treatment_col],data[,subject_col]),]
        diff=data_diff[,value_col][1:(nrow(data_diff)/2)]-data_diff[,value_col][(nrow(data_diff)/2+1):nrow(data_diff)]
        normality=shapiro.test(diff)
      }
    }else
      if(ntreat>2){
        if(paired==FALSE){
          model=lm(data[,value_col]~data[,treatment_col])
          equal_variance=hov(data[,value_col]~data[,treatment_col])
        }else{
          #paired differences##
          data_diff=data[order(data[,treatment_col],data[,subject_col]),]
          diff=as.numeric()
          for(i in unique(data_diff[,treatment_col])){
            diff_noi=data_diff[which(data_diff[,treatment_col]!=i),]
            diff1=data_diff[,value_col][which(data_diff[,treatment_col]==i)]*(ntreat-1)-aggregate(diff_noi[,value_col],by=list(diff_noi[,subject_col]),FUN=sum)
            diff1= diff1[,2]
            diff=c(diff,diff1)
          }
          data_diff$diff=diff
          model=lm(data_diff$diff~data_diff[,treatment_col])
          equal_variance=hov(data_diff$diff~data_diff[,treatment_col])
        }
        normality=shapiro.test(residuals(model))
      }
  ##output report part1####
  if(report==TRUE){
    if(normality$p.value>0.05){
      cat("Normality Test (Shapiro-Wilk): Passed (P = ",round(normality$p.value,3),")\n\n")}else{
        cat("Normality Test (Shapiro-Wilk): Failed (P = ",normality$p.value,")\n")
      }
    if(ntreat!=2|paired==FALSE){
      if(equal_variance$p.value>0.05){
        cat("Equal Variance Test (Brown-Forsythe):	Passed	(P = ",round(equal_variance$p.value,3),")\n")}else{
          cat("Equal Variance Test (Brown-Forsythe):	Failed	(P = ",equal_variance$p.value,")\n")
        }
    }
    }
  ####differences####
  if(paired==TRUE){
    if(ntreat==1){warnings("Only one treatment detected!!!")}else
      if(ntreat==2){
        if(normality$p.value>0.05){
          t_test_results<-t_test_report(data = data,treatment_col = treatment_col,value_col = value_col,paired = TRUE,subject_col = subject_col)
          message("Analysis finished")
          if(output=="file"){sink()}
          return(t_test_results)
        }else{
          wilcox_result<-wilcox_test_report(data = data,treatment_col = treatment_col,value_col = value_col,paired = TRUE,subject_col = subject_col)
          message("Analysis finished")
          if(output=="file"){sink()}
          return(wilcox_result)
        }
      }else{
        if(normality$p.value>0.05&equal_variance$p.value>0.05){
          message("Recommendation: use 'Repeated Measures Analysis of Variance'")
        }else{
          message("Recommendation: use 'Friedman Repeated Measures Analysis of Variance on Ranks'")
        }
      }
  }else{
    ##paired==FALSE####
    if(ntreat==1){warnings("Only one treatment detected!!!")}else
      if(ntreat==2){
        if(normality$p.value>0.05&equal_variance$p.value>0.05){
          t_test_results<-t_test_report(data = data,treatment_col = treatment_col,value_col = value_col)
          message("Analysis finished")
          if(output=="file"){sink()}
          return(t_test_results)
        }else{
          wilcox_results<-wilcox_test_report(data = data,treatment_col = treatment_col,value_col = value_col,paired = FALSE)
          message("Analysis finished")
          if(output=="file"){sink()}
          return(wilcox_results)
        }
      }else{
        ###ntreat>3###
        if(normality$p.value>0.05&equal_variance$p.value>0.05){
          aov_results<-anova_report(data = data,treatment_col = treatment_col,value_col = value_col,prior=prior,comparison_method=comparison_method,equally_rep=equally_rep)
          message("Analysis finished")
          if(output=="file"){sink()}
          return(aov_results)
        }else{
          kruskal_results<-kruskal_report(data = data,treatment_col = treatment_col,value_col = value_col,prior=prior,comparison_method=comparison_method)
          message("Analysis finished")
          if(output=="file"){sink()}
          return(kruskal_results)
        }
      }
  }
}
