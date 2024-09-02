#' Print Wilcoxon Rank Sum and Signed Rank Tests reprot
#'
#' @param data Data frame containing the treatment, value and other information.
#' @param treatment_col Numeric indicating where treatment locates (column number) in data.
#' @param value_col Numeric indicating where treatment value (column number) in data.
#' @param paired  Logical indicating whether you want a paired test.
#' @param subject_col Only meaningful when Pair is ture. Numeric indicating where subject of treatment (column number) in data.
#' @param report Logical. If print report to console. Default:TRUE
#'
#' @return wilcox_test_report returns data frame of basic data description.
#'
#' @export
#'
#' @importFrom stats wilcox.test
#' @importFrom coin wilcoxsign_test
#'
#' @examples
#' {
#'   # Data preparation
#'   testdata <- data.frame(
#'     treatment = c(rep("A", 6), rep("B", 6)),
#'     subject = rep(1:6, 2),
#'     value = c(rnorm(6, 2), rnorm(6, 1))
#'   )
#'
#'   # Wilcoxon test (unpaired)
#'   wilcox_result <- wilcox_test_report(
#'     data = testdata,
#'     treatment_col = 1,
#'     value_col = 3
#'   )
#'
#'   # Wilcoxon signed rank test (paired)
#'   wilcox_result <- wilcox_test_report(
#'     data = testdata,
#'     treatment_col = 1,
#'     value_col = 3,
#'     paired = TRUE,
#'     subject_col = 2
#'   )
#'
#'   ### Basic data description ###
#'   print(wilcox_result)
#' }
wilcox_test_report=function(data,treatment_col,value_col,paired=FALSE,subject_col=NULL,report=TRUE){
  firstq=function(x){
    summaryx=summary(x)
    summaryx[2]%>% return()
  }
  thirdq=function(x){
    summaryx=summary(x)
    summaryx[4]%>% return()}
  Median=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN="median")
  Q1=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN=firstq)
  Q1=Q1[,2]
  Q3=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN=thirdq)
  Q3=Q3[,2]
  N=table(data[,treatment_col]) %>% as.numeric()
  basic_data=data.frame(Treatment_Name=Median[,1],N,Median=Median[,2],Q1,Q3)
  print(basic_data)
  if(report==TRUE){cat("\n\n###Dependent Variable:",colnames(data)[value_col],"####\n\n")}
  if(paired==TRUE){

    data_diff=data[order(data[,treatment_col],data[,subject_col]),]
    rank1=data_diff[,value_col][1:(nrow(data_diff)/2)]
    rank2=data_diff[,value_col][(nrow(data_diff)/2+1):nrow(data_diff)]
    wilcox_results=wilcox.test(rank1,rank2,paired = TRUE,correct = FALSE,exact = FALSE)
    p_exact=wilcox_results$p.value
    p_estimate=wilcox.test(rank1,rank2,paired = TRUE,correct = TRUE,exact = FALSE)$p.value
    T1=wilcox.test(rank1,rank2,paired = TRUE,correct = FALSE,exact = FALSE)$statistic
    T2=wilcox.test(rank2,rank1,paired = TRUE,correct = FALSE,exact = FALSE)$statistic
    Zsummary=coin::wilcoxsign_test(rank1 ~ rank2, distribution = "exact")
    Z= Zsummary@statistic
    Z= Z@standardizedlinearstatistic
    if(report==TRUE){
      cat("###Wilcoxon Signed Rank Test begin#### \n\n")
      cat("###Statistics#### \n\n")
      cat("W = ",abs(T1-T2)," T+ = ",max(c(T1,T2))," T- = ",min(c(T1,T2))," \n\n")
      cat("Z-Statistic (based on positive ranks) = ",abs(Z),"\n\n")
      cat("Yates continuity correction option applied to calculations.\n\n")
      cat("P_estimate = ",p_estimate," ,P_exact = ",p_exact,"\n\n")
      cat("###Conclusion####\n")
      if(p_exact>0.05){
        cat("The change that occurred with the treatment is greater than would be expected by chance; there is a statistically significant difference  (P = ",p_exact,").\n\n")
      }else{
          cat("The change that occurred with the treatment is not great enough to exclude the possibility that it is due to chance  (P =",p_exact,").\n\n")
        }
    }
  }else{
    wilcox_results=wilcox.test(get(colnames(data)[value_col])~get(colnames(data)[treatment_col]),data=data,correct = FALSE,exact = FALSE)
    p_exact=wilcox_results$p.value
    p_estimate=wilcox.test(get(colnames(data)[value_col])~get(colnames(data)[treatment_col]),data=data,correct = TRUE,exact = FALSE)$p.value
    U=wilcox_results$statistic %>% as.numeric()
    if(report==TRUE){
      cat("###Mann-Whitney Rank Sum Test#### \n\n")
      cat("###Statistics#### \n\n")
      cat("Mann-Whitney U Statistic= ",U," \n\n")
      cat("n(small)= ",min(basic_data$N),", n(big)= ",max(basic_data$N),"\n")
      cat("P_estimate = ",p_estimate," ,P_exact = ",p_exact,"\n\n")
      cat("###Conclusion####\n")
      if(p_exact<0.05){
        cat("The difference in the median values between the two groups is greater than would be expected by chance; there is a statistically significant difference  (P = ",p_exact,").\n\n")
      }else{
        cat("The difference in the median values between the two groups is not great enough to exclude the possibility that the difference is due to random sampling variability; there is not a statistically significant difference  (P = ",p_exact,").\n\n")
        }
      }
    }
outlist=c(list(basic_data),list(wilcox_results))
names(outlist)=c("Data_description","Statistics")
return(outlist)
}
