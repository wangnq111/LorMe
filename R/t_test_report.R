#' Print Student's t-Test report
#'
#' @param data Data frame containing the treatment, value and other information.
#' @param treatment_col Numeric indicating where treatment locates (column number) in data.
#' @param value_col Numeric indicating where treatment value (column number) in data.
#' @param paired Logical indicating whether you want a paired t-test.
#' @param subject_col Only meaningful when Pair is ture. Numeric indicating where subject of treatment (column number) in data.
#' @param report Logical. If print report to console. Default:TRUE
#'
#' @return
#' t_test_report returns list containing:
#' 1) data frame of basic data descrption
#' 2) results of student's t-Test
#'
#' @import magrittr
#' @export
#'
#' @examples
#' {
#'   ### Data preparation ###
#'   testdata <- data.frame(
#'     treatment = c(rep("A", 6), rep("B", 6)),
#'     subject = rep(c(1:6), 2),
#'     value = c(rnorm(6, 2), rnorm(6, 1))
#'   )
#'
#'   # Perform t-test (unpaired)
#'   t_test_result <- t_test_report(
#'     data = testdata,
#'     treatment_col = 1,
#'     value_col = 3
#'   )
#'
#'   # Perform paired t-test
#'   t_test_result <- t_test_report(
#'     data = testdata,
#'     treatment_col = 1,
#'     value_col = 3,
#'     paired = TRUE,
#'     subject_col = 2
#'   )
#'
#'   ### Basic data description ###
#'   print(t_test_result[[1]])
#'   print(t_test_result$basicdata)
#'
#'   ### T-test results ###
#'   print(t_test_result[[2]])
#'   print(t_test_result$t.test_results)
#' }
t_test_report=function(data,treatment_col,value_col,paired,subject_col,report=TRUE){
  if(missing(paired)){paired=FALSE}
  if(missing(subject_col)){subject_col=NULL}
  if((table(data[,treatment_col]) %>% min())==1){
    stop("Observation not enough!")
  }else{
    if(paired==TRUE){

      data_diff=data[order(data[,treatment_col],data[,subject_col]),]
      diff=data_diff[,value_col][1:(nrow(data_diff)/2)]-data_diff[,value_col][(nrow(data_diff)/2+1):nrow(data_diff)]
      mean_frame=aggregate(data_diff[,value_col],by=list(data_diff[,treatment_col]),FUN=mean)
      sd=aggregate(data_diff[,value_col],by=list(data_diff[,treatment_col]),FUN=sd)
      sd=sd[,2]
      Treatment_Name=c(as.character(mean_frame$Group.1),"Difference")
      N=rep(nrow(data)/2,3)
      Mean=c(mean_frame[,2],mean(diff))
      Sd=c(sd,sd(diff))
      SEM=Sd/(N^0.5)
      results=t.test(Pair(data_diff[,value_col][1:(nrow(data_diff)/2)], data_diff[,value_col][(nrow(data_diff)/2+1):nrow(data_diff)]) ~ 1)
      if(report==TRUE){cat("###Paired t-test begin ####\n\n")}
    }else{
      mean_frame=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN=mean)
      Sd=aggregate(data[,value_col],by=list(data[,treatment_col]),FUN=sd)
      Sd=Sd[,2]
      Treatment_Name=mean_frame$Group.1
      N=table(data[,treatment_col]) %>% as.numeric()
      Mean=mean_frame[,2]
      SEM=Sd/(N^0.5)
      results=t.test(get(colnames(data)[value_col])~get(colnames(data)[treatment_col]),data=data)
      if(report==TRUE){cat("###T-test begin ####\n\n")}
    }
    basic_data=data.frame(Treatment_Name,N,Mean,Sd,SEM)
    if(report==TRUE){
      cat("###Dependent Variable:",colnames(data)[value_col],"####\n\n")
      print(basic_data)
      results$data.name=colnames(data)[value_col]
      results$method=paste0("Method: ",results$method," (Two-tailed)")
      cat("###Statistics####\n")
      print(results)
      cat("###Conclusion####\n")
      if(paired==TRUE){
        if(results$p.value<0.05){cat("The change that occurred with the treatment is greater than would be expected by chance;\n There is a statistically significant change  (P = ",results$p.value,") \n")
        }else{cat("The change that occurred with the treatment is not great enough to exclude the possibility that the difference is due to chance  (P = ",results$p.value,") \n")}
      }else{
        if(results$p.value<0.05){cat("The difference in the mean values of the two groups is greater than would be expected by chance;\n There is a statistically significant difference between the input groups (P = ",results$p.value,") \n")
        }else{cat("The difference in the mean values of the two groups is not great enough to reject the possibility that the difference is due to random sampling variability. There is not a statistically significant difference between the input groups (P = ",results$p.value,") \n")}
      }
    }
  }
  returnlist=list(basic_data,results)
  names(returnlist)=c("basicdata","t.test_results")
  return(returnlist)
}
