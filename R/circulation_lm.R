#version 1.1.0
#' Circulation of fitting Linear Models
#' @description  Using circulation to fit linear models between one dependent variable and series of independent variable
#' @param y Dependent variable
#' @param xframe Matrix or data frame of independent variable
#' @param margin A vector of 1 or 2 indicates arrangement of xframe. 1:by rows 2:by columns
#'
#' @return Data frame contains lm statistics of all Independent Variable
#' @export
#' @details if row names(for margin 1) and column names(for margin 2) are not given, ID column of return data frame will be row/column numbers.
#' @note Other arguments used in function lm were set as default. See in  \code{\link{lm}}.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats lm
#' @author  Wang Ningqi<2434066068@qq.com>
#'
#' @examples
#' data(testotu)
#'
#' ###using margin 1, arrange by rows##
#' dep=testotu[1,2:21]
#' in_dep=testotu[-1,2:21]
#' lm_stat<-circulation_lm(y = dep,xframe = in_dep,margin = 1)
#' lm_stat
#'
#' ###using margin 2, arrange by column##
#' dep=testotu[,2]
#' in_dep=testotu[,3:21]
#' lm_stat<-circulation_lm(y = dep,xframe = in_dep,margin = 2)
#' lm_stat
#'
circulation_lm=function(y,xframe,margin){
   y=as.numeric(y)
   pb <- txtProgressBar(style=3)
   i=1;ID=as.character()
   pvalue=as.numeric();Rsquared=as.numeric();adjRsquared=as.numeric();Fstatistic=as.numeric();DF=as.numeric()
   Coefficients_x=as.numeric();Coefficients_Intercept=as.numeric()
   star_time <- Sys.time()
   if(margin==1){
    for (i in 1:nrow(xframe)){
      if(sum(as.numeric(xframe[i,]),na.rm = TRUE)==0){
        pvalue=c(pvalue,NA)
        Rsquared=c(Rsquared,NA)
        adjRsquared=c(adjRsquared,NA)
        Fstatistic=c(Fstatistic,NA)
        DF=c(DF,NA)
        Coefficients_x=c(Coefficients_x,NA)
        Coefficients_Intercept=c(Coefficients_Intercept,NA)
      }else{
      lm_summary=lm(y~as.numeric(xframe[i,])) %>% summary()
      pvalue1=lm_summary$coefficients[2,4]
      Rsquared1=lm_summary$r.squared
      adjRsquared1=lm_summary$adj.r.squared
      Fstatistic1=lm_summary$fstatistic[1]
      DF1=lm_summary$df[2]
      Coefficients_x1=lm_summary$coefficients[2,1]
      Coefficients_Intercept1=lm_summary$coefficients[1,1]
      pvalue=c(pvalue,pvalue1)
      Rsquared=c(Rsquared,Rsquared1)
      adjRsquared=c(adjRsquared,adjRsquared1)
      Fstatistic=c(Fstatistic,Fstatistic1)
      DF=c(DF,DF1)
      Coefficients_x=c(Coefficients_x,Coefficients_x1)
      Coefficients_Intercept=c(Coefficients_Intercept,Coefficients_Intercept1)
      }
      ID=c(ID,rownames(xframe)[i])
      setTxtProgressBar(pb, i/nrow(xframe))
      end_time <- Sys.time()
    }
  }else if(margin==2){
    while(i<=ncol(xframe)){
      if(sum(as.numeric(xframe[,i]),na.rm = TRUE)==0){
        pvalue=c(pvalue,NA)
        Rsquared=c(Rsquared,NA)
        adjRsquared=c(adjRsquared,NA)
        Fstatistic=c(Fstatistic,NA)
        DF=c(DF,NA)
        Coefficients_x=c(Coefficients_x,NA)
        Coefficients_Intercept=c(Coefficients_Intercept,NA)
      }else{
        lm_summary=lm(y~as.numeric(xframe[,i])) %>% summary()
        pvalue1=lm_summary$coefficients[2,4]
        Rsquared1=lm_summary$r.squared
        adjRsquared1=lm_summary$adj.r.squared
        Fstatistic1=lm_summary$fstatistic[1]
        DF1 = lm_summary$df[2]
        Coefficients_x1=lm_summary$coefficients[2,1]
        Coefficients_Intercept1=lm_summary$coefficients[1,1]
        pvalue=c(pvalue,pvalue1)
        Rsquared=c(Rsquared,Rsquared1)
        adjRsquared=c(adjRsquared,adjRsquared1)
        Fstatistic=c(Fstatistic,Fstatistic1)
        DF=c(DF,DF1)
        Coefficients_x=c(Coefficients_x,Coefficients_x1)
        Coefficients_Intercept=c(Coefficients_Intercept,Coefficients_Intercept1)
      }
      ID=c(ID,colnames(xframe)[i])
      setTxtProgressBar(pb, i/ncol(xframe))
      end_time <- Sys.time()
      i=i+1
    }
  }else{warning("marin must be 1 or 2")}
   Sys.sleep(1)
   close(pb)
   data.frame(ID,pvalue,Rsquared,adjRsquared,Fstatistic,DF,Coefficients_x,Coefficients_Intercept)  %>% return()
}
