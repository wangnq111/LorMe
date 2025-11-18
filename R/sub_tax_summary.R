#' Subsetting tax summary objects
#'
#' @param taxobj tax summary objects computed by \code{\link{tax_summary}}.
#' @param ... 	logical expression that are are defined in terms of the variables in Groupfile of tax summary objects. See details in \code{\link{subset}}.
#' @param specificnum specific numbers indicating samples to keep based on Groupfile of tax summary objects.
#' @param taxnum specific numbers indicating taxonomy to keep based on Base file
#'
#' @return Subset of tax summary objects.Same as \code{\link{tax_summary}}.
#' @author Wang Ningqi <2434066068@qq.com>
#' @export
#'
#' @importFrom magrittr %T>%
#'
#' @examples
#'   data("Three_group")
#'
#'   # Check meta file
#'   print(Three_group$groupfile)
#'
#'   # Subsetting tax summary objects
#'
#'   # Select BF and OF groups
#'   sub_testtax_summary <- sub_tax_summary(Three_group, Group %in% c("BF", "OF"))
#'   print(sub_testtax_summary$groupfile)
#'

sub_tax_summary=function(taxobj,...,specificnum=NULL,taxnum=NULL){
  ## 1. data loading
  if(is.null(specificnum)){
    sub_group=subset(methods::slot(taxobj,"groupfile"),...)
    sub_rownames=rownames(sub_group)
    select_num=which(rownames(methods::slot(taxobj,"groupfile")) %in% sub_rownames)
  }else{select_num=specificnum}

  ## 2. select function
  select_temp=function(x){
    if(is.numeric(x[,-1]%>% as.matrix())==TRUE){x[,c(1,select_num+1)]}else(x=x)
  }

  ## 3. judgements
  if(length(which(names(methods::slot(taxobj,"configuration"))=="configuration"))==0){
    sub_obj=lapply(methods::slot(taxobj, "data")[which(!names(methods::slot(taxobj,"data"))%in% c("configuration"))],select_temp)
  }else{
    sub_obj=lapply(methods::slot(taxobj, "data")[which(!names(methods::slot(taxobj,"data"))%in% c("configuration","configuration"))],select_temp)
    sub_obj$configuration=methods::slot(taxobj,"configuration")$configuration
  }

  ## 4.
  sub_obj$groupfile=methods::slot(taxobj,"groupfile") %>% .[select_num,]
  sub_obj$configuration=methods::slot(taxobj,"configuration")

  ## 5. Re-summary
  if(!is.null(taxnum)){
    if(length(taxnum)>nrow(methods::slot(taxobj,"data")$Base_percent)){warning("taxnum does not match Base table, please check!!!")}
    groupfiletemp=sub_obj$groupfile
    select_Base_pct=methods::slot(taxobj,"data")$Base_percent[taxnum,]
    select_Base_tax=methods::slot(taxobj,"data")$Base_taxonomy[taxnum,]
    sink("./sub_tax_summary_temp.txt")
    temptax_summary=tax_summary(groupfile =groupfiletemp,
                                inputtable =select_Base_pct[,-1],
                                reads = FALSE,
                                taxonomytable = select_Base_tax[,c(2,3)],
                                into = sub_obj$configuration$into,
                                sep = sub_obj$configuration$sep,
                                outputtax = sub_obj$configuration$outputtax) %>%
      suppressMessages()
    sink()
    file.remove("./sub_tax_summary_temp.txt")
    methods::slot(temptax_summary,"data")$Base=temptax_summary$Base_percent %T>%
      {.[,-1]=sweep(.[,-1],colSums(sub_obj$Base[,-1])%>%as.numeric,"*",MARGIN=2)%>% round(.,0)}
    if("standard" %in% temptax_summary$configuration$outputtax){
      outputtax=c("Domain","Phylum","Class","Order","Family","Genus","Species")
    }else if("complete" %in% temptax_summary$configuration$outputtax){
      outputtax=c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    }else{outputtax=temptax_summary$configuration$outputtax}
    for(i in outputtax){
      add=eval(parse(text=paste0("temptax_summary","$",i,"_percent"))) %T>%
        {.[,-1]=sweep(.[,-1],colSums(sub_obj$Base[,-1])%>%as.numeric,"*",MARGIN=2)%>% round(.,0)}
      methods::slot(temptax_summary,"data")=c(methods::slot(temptax_summary,"data"),list(add))
      names(methods::slot(temptax_summary,"data"))[length(names( methods::slot(temptax_summary,"data")))]=i
    }
    methods::slot(temptax_summary,"configuration")=taxobj$configuration
    methods::slot(temptax_summary,"groupfile")=groupfiletemp
    sub_obj=methods::slot(temptax_summary,"data")[names(methods::slot(temptax_summary,"data")) %in% names(sub_obj)]
    sub_obj$configuration=taxobj$configuration
    sub_obj$groupfile=groupfiletemp
    rm(temptax_summary)
  }

  output_obj <- methods::new("LorMe",
                             groupfile = sub_obj$groupfile,
                             data      = sub_obj[setdiff(names(sub_obj), c("groupfile", "parameters","configuration"))],
                             configuration    = sub_obj$configuration)

  sink()%>% suppressWarnings()
  return(output_obj)
}
