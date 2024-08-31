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
#' @examples
#'   data("Three_group")
#'
#'   # Check meta file
#'   print(Three_group$Groupfile)
#'
#'   # Subsetting tax summary objects
#'
#'   # Select BF and OF groups
#'   sub_testtax_summary <- sub_tax_summary(Three_group, Group %in% c("BF", "OF"))
#'   print(sub_testtax_summary$Groupfile)
#'
#'   # Subsetting according to taxonomy
#'
#'   Proteo <- sub_tax_summary(
#'     Three_group,
#'     taxnum = which(Three_group$Base_taxonomy$Phylum == "p__Proteobacteria")
#'   )
#'   print(Proteo$Phylum_percent)  # Check phylum table
#'   print(Proteo$Genus_percent)   # Check genus table
sub_tax_summary=function(taxobj,...,specificnum=NULL,taxnum=NULL){
  if(is.null(specificnum)){
    sub_group=subset(taxobj$Groupfile,...)
    sub_rownames= rownames(sub_group)
    select_num=which(rownames(taxobj$Groupfile) %in% sub_rownames)
  }else{select_num=specificnum}
  select_temp=function(x){if(is.numeric(x[,-1]%>% as.matrix())==TRUE){x[,c(1,select_num+1)]}else(x=x)}
  if(length(which(names(taxobj)=="configuration"))==0){
    sub_obj=lapply(taxobj[which(!names(taxobj)%in% c("parameters"))],select_temp)
  }else{
    sub_obj=lapply(taxobj[which(!names(taxobj)%in% c("parameters","configuration"))],select_temp)
    sub_obj$configuration=taxobj$configuration
  }
  sub_obj$Groupfile=taxobj$Groupfile %>% .[select_num,]
  sub_obj$parameters=taxobj$parameters
  if(is.null(taxnum)==FALSE){
    if(length(taxnum)>nrow(sub_obj$Base_percent)){warning("taxnum does not match Base table, please check!!!")}
    Groupfiletemp=sub_obj$Groupfile
    select_Base_pct=sub_obj$Base_percent[taxnum,]
    select_Base_tax=sub_obj$Base_taxonomy[taxnum,]
    sink("./sub_tax_summary_temp.txt")
    temptax_summary=tax_summary(groupfile =Groupfiletemp,inputtable =select_Base_pct[,-1],reads = FALSE,taxonomytable = select_Base_tax[,c(2,3)],into = sub_obj$parameters$into,sep = sub_obj$parameters$sep,outputtax = sub_obj$parameters$outputtax)%>% suppressMessages()
    sink()
    file.remove("./sub_tax_summary_temp.txt")
    temptax_summary$Base=temptax_summary$Base_percent %T>% {.[,-1]=sweep(.[,-1],colSums(sub_obj$Base[,-1])%>%as.numeric,"*",MARGIN=2)%>% round(.,0)}
    if("standard" %in% temptax_summary$parameters$outputtax){
      outputtax=c("Domain","Phylum","Class","Order","Family","Genus","Species")
    }else if("complete" %in% temptax_summary$parameters$outputtax){
      outputtax=c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    }else{outputtax=temptax_summary$parameters$outputtax}
    for(i in outputtax){
      add=eval(parse(text=paste0("temptax_summary","$",i,"_percent"))) %T>%
        {.[,-1]=sweep(.[,-1],colSums(sub_obj$Base[,-1])%>%as.numeric,"*",MARGIN=2)%>% round(.,0)}
      temptax_summary=c(temptax_summary,list(add))
      names(temptax_summary)[length(names(temptax_summary))]=i
    }
    temptax_summary$configuration=taxobj$configuration
    temptax_summary$Groupfile=Groupfiletemp
    sub_obj=temptax_summary[names(temptax_summary) %in% names(sub_obj)]
    sub_obj$parameters=taxobj$parameters
    rm(temptax_summary)
  }
  sink()%>% suppressWarnings()
  return(sub_obj)
}
