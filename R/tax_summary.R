####Version1.0.0###
###Author Wangningqi####
#' @title Encapsulate meta file, feature tables and taxonomy annotation into tax summary object
#' @description The function packages meta file, feature tables and taxonomy annotation into tax summary object
#' @param groupfile A data frame containing treatment information
#' @param inputtable OTU/ASV/species data frame with all numeric. Samples ID should be in column names.
#' @param reads Logical.True for reads table and FALSE for percentage table. Default: TRUE
#' @param taxonomytable Taxonomy annotation data frame,with first column OTU/ASV/TAX number ID and second column taxonomy annotation. See details in example.
#' @param into Names of separated taxonomy to create as character vector. Must select from c("Domain","Phylum","Class","Order","Family","Genus","Species").
#'             Shortcut input:1)By default."standard":c("Domain","Phylum","Class","Order","Family","Genus","Species"). Used for standard taxonomy annotation to OTU/ASV table.
#'                            2)"complete":c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"). Used for complete taxonomy annotation to meta genomic table.
#' @param sep Separator of taxonomy table.Default: ";".
#' @param outputtax Names of output taxonomy level table. Default:c("Phylum","Genus"). Shortcut input is available with 'standard' and 'complete' same as above.
#'
#' @author Wang Ningqi <2434066068@qq.com>
#' @return One list containing taxonomy table data frame,containing reads and percentage table for each specified output. Full taxonomy annotation data frame is output in global environment.
#' @export
#'
#' @note
#' For taxonomy annotation with 'Kingdom' level, please set 'into' parameter as 'complete'!!!
#'
#' @import magrittr
#' @importFrom tidyr separate
#' @importFrom utils tail
#' @importFrom stats aggregate
#' @examples
#' {
#'   # Load data
#'   data(testotu)
#'
#'   # Create group information data frame
#'   groupinformation <- data.frame(
#'     group = c(rep("a", 10), rep("b", 10)),
#'     factor1 = rnorm(10),
#'     factor2 = rnorm(mean = 100, 10),
#'     subject = factor(c(1:10, 1:10))
#'   )
#'
#'   # Packaging data into a taxonomy summary object
#'   test_object <- tax_summary(
#'     groupfile = groupinformation,
#'     inputtable = testotu[, 2:21],
#'     reads = TRUE,
#'     taxonomytable = testotu[, c(1, 22)]
#'   )
#'
#'   # Check integrated object
#'   print(test_object)
#'
#'   # Extract genus relative abundance table
#'   test_Genus <- test_object$Genus_percent
#'   head(test_Genus)
#'
#'   # Check corresponding taxonomy information of genus table
#'   test_Genus_tax <- test_object$Genus_taxonomy
#'   head(test_Genus_tax)
#'
#'   # Summary base table into all taxonomy levels with standard output
#'   test_object <- tax_summary(
#'     groupfile = groupinformation,
#'     inputtable = testotu[, 2:21],
#'     reads = TRUE,
#'     taxonomytable = testotu[, c(1, 22)],
#'     outputtax = "standard"
#'   )
#'   head(test_object$Species_percent)  # View first 10 rows of species percentage
#'   head(test_object$Genus)  # View first 10 rows of genus table
#' }
tax_summary=function(groupfile,inputtable,reads=TRUE,taxonomytable,into="standard",sep=";",outputtax=c("Phylum","Genus")){
  parameters=list(into,sep,outputtax)
  if(into=="standard"){
    into=c("Domain","Phylum","Class","Order","Family","Genus","Species")
  }else if(into=="complete"){
    into=c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  }else{into=into}
  comp_tax=c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  if((which(into %in% comp_tax)%>% length())!=length(into)){warning("illegal input characters in parameter `into`")}
  if("standard" %in% outputtax){
    outputtax=c("Domain","Phylum","Class","Order","Family","Genus","Species")
  }else if("complete" %in% outputtax){
    outputtax=c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  }else{outputtax=outputtax}
  if((which(outputtax %in% comp_tax)%>% length())!=length(outputtax)){warning("Illegal input characters in parameter `outputtax`")}
  taxonomy=data.frame(taxonomytable,full_taxonomy=taxonomytable[,2]) %>%
    separate(.,col=3,into=into,sep=sep)
  taxonomy[is.na(taxonomy)]="Unassigned"
  taxonomy_update=taxonomy
  message("Taxonomy tidying...")
  for(j in into[!into %in% c("Domain","Kingdom")]){ #re_assign of taxonomy
    rp_obj=c("norank","uncultured","metagenome","Unassigned")
    for(rp in rp_obj){
      if(length(grep(rp,taxonomy_update[,j]))==0){
        taxonomy_update=taxonomy_update}else{
          select_frame=taxonomy_update[grep(rp,taxonomy_update[,j]),]
          fuc=function(x){
            temp_colname=colnames(select_frame)[3:which(colnames(select_frame)==j)]
            select_tempframe=x[3:which(colnames(select_frame)==j)]
            rp_num=c(grep(rp_obj[1],select_tempframe),grep(rp_obj[2],select_tempframe),grep(rp_obj[3],select_tempframe),grep(rp_obj[4],select_tempframe))
            prefix_temp=c("d__","k__","p__","c__","o__","f__","g__","s__")
            prefix_temp=prefix_temp[which(comp_tax==j)]
            rpname=paste0(select_tempframe[-rp_num] %>% as.character() %>% tail(1),";",prefix_temp,rp) %>% return()
          }
          taxonomy_update[grep(rp,taxonomy_update[,j]),j]=apply(select_frame,1,fuc)
        }
    }
  }
  message("Genarating tables...")
  output=list()
  output=c(output,list(groupfile))
  names(output)[1]="Groupfile"
  base_table=data.frame(ID=taxonomytable[,1],inputtable)
  base_table_pct=data.frame(ID=taxonomytable[,1],sweep(inputtable,colSums(inputtable),"/",MARGIN = 2))
  base_taxonomy=data.frame(ID=taxonomytable[,1],taxonomy_update)
  output_temp=list(base_table,base_table_pct,base_taxonomy)
  names(output_temp)[1]="Base"
  names(output_temp)[2]="Base_percent"
  names(output_temp)[3]="Base_taxonomy"
  output=c(output,output_temp)
  temp=function(x){
    as.character(x) %>%paste0(.,collapse =":")
  }
  for(j in outputtax){
    if(isTRUE(reads)){
      if(j=="Domain"){
        temp_taxonomy=data.frame(taxonomy=taxonomy_update[,3:which(colnames(taxonomy_update)==j)])
      }else{
        temp_taxonomy=data.frame(taxonomy=apply(taxonomy_update[,3:which(colnames(taxonomy_update)==j)],1,temp))
      }
      sumtable=aggregate(inputtable,by=list(temp_taxonomy[,1]),FUN=sum)
      colnames(sumtable)[1]=j
      sumtable_pct=data.frame(sumtable[,1],sweep(sumtable[,-1],colSums(sumtable[,-1]),"/",MARGIN = 2))
      colnames(sumtable_pct)[1]=j
      if(j=="Domain"){
        tax_table=data.frame(taxID=paste0(j,1:nrow(sumtable_pct)),tax=sumtable_pct[,1])
      }else{
        tax_table=data.frame(taxID=paste0(j,1:nrow(sumtable_pct)),tax=sumtable_pct[,1]) %T>%
          {colnames(.)=c(paste0(j,"ID"),"temp")} %>%
          separate(.,col=temp,into=colnames(taxonomy_update)[3:which(colnames(taxonomy_update)==j)],sep=":")
        sumtable[,1]=sumtable_pct[,1]=tax_table[,j]
      }

      output_temp=list(sumtable,sumtable_pct,tax_table)
      names(output_temp)[1]=j
      names(output_temp)[2]=paste0(j,"_percent")
      names(output_temp)[3]=paste0(j,"_taxonomy")
    }else{
      sumtable_pct=aggregate(inputtable,by=list(taxonomy_update[,j]),FUN=sum)
      colnames(sumtable_pct)[1]=j
      if(j=="Domain"){
        tax_table=data.frame(taxID=paste0(j,1:nrow(sumtable_pct)),tax=sumtable_pct[,1])
      }else{
        tax_table=data.frame(taxID=paste0(j,1:nrow(sumtable_pct)),tax=sumtable_pct[,1]) %T>%
          {colnames(.)=c(paste0(j,"ID"),j)} %>%
          left_join(.,unique(taxonomy_update[,3: which(colnames(taxonomy_update)==j)]))%>%suppressMessages()
        if(length(tax_table[,1])!=nrow(sumtable_pct)){
          repeat_name=which(table(tax_table[,1])!=1) %>% names()
          unique_name=which(table(tax_table[,1])==1) %>% names()
          keep_rownames=tax_table[which(tax_table[,1] %in% unique_name),] %>% rownames()
          for(i in repeat_name){
            keep_rownames=c(keep_rownames,tax_table[which(tax_table[,1] %in% i),] %>% rownames() %>% .[1])
          }
          tax_table=tax_table[which(rownames(tax_table) %in% keep_rownames),]
        }
        tax_table=tax_table[,c(1,3:ncol(tax_table),2)]}
      output_temp=list(sumtable_pct,tax_table)
      names(output_temp)[1]=paste0(j,"_percent")
      names(output_temp)[2]=paste0(j,"_taxonomy")
    }
    output=c(output,output_temp)
  }
  if(reads==FALSE){
    output$Base_percent=output$Base
    output$Base=NULL
  }
  output$parameters=parameters
  names(output$parameters)=c("into","sep","outputtax")
  message("Done")
  cat("\nContained elements:\n")
  names(output) %>% print()
  return(output)
}
