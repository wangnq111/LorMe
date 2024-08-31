#' Generate Differential Bar Plot and Error bar Plot
#'
#' @param taxobj Tax summary objects computed by \code{\link{tax_summary}}. Default: NULL.
#' @param taxlevel Taxonomy levels used for visualization. Must be one of
#'   c("Domain","Phylum","Class","Order","Family","Genus","Species","Base"). Default: NULL.
#' @param comparison A vector of conditions to compare. Default: NULL, all unique conditions are compared (only for Two groups).
#' @param rel_threshold Threshold filtering taxa for differential analysis.
#' @param anno_row Default: 'taxonomy'. Rownames for visualization. Options are 'taxonomy' for showing taxonomic information and 'ID' for showing taxonomic ID.
#' @param aes_col A named vector of colors to be used in the plots.
#' @param limit_num Numeric. The maximum number of significant results to display. Default: NULL, showing all differential taxa.
#'
#' @note The differential analysis is performed using two-sided Welch's t-test.
#'   The p-values are adjusted using the 'BH' (i.e., FDR) method.
#' @return A list containing the bar plot, source data for the bar plot, difference plot, and source data for the difference plot.
#' @export
#'
#' @importFrom stats t.test p.adjust
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar geom_point labs theme element_rect element_blank scale_x_discrete coord_flip margin annotate position_dodge
#' @importFrom scales percent
#'
#' @examples
#' {
#'   # Data preparation
#'   data("Two_group")
#'
#'   # Simple mode
#'   diff_results <- differential_bar(
#'     taxobj = Two_group,
#'     taxlevel = "Genus"
#'   )
#'   print(diff_results$Barplot)  # Print Barplot
#'   head(diff_results$Barplot_sourcedata)  # Show source data of barplot
#'   print(diff_results$Differenceplot)  # Print Differential errorbar plot
#'   head(diff_results$Differenceplot_sourcedata)  # Show source data of Differential errorbar plot
#'
#'   require(patchwork)
#'   diff_results$Barplot|diff_results$Differenceplot
#'   # Displaying ID
#'   diff_results <- differential_bar(
#'     taxobj = Two_group,
#'     taxlevel = "Base",
#'     anno_row = "ID"
#'   )
#'   print(diff_results$Barplot)
#'
#'   # Threshold adjustment
#'   diff_results <- differential_bar(
#'     taxobj = Two_group,
#'     taxlevel = "Base",
#'     rel_threshold = 0.001
#'   )
#'   print(diff_results$Barplot)
#'
#'   # Limit the displaying number
#'   diff_results <- differential_bar(
#'     taxobj = Two_group,
#'     taxlevel = "Base",
#'     rel_threshold = 0.001,
#'     limit_num = 10
#'   )
#'   print(diff_results$Barplot)
#'
#'   # For object with more than two groups
#'   # Data preparation
#'   data("Three_group")
#'
#'   # Specific comparison
#'   Three_group_col <- Three_group$configuration$treat_col
#'   diff_results <- differential_bar(
#'     taxobj = Three_group,
#'     taxlevel = "Genus",
#'     comparison = c("BF", "CF"),
#'     aes_col = Three_group_col
#'   )
#'   print(diff_results$Barplot)
#' }
differential_bar=function(taxobj,
                          taxlevel,
                          comparison = NULL,
                          rel_threshold=0.005,
                          anno_row="taxonomy",
                          aes_col=NULL,
                          limit_num=NULL){
  if (is.null(taxobj$configuration)) {
    stop("taxonomic summary object not configured yet, call '?object_config' for configuration!")
    return(NULL)
  }
  if (is.null(eval(parse(text = paste0("taxobj", "$", taxlevel))))) {
    warning("Illegal 'taxlevel'!")
    return(NULL)
  }
  groupfile = eval(parse(text = paste0("taxobj$Groupfile")))
  condition = groupfile[eval(parse(text = paste0("taxobj$configuration$treat_location")))]
  condition = condition[, 1]
  if (is.null(comparison)) {
    comparison = unique(condition)
  } else {
    if ((which(condition %in% comparison) %>% length()) ==
        0) {
      stop("comparision does not match, Please check 'comparision'")
    }
    taxobj = sub_tax_summary(taxobj = taxobj, specificnum = which(condition %in% comparison))
    groupfile= eval(parse(text = paste0("taxobj$Groupfile")))
    condition= groupfile[eval(parse(text = paste0("taxobj$configuration$treat_location")))]
    condition= condition[, 1]
  }
  if(is.null(aes_col)){
    aes_col=taxobj$configuration$treat_col
  }
  if (is.null(names(aes_col)[1])) {
    names(aes_col)=unique(taxobj$Groupfile[,taxobj$configuration$treat_location]) %>% sort()
    warning("Color names not assigned,use automatic assignment")
  }
  aes_col=aes_col[names(aes_col) %in% unique(condition)]
  taxonomy = eval(parse(text = paste0("taxobj", "$", taxlevel,
                                      "_taxonomy")))
  input0=eval(parse(text=paste0("taxobj","$",taxlevel,"_percent")))
  inputframe = data.frame(input0[, -1], row.names = taxonomy[,1])
  inputframe = inputframe[which(rowMeans(input0[,-1])>rel_threshold),]
  cal=function(x){t.test(x~condition) %>% return()}
  t_list=apply(inputframe,1,cal)
  i=1
  t_summary=data.frame()
  while(i<=length(t_list)){
    t_results=t_list[[i]]
    ID=rownames(inputframe)[i]
    p.value=t_results$p.value
    meanall=t_results$estimate%>% as.numeric()
    mean1=meanall[1]
    mean2=meanall[2]
    conf.low=t_results$conf.int[1]
    conf.high=t_results$conf.int[2]
    estimate=mean(c(conf.low,conf.high))
    t_summary=rbind(t_summary,data.frame(ID,p.value,estimate,conf.low,conf.high,mean1,mean2))
    i=i+1
  }
  colnames(t_summary)[6:7]=unique(taxobj$Groupfile$Group) %>% sort()
  t_summary$padj=p.adjust(t_summary$p.value,"BH")
  colnames(t_summary)[1]=colnames(taxonomy)[1]
  t_summary_anno=left_join(t_summary,taxonomy) %>% suppressMessages()
  t_summary_anno_sig=t_summary_anno[t_summary_anno$padj<0.05,]
  if(!is.null(limit_num)){
    if(nrow(t_summary_anno_sig)>limit_num){
      t_summary_anno_sig=t_summary_anno_sig[order(t_summary_anno_sig$p.value),]
      t_summary_anno_sig=t_summary_anno_sig[1:limit_num,]
    }
  }
  if(nrow(t_summary_anno_sig)==0){
    warning("No differences founded! Please check you data or decrease 'rel_threshold'")
    return(NULL)
  }
  if(taxlevel=="Base"){taxlevel="ID"}
  if(anno_row=="taxonomy"){
    anno_order=t_summary_anno_sig[,taxlevel]
    anno_order=anno_order[order(t_summary_anno_sig$estimate)]
    barframe=t_summary_anno_sig[,which(colnames(t_summary_anno_sig) %in% c(unique(condition),taxlevel))] %>%
      gather("Treatment","Mean",-taxlevel) %>% suppressMessages()
    barframe[,1]=factor(barframe[,1],levels=anno_order)
    t_summary_anno_sig[,taxlevel]=factor(t_summary_anno_sig[,taxlevel],levels=anno_order)
  }else if( anno_row=="ID"){
    anno_order=t_summary_anno_sig[,1]
    anno_order=anno_order[order(t_summary_anno_sig$estimate)]
    barframe=t_summary_anno_sig[,which(colnames(t_summary_anno_sig) %in% c(unique(condition),colnames(t_summary_anno_sig)[1]))]
    barframe= gather(barframe,"Treatment","Mean",-colnames(t_summary_anno_sig)[1]) %>% suppressMessages()
    barframe[,1]=factor(barframe[,1],levels=anno_order)
    t_summary_anno_sig[,1]=factor(t_summary_anno_sig[,1],levels=anno_order)
  }else{
    warning("Invalid 'anno_row',please select 'taxonomy' or 'ID'")
    return(NULL)
  }
  barplot=ggplot(barframe,aes(barframe[,1],barframe[,'Mean'],fill = barframe[,'Treatment'])) +
    scale_fill_manual(values=aes_col)+
    scale_x_discrete(limits = anno_order) +
    scale_y_continuous(expand = c(0,0),labels = scales::percent)+
    coord_flip() +#翻转XY轴
    labs(x="",y="Mean proportion")+
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=10,face = "bold"),
          axis.text=element_text(colour='black',size=8,face = "bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=10,face = "bold",colour = "black",
                                   margin = margin(r = 20)),
          legend.position = c(-1,-0.1),
          legend.direction = "horizontal",
          legend.key.width = unit(0.8,"cm"),
          legend.key.height = unit(0.5,"cm"))
  for (i in 1:(nrow(t_summary_anno_sig) - 1)){
    barplot <- barplot + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                                  fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }
  barplot = barplot+geom_bar(stat = "identity",position = "dodge",
                             width = 0.7,colour = "black")

  t_summary_anno_sig$tag[t_summary_anno_sig$estimate<0]=colnames(t_summary_anno_sig)[7]
  t_summary_anno_sig$tag[t_summary_anno_sig$estimate>0]=colnames(t_summary_anno_sig)[6]
  if(anno_row=="taxonomy"){
    diffplot=ggplot(t_summary_anno_sig,aes(t_summary_anno_sig[,taxlevel],
                                           t_summary_anno_sig[,'estimate'],
                                           fill = t_summary_anno_sig[,'tag']))
  }else{
    diffplot=ggplot(t_summary_anno_sig,aes(t_summary_anno_sig[,1],
                                           t_summary_anno_sig[,'estimate'],
                                           fill = t_summary_anno_sig[,'tag']))
  }

  diffplot=diffplot+
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=10,face = "bold"),
          axis.text=element_text(colour='black',size=8,face = "bold"),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 10,face = "bold",colour = "black",hjust = 0.5)) +
    scale_x_discrete(limits = anno_order) +
    scale_y_continuous(expand = c(0,0),labels = scales::percent)+
    coord_flip() +
    xlab("") +
    ylab("Difference in mean proportions") +
    labs(title="95% confidence intervals")

  for (i in 1:(nrow(t_summary_anno_sig) - 1)) {
    diffplot <- diffplot + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                                    fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }
  diffplot <- diffplot +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(0.8), width = 0.5, linewidth = 0.5) +#误差线
    geom_point(shape = 21,size = 3) +#散点图参数
    scale_fill_manual(values=aes_col) +#点颜色
    geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')#加虚线
  outlist=c(list(barplot),list(barframe),list(diffplot),list(t_summary_anno_sig))
  names(outlist)=c("Barplot","Barplot_sourcedata","Differenceplot","Differenceplot_sourcedata")
  return(outlist)
}
