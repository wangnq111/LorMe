#' @title One-stop microbial analysis pipeline for LorMe objects
#'
#' @description
#' `LorMe_pipeline()` performs a complete microbial ecology analysis workflow
#' for a configured **LorMe** object, including:
#' community profiling, differential analysis, indicator species analysis,
#' treatment-specific sub-network construction, and global network construction.
#'
#' The function is modular and allows users to execute only selected steps
#' through the `step` argument.
#'
#' @param taxobj A configured **LorMe** object created by
#'   \code{\link{object_config}}.
#'   The object must contain preprocessed taxonomy tables, metadata,
#'   and analysis configurations.
#'
#' @param step Character vector specifying which analysis modules to run.
#'   Must be one or more of:
#'   \itemize{
#'     \item `"all"` – run the entire pipeline (default)
#'     \item `"profile"` – alpha/beta diversity and composition analysis
#'     \item `"diff"` – differential abundance (DESeq2, differential barplot)
#'     \item `"sub_net"` – treatment-specific subnetworks
#'     \item `"all_net"` – combined/co-occurrence network across all samples
#'   }
#'
#' @details
#' The function automatically respects global options set via
#' \code{\link{LorMe_options}}.
#' These options control:
#'   * analysis taxonomic level
#'   * palettes and plotting parameters
#'   * DESeq2 parameters
#'   * network analysis thresholds
#'
#' Internal failures in any analysis sub-module do **not** stop the pipeline:
#' the corresponding output is returned as `NULL`, and the function prints
#' a summary of failed steps on completion.
#'
#' @return
#' A named list containing the results of all requested modules.
#' Depending on `step`, the list may include:
#' \describe{
#'   \item{\code{alpha_results}}{Alpha diversity tables and plots}
#'   \item{\code{beta_results}}{Beta diversity ordinations (e.g., PCoA)}
#'   \item{\code{composition_results}}{Community composition barplots}
#'   \item{\code{diffbar_result}}{Differential barplot results}
#'   \item{\code{Deseq_results}}{DESeq2 results for each comparison}
#'   \item{\code{Deseq_volcano}}{Volcano plots for differential analysis}
#'   \item{\code{Deseq_manhattan}}{Manhattan plots for differential features}
#'   \item{\code{indicator_results}}{Indicator species results}
#'   \item{\code{indic_volcano}}{Volcano plot of indicator species analysis}
#'   \item{\code{indic_manhattan}}{Manhattan plot of indicator analysis}
#'   \item{\code{sub_network_results}}{Treatment-specific subnetworks}
#'   \item{\code{combine_network_results}}{Global co-occurrence network}
#' }
#'
#' Failed modules return as `NULL`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## View current global analysis options
#' getOption("LorMe")
#'
#' ## Set analysis options
#' LorMe_options(
#'   global = list(Analysis_level = "Genus"),
#'   sub_net = list(threshold = 0.7),
#'   all_net = list(threshold = 0.7)
#' )
#'
#' ## Run pipeline (time-consuming)
#' Two_group_analysis <- LorMe_pipeline(Two_group)
#'
#' ## Access results:
#' # Alpha diversity
#' Two_group_analysis$alpha_results$plotlist$Plotobj_Shannon$Boxplot
#'
#' # Beta diversity
#' Two_group_analysis$beta_results$PCoA_Plot
#'
#' # Community composition
#' Two_group_analysis$composition_results$barplot
#'
#' # Differential analysis
#' Two_group_analysis$Deseq_volcano$FC_FDR
#' Two_group_analysis$Deseq_manhattan$manhattan
#'
#' # Differential barplot
#' library(patchwork)
#' Two_group_analysis$diffbar_result$Barplot |
#'   Two_group_analysis$diffbar_result$Differenceplot
#'
#' # Subnetworks
#' require(magrittr)
#' Two_group_analysis$sub_network_results$Treatment_sub_network %>%
#'   network_visual()
#'
#' # Combined network
#' Two_group_analysis$combine_network_results %>% network_visual()
#'
#' ## Reset to default options
#' LorMe_defaults()
#'
#' ## Example: three-group comparison with custom options
#' LorMe_options(
#'   global = list(
#'     Analysis_level = "Species",
#'     compare_list = c("CF_OF", "CF_BF")
#'   ),
#'   all_net = list(threshold = 0.95, method = "pearson")
#' )
#'
#' Three_group_analysis <- LorMe_pipeline(Three_group)
#' }
LorMe_pipeline=function(taxobj,step="all"){
  if(length(which(step %in%  c("all", "profile", "diff", "sub_net", "all_net")==TRUE))!=length(step)){
    stop("Invalid characters in parameter 'step', please choose among  c('all', 'profile', 'diff', 'sub_net', 'all_net')")
    return()
  }
  ##set saferun####
  safe_run <- function(func) {
    tryCatch(
      func,
      error = function(e) {
        message(conditionMessage(e))
        NULL
      }
    )
  }
  ##Configuration check####
  if(is.null(taxobj$configuration$treat_location)){
    stop("LorMe object not configured yet, please try 'object_config'")
  }
  opt <- getOption("LorMe")
  A_level=opt$global$Analysis_level
  treat_stat=methods::slot(taxobj,"groupfile")[,taxobj$configuration$treat_location] %>% as.character()
  uni_treat= unique(treat_stat)
  fail_list=as.character()
  all_results=list()
  ##community profile####
  if("all" %in% step | "profile" %in% step){
    message("#Community profile analysis")
    alpha_results=safe_run(Alpha_diversity_calculator(taxobj,A_level))
    if(is.null(alpha_results)){fail_list=c(fail_list,"Alpha diversity")}
    all_results=c(all_results,list(alpha_results))
    names(all_results)[length(all_results)]="alpha_results"
    beta_results=safe_run(structure_plot(taxobj,A_level,
                                         ptsize = opt$beta$ptsize,
                                         diagram = opt$beta$diagram,
                                         ellipse.level = opt$beta$ellipse.level,
                                         facet_row = opt$beta$facet_row))
    all_results=c(all_results,list(beta_results))
    if(is.null(beta_results)){fail_list=c(fail_list,"Beta diversity")}
    names(all_results)[length(all_results)]="beta_results"
    composition_results=safe_run(community_plot(taxobj,taxlevel =opt$comp$taxlevel,
                                                n = opt$comp$n,
                                                palette = opt$comp$palette,
                                                nrow = opt$comp$nrow,
                                                rmprefix = opt$comp$rmprefix ))
    all_results=c(all_results,list(composition_results))
    if(is.null(composition_results)){fail_list=c(fail_list,"Community composition")}
    names(all_results)[length(all_results)]="composition_results"
  }
  ##Differential analysis####
  if("all" %in% step | "diff" %in% step){
    message("#")
    message("#Differential analysis")

    if(opt$manh$taxlevel==opt$comp$taxlevel){
      man_palette=composition_results$filled_color
    }else{
      man_palette=opt$manh$palette
    }
    run_diff=TRUE
    if(!is.null(methods::slot(taxobj,"configuration")$facet_order)&is.null(opt$global$compare_list)){
      run_diff=FALSE
      warning("You set a facet parameters but not assigned pairwise comparision yet,'DESeq'and 'Differential bar' will not run!")
    }
    if(length(uni_treat)>3){
      if(is.null(opt$global$compare_list)){
        warning("You have more than Four treatments but not assigned pairwise comparision yet,'DESeq'and 'Differential bar' will not run!")
        run_diff=FALSE
      }
    } else if(length(uni_treat)==3){
      if(is.null(opt$global$compare_list)){
        opt$global$compare_list=c(
          paste(uni_treat[1],uni_treat[2],sep="_"),
          paste(uni_treat[2],uni_treat[3],sep="_"),
          paste(uni_treat[1],uni_treat[3],sep="_")
        )
      }
    }

    if(run_diff==FALSE){
      message("#")
      message("#Differential bar analysis and DESeq not run")
      diffbar_result=NULL
      fail_list=c(fail_list,"Differential bar")
      Deseq_results=NULL
      fail_list=c(fail_list,"DESeq")
    }else{
      if(length(uni_treat)==2){
        message("#")
        message("#Differential bar analysis start")
        diffbar_result=safe_run(differential_bar(taxobj,taxlevel =A_level,
                                                 comparison =NULL,
                                                 rel_threshold = opt$diff_bar$rel_threshold,
                                                 anno_row = opt$diff_bar$anno_row,
                                                 aes_col = opt$diff_bar$aes_col,
                                                 limit_num = opt$diff_bar$limit_num ))
        if(is.null(diffbar_result)){fail_list=c(fail_list,"Differential bar")}
        if(is.null(opt$deseq$control_name)){
          deseq_control=uni_treat[1]
        }else{
          deseq_control=opt$deseq$control_name
        }
        message("#")
        message("#DESeq start")
        Deseq_results=safe_run(Deseq_analysis(taxobj,A_level,
                                              cutoff = opt$deseq$cutoff,
                                              control_name =  deseq_control,
                                              paired =  opt$deseq$paired,
                                              subject =  opt$deseq$subject))%>% suppressMessages()
        if(is.null(diffbar_result)){fail_list=c(fail_list,"Differential bar")}
        if(is.null(Deseq_results)){
          message("Deseq returns NULL, please note")
          fail_list=c(fail_list,"DESeq")
          Deseq_manhattan=NULL
          Deseq_volcano=NULL
        }else{
          Deseq_volcano=safe_run(volcano_plot(Deseq_results,cutoff = opt$deseq$cutoff,aes_col = methods::slot(taxobj,"configuration")$treat_col))
          if(is.null(opt$manh$controlname)){
            manh_control=uni_treat[1]
          }else{
            manh_control=opt$manh$controlname
          }
          Deseq_manhattan=safe_run(manhattan(Deseq_results,
                                             taxlevel = opt$manh$taxlevel,
                                             control_name =manh_control,
                                             mode = opt$manh$mode,
                                             top_n = opt$manh$top_n,
                                             palette = man_palette,
                                             select_tax = opt$manh$select_tax,
                                             rmprefix = opt$manh$rmprefix))
        }
      }else{
        Deseq_results=diffbar_result=Deseq_volcano=Deseq_manhattan=list()
        for(i in opt$global$compare_list){
          comp=strsplit(i,"_")[[1]]
          con_name=comp[1]
          message("#")
          message("#Differential bar analysis for ",i," start")
          diffbar_result_sub=safe_run(differential_bar(taxobj,taxlevel =A_level,
                                                       comparison =comp,
                                                       rel_threshold = opt$diff_bar$rel_threshold,
                                                       anno_row = opt$diff_bar$anno_row,
                                                       aes_col = opt$diff_bar$aes_col,
                                                       limit_num = opt$diff_bar$limit_num ))
          if(is.null(diffbar_result_sub)){fail_list=c(fail_list,paste0(i,"-Differential bar"))}
          diffbar_result=c(diffbar_result,list(diffbar_result_sub))
          names(diffbar_result)[length(diffbar_result)]=paste0(i,"_diffbar_result")
          ###
          message("#")
          message("#DESeq for ",i," start")
          Deseq_results_sub=safe_run(Deseq_analysis(taxobj,A_level,
                                                    cutoff = opt$deseq$cutoff,
                                                    comparison = comp,
                                                    control_name =  con_name,
                                                    paired =  opt$deseq$paired,
                                                    subject =  opt$deseq$subject))%>% suppressMessages()
          if(is.null(Deseq_results_sub)){
            message("Deseq for pairwise comparision ",i," returns NULL, please note")
            fail_list=c(fail_list,paste0(i,"-DESeq"))
          }else{
            Deseq_results=c(Deseq_results,list(Deseq_results_sub))
            names(Deseq_results)[length(Deseq_results)]=paste0(i,"_Deseq_results")
            ###
            Deseq_volcano_sub=safe_run(volcano_plot(Deseq_results_sub,
                                                    cutoff = opt$deseq$cutoff,
                                                    aes_col = methods::slot(taxobj,"configuration")$treat_col))
            Deseq_volcano=c(Deseq_volcano,list(Deseq_volcano_sub))
            names(Deseq_volcano)[length(Deseq_volcano)]=paste0(i,"_Deseq_volcano")
            ###
            Deseq_manhattan_sub=safe_run(manhattan(Deseq_results_sub,
                                                   taxlevel = opt$manh$taxlevel,
                                                   control_name =con_name ,
                                                   mode = opt$manh$mode,
                                                   top_n = opt$manh$top_n,
                                                   palette = man_palette,
                                                   select_tax = opt$manh$select_tax,
                                                   rmprefix = opt$manh$rmprefix))
            Deseq_manhattan=c(Deseq_manhattan,list(Deseq_manhattan_sub))
            names(Deseq_manhattan)[length(Deseq_manhattan)]=paste0(i,"_Deseq_manhattan")
          }
        }
      }
    }
    all_results=c(all_results,list(diffbar_result))
    names(all_results)[length(all_results)]="diffbar_result"
    all_results=c(all_results,list(Deseq_results))
    names(all_results)[length(all_results)]="Deseq_results"
    all_results=c(all_results,list(Deseq_volcano))
    names(all_results)[length(all_results)]="Deseq_volcano"
    all_results=c(all_results,list(Deseq_manhattan))
    names(all_results)[length(all_results)]="Deseq_manhattan"
    #
    message("#")
    message("#Indicator analysis start")
    indicator_results=indicator_analysis(taxobj,A_level,func =opt$indic$func,reads = opt$indic$reads )
    if(is.null(indicator_results)){fail_list=c(fail_list,"Indicator")}
    all_results=c(all_results,list(indicator_results))
    names(all_results)[length(all_results)]="indicator_results"
    if(length(uni_treat)==2){
      indic_volcano=safe_run(volcano_plot(indicator_results,
                                          cutoff = opt$deseq$cutoff,
                                          aes_col = methods::slot(taxobj,"configuration")$treat_col))
      all_results=c(all_results,list(indic_volcano))
      names(all_results)[length(all_results)]="indic_volcano"
      ###
      if(is.null(opt$deseq$control_name)){
        indic_control=uni_treat[1]
      }else{
        indic_control=opt$deseq$control_name
      }
      indic_manhattan=safe_run(manhattan(indicator_results,
                                         taxlevel = opt$manh$taxlevel,
                                         control_name =indic_control ,
                                         mode = opt$manh$mode,
                                         top_n = opt$manh$top_n,
                                         palette = man_palette,
                                         select_tax = opt$manh$select_tax,
                                         rmprefix = opt$manh$rmprefix))
      all_results=c(all_results,list(indic_manhattan))
      names(all_results)[length(all_results)]="indic_manhattan"

    }
  }
  ##sub network####
  if("all" %in% step | "sub_net" %in% step){
    ###Pre-judgment####
    message("#")
    message("#Network analysis for each treatment")
    run_subnet=TRUE
    if(max(treat_stat)<=5){
      warning("Your replications are not more than 5, sub network for each treatment will not run!")
      run_subnet=FALSE
    }else if (max(treat_stat)>5&max(treat_stat)<=7 ){
      warning("Your replications are not more than 8, please pay attention to false positive in networks!")
    }
    ###
    if(run_subnet==FALSE){
      sub_network_results=NULL
    }else{
      sub_network_results=list()
      for(i in unique(treat_stat)){
        num=which(treat_stat==i)
        sub_obj=sub_tax_summary(taxobj,specificnum = num)
        if(as.numeric(table(treat_stat)[i])<=5){
          warning("Replications of treatment '",i,"' are not more than 5, sub network will not run!")
          fail_list=c(fail_list,paste0("Treatment ",i,"-sub network"))
        }else if(as.numeric(table(treat_stat)[i])>5&as.numeric(table(treat_stat)[i])<=7 ){
          warning("Replications of treatment '",i,"'  are not more than 8, please pay attention to false positive in network!")
        }else{
          if(is.null(opt$sub_net$n)){
            nrep=as.numeric(table(treat_stat)[i])
          }  else {
            nrep= opt$sub_net$n
          }
          if (nrep>as.numeric(table(treat_stat)[i])){
            nrep=as.numeric(table(treat_stat)[i])
          }
          message("")
          message("##Running network for treatment '",i,"'")
          sub_network_temp=safe_run(network_analysis(sub_obj,A_level,
                                                     reads = opt$sub_net$reads,
                                                     n =nrep,
                                                     threshold =opt$sub_net$threshold,
                                                     rel_threshold =opt$sub_net$rel_threshold,
                                                     method = opt$sub_net$method,
                                                     display = opt$sub_net$display))
          if(is.null(sub_network_temp)){fail_list=c(fail_list,paste0("Treatment ",i,"-sub network"))}
          sub_network_results=c(sub_network_results,list(sub_network_temp))
          if(is.null(sub_network_results[[1]])){
            warning("Return NULL in running sub network for treatment '",i,"'")
          }else{
            names(sub_network_results)[length(sub_network_results)]=paste0(i,"_sub_network")
          }
        }
      }
    }
    all_results=c(all_results,list(sub_network_results))
    names(all_results)[length(all_results)]="sub_network_results"
  }
  ##combined network####
  if("all" %in% step | "all_net" %in% step){
    message("#")
    message("#Combined network analysis")
    run_comb_net=TRUE
    if(length(treat_stat)<=5){
      warning("Your replications are not more than 5, sub network for each treatment will not run!")
      run_comb_net=FALSE
    }else if (length(treat_stat)>5&length(treat_stat)<=7 ){
      warning("Your replications are not more than 8, please pay attention to false positive in networks!")
    }
    if(run_comb_net==TRUE){
      if(is.null(opt$all_net$n)){
        nrep=round(0.5*length(treat_stat),0)
      }  else {
        nrep=opt$all_net$n
      }
      if (nrep>as.numeric(table(treat_stat)[i])){
        nrep=length(treat_stat)
      }
      combine_network_results=safe_run(network_analysis(taxobj ,A_level,
                                                        reads = opt$all_net$reads,
                                                        n =nrep,
                                                        threshold =opt$all_net$threshold,
                                                        rel_threshold =opt$all_net$rel_threshold,
                                                        method = opt$all_net$method,
                                                        display = opt$all_net$display))
      if(is.null(combine_network_results)){fail_list=c(fail_list,"Combine network")}
    }
    all_results=c(all_results,list(combine_network_results))
    names(all_results)[length(all_results)]="combine_network_results"
  }

  message("#All analysis done")
  if(length(fail_list)!=0){
    message("##Notion: Following were Failed analysis")
    for(fail in fail_list){
      message(paste0("### ",fail," analysis"))
    }
  }
  return(all_results)
}
