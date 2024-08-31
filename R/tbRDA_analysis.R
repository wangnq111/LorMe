###version 1.2.1###
###Auther Wang Ningqi###
#' tbRDA analysis
#' @description RDA analysis including co-linearity diagnostics and necessary statistics.
#' @param otudata Feature table of all numeric variable, with annotation in row names
#' @param envdata Environmental factor of all numeric variable,with sample-ID in row names and environmental factor in column names
#' @param collinearity If done collinearity diagnostics. Default,TRUE.
#' @param perm.test Logical. If conduct permutation test. Default:TRUE.
#'
#' @note
#' 1. When Axis length in first axis more than 4, you should choose CCA instead of RDA.
#' @author  Wang Ningqi <2434066068@qq.com>
#' @return  Three permutation test result print ,one preview plot ,a RDA object(default name:otu.tab.1) and a summary of RDA object
#' @export
#'
#' @importFrom stats na.omit
#' @importFrom vegan decostand decorana rda vif.cca anova.cca RsquareAdj
#'
#' @examples
#'   ### Data preparation ###
#'   library(vegan)
#'   data(varechem)
#'   head(varechem)
#'   data(testotu)
#'   require(tidyr); require(magrittr)  ## Or use pipe command in "dplyr"
#'
#'   sep_testotu <- Filter_function(
#'     input = testotu,
#'     threshold = 0.0001,
#'     format = 1
#'   ) %>%
#'   separate(
#'     ., col = taxonomy,
#'     into = c("Domain", "Phylum", "Order", "Family", "Class", "Genus", "Species"),
#'     sep = ";"
#'   )
#'
#'   top10phylum <- aggregate(
#'     sep_testotu[, 2:21],
#'     by = list(sep_testotu$Phylum),
#'     FUN = sum
#'   ) %>%
#'   Top_taxa(
#'     input = .,
#'     n = 10,
#'     inputformat = 2,
#'     outformat = 1
#'   )
#'   rownames(top10phylum) <- top10phylum[, 1]
#'   top10phylum <- top10phylum[, -1]
#'
#'   group <- data.frame(
#'     group = c(rep("a", 10), rep("b", 10)),
#'     factor1 = rnorm(10),
#'     factor2 = rnorm(mean = 100, 10)
#'   )
#'
#'   ### RDA analysis ###
#'   set.seed(999)
#'   RDAresult <- tbRDA_analysis(
#'     top10phylum,
#'     varechem[1:20, ],
#'     TRUE
#'   )
#'
#'   # Environmental statistics
#'   print(RDAresult$envstat)
#'
#'   # Visualization using ggplot
#'   rda_object <- RDAresult$rda_object
#'   rda_summary <- RDAresult$rdasummary
#'   rda_env <- as.data.frame(rda_summary$biplot)
#'   rda_sample <- as.data.frame(rda_summary$sites)
#'   rda_otu <- as.data.frame(rda_summary$species)
#'
#'   xlab <- paste0("RDA1:", round(RDAresult$rdasummary$concont$importance[2, 1], 4) * 100, "%")
#'   ylab <- paste0("RDA2:", round(RDAresult$rdasummary$concont$importance[2, 2], 4) * 100, "%")
#'
#'   library(ggrepel)
#'   # Create a sample RDA plot
#'   RDAplot <- ggplot(data = rda_sample, aes(RDA1, RDA2)) +
#'     geom_point(aes(color = group$group), size = 2) +
#'     geom_point(data = rda_otu, pch = "+", color = "orange", size = 4) +
#'     geom_hline(yintercept = 0) +
#'     geom_vline(xintercept = 0) +
#'     geom_segment(data = rda_env, aes(x = 0, y = 0, xend = RDA1 * 0.8, yend = RDA2 * 0.8),
#'                  arrow = arrow(angle = 22.5, length = unit(0.35, "cm")),
#'                  linetype = 1, size = 0.6, colour = "red") +
#'     geom_text_repel(color = "red", data = rda_env,
#'                     aes(RDA1, RDA2, label = row.names(rda_env))) +
#'     labs(x = xlab, y = ylab, color = "Treatment",
#'          title = paste0("p = ", anova.cca(rda_object)["Model", "Pr(>F)"])) +
#'     stat_ellipse(aes(color = group$group), level = 0.95) +
#'     geom_text_repel(size = 3, color = "orange",
#'                     data = subset(rda_otu, RDA1 > 0.1 | RDA1 < (-0.1)),
#'                     aes(RDA1, RDA2, label = rownames(subset(rda_otu, RDA1>0.1|RDA1<(-0.1))))) +
#'     theme_zg()
#'
#'   # Print the RDA plot
#'   print(RDAplot)

tbRDA_analysis<-function(otudata,envdata,collinearity,perm.test=TRUE){
  env <- stats::na.omit(log1p(envdata))
  otu.hell <- decostand(t(otudata), "hellinger")
  if(max(decorana(otu.hell)$rproj[,1])>4){warning("Axis length in first axis more than 4,please use CCA analysis!!!")}
  otu.tab.1<- vegan::rda(otu.hell ~ ., env)
  Diagnostic<-vif.cca(otu.tab.1)
  Diagnostic=Diagnostic[!is.na(Diagnostic)]
  if(collinearity==TRUE){
    orgi<-length(Diagnostic)
    for (i in 1:orgi){             ##delate max until all value less than 10###
      if(max(Diagnostic)>10){
        env=env[,-(which(Diagnostic==max(Diagnostic)))]
        otu.tab.1<- rda(otu.hell ~ ., env)
        Diagnostic<-vif.cca(otu.tab.1)
        Diagnostic=Diagnostic[!is.na(Diagnostic)]
        i=i-1}}}else if (collinearity==FALSE){} else{stop("Please choose TRUE/FALSE on collinearity");return()}
  summaryrda<-summary(otu.tab.1)
  if(perm.test==TRUE){
    message("###Global permutation test###","\n")
    print(anova.cca(otu.tab.1))
    message("###permutation test by term###","\n")
    print(anova.cca(otu.tab.1, by = "term"))
    message("###permutation test by axis###","\n")
    print(anova.cca(otu.tab.1, by = "axis"))
    plot(otu.tab.1)
  }
  terms=anova.cca(otu.tab.1, by = "term")
  env.adj.r.squared=as.numeric()
  for (i in 1:ncol(env)){
    env.adj.r.squared=c(env.adj.r.squared,RsquareAdj(rda(otu.hell ~ env[,i], env))$r.squared)
  }
  env.adj.r.squared=c(env.adj.r.squared,RsquareAdj(rda(otu.hell , env))$r.squared)
  envstat=data.frame(terms,env.adj.r.squared,ID=rownames(terms))
  rownames(envstat)[nrow(envstat)]<- envstat[nrow(envstat),6]<-"Total"
  envstat[nrow(envstat),1:4]=NA
  envstat<-envstat
  outlist=c(list(summaryrda),list(envstat),list(otu.tab.1))
  names(outlist)=c("rdasummary","factor_statistics","rda_object")
  return(outlist)
}
