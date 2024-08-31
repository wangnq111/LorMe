#' @_PACKAGE
#' @title LorMe package: Lightening One-Code Resolving Microbial Ecology Program
#' @name LorMe
#' @description LorMe package summarizes a series of functions normally used in microbiome analysis
#' analysis.
#' @details
#' #Basic functions####
#'
#' \code{\link{auto_signif_test}} Automatically conduct significance testing
#'
#' \code{\link{compare_plot}} Comparison plot generator
#'
#' \code{\link{Filter_function}} Filter OTU/ASV/metagenomic profile/gene profile by threshold
#'
#' \code{\link{tax_summary}} Encapsulate meta file, feature tables and taxonomy annotation into tax summary object
#'
#' \code{\link{sub_tax_summary}} subsets tax summary objects according to meta file
#'
#' \code{\link{combine_and_translate}} Combine feature table with meta file and transform into a recognizable data frame for visualization.
#'
#' \code{\link{color_scheme}} generate color scheme from nine color scheme database and expand into colorRamp
#'
#' \code{\link{theme_zg}} A classic theme for ggplot.
#'
#' #Community features####
#'
#' \code{\link{Alpha_diversity_calculator}} Calculator for alpha diversity of each sample.
#'
#' \code{\link{Dimension_reduction}} Dimension reduction analysis including PCA,PCOA and NMDS
#'
#' \code{\link{structure_plot}} A fast view of microbial structure with PCA plot,PCOA plot and NMDS plot.
#'
#' \code{\link{Top_taxa}} Calculate most abundant taxon
#'
#' \code{\link{community_plot}} A fast view of microbial community with bar plot,alluvial plot and area plot.
#'
#' #Differential analysis####
#'
#' \code{\link{Deseq_analysis}} Performs a differential expression analysis
#'
#' \code{\link{indicator_analysis}} Performs the indicator analysis based on taxonomic summary object
#'
#' \code{\link{differential_bar}} Generate Differential Bar Plot and errorbar plot
#'
#' \code{\link{volcano_plot}} Generate volcano plot base on Deseq_analysis or indicator_analysis results
#'
#' \code{\link{manhattan}} Generate Manhattan Plot  base on Deseq_analysis or indicator_analysis results
#'
#' #Network analysis####
#'
#' \code{\link{network_analysis}} A convenient and fast network analysis function, with output results suitable for cytoscape and gephi
#'
#' \code{\link{network_withdiff}} Meta network analysis integrating differential taxon into a network analysis
#'
#' \code{\link{network_visual}} Visualizes a network based on network object from \code{\link{network_analysis}}
#'
#' \code{\link{network_visual_re}}  Re-visualize or adjust network plot from \code{\link{network_visual}} or \code{\link{network_withdiff}}
#'
#' \code{\link{Module_composition}} Pie chart for network module composition
#'
#' \code{\link{Module_abundance}} Calculate network module abundance for each sample
#'
#' \code{\link{nc}} Calculate network Natural Connectivity
#'
#' \code{\link{NC_remove}} Conduct natural connectivity analysis
#'
#' #Correlation analysis####
#'
#' \code{\link{circulation_lm}} Quick test using circulation to fit linear models between one dependent variable and series of independent variable
#'
#' \code{\link{tbRDA_analysis}} RDA analysis including co-linearity diagnostics and necessary statistics.
#'
#' @author Wang Ningqi
NULL
