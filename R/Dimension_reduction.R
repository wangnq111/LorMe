###version1.0.2 ###
###author Wangningqi####
#' Dimension_reduction: PCA, PCOA, and NMDS Analysis
#'
#' @description Performs dimension reduction analysis using PCA, PCOA, or NMDS.
#' @param inputframe An OTU/gene/taxa table with all numeric variables and no NA/NAN/inf values.
#' @param group Group information with the sample order the same as in inputframe.
#' @param format The format of analysis: 1 for PCA, 2 for PCOA, 3 for NMDS.
#' @return A list containing data frames and other statistics for dimension reduction analysis.
#' @importFrom stats prcomp
#' @importFrom vegan metaMDS vegdist
#' @importFrom ape pcoa
#'
#' @export
#' @note Inputframe should be a numeric matrix without NA/NAN/inf values.
#'
#'   The row names of inputframe should be set as OTU/gene/taxa annotations for further analysis.
#'
#'   The results are combined into a list for output. Use `as.data.frame(result[[1]])` to extract the data frame, and `$result$` to extract other statistics. See examples for details.
#' @author Wang Ningqi <2434066068@qq.com>
#'
#' @examples
#' ### Data preparation ###
#' data(testotu)
#' rownames(testotu) <- testotu[, 1]
#' inputotu <- testotu[, -c(1, ncol(testotu))]
#' head(inputotu)
#'
#' groupinformation1 <- data.frame(
#'   group = c(rep("a", 10), rep("b", 10)),
#'   factor1 = rnorm(10),
#'   factor2 = rnorm(mean = 100, 10)
#' )
#'
#' ### PCA ###
#' PCAresult <- Dimension_reduction(inputotu, groupinformation1, 1)
#' PCAframe <- PCAresult$outframe   # Extract data for visualization
#' head(PCAresult$data.pca$rotation,5)  # OTU coordinates
#'
#' ### PCOA ###
#' PCOAresult <- Dimension_reduction(inputotu, groupinformation1, 2)
#' PCOAframe <- PCOAresult$outframe  # Extract data for visualization
#' head(PCOAresult$PCOA$values,2)  # Explanation of first two axis
#'
#' ### NMDS ###
#' NMDSresult <- Dimension_reduction(inputotu, groupinformation1, 3)
#' NMDSframe <- NMDSresult$outframe  # Extract data for visualization
#' # Here we got a warning of `stress is (nearly) zero: you may have insufficient data`,
#' # so make sure you have sufficient data for NMDS
#' print(NMDSresult$NMDSstat$stress)  # Extract stress of NMDS
Dimension_reduction<-function(inputframe,group,format){
  if(format==2){
    distance_matrix=vegdist(t(inputframe),method = 'bray')
    PCOA<- pcoa(as.matrix(distance_matrix), correction="none", rn=NULL)
    result<- PCOA$values[,"Relative_eig"]
    pc<-as.data.frame(PCOA$vectors)
    xlab<-paste("PCOA1(",as.numeric(sprintf("%.3f",result[1]))*100,"%)",sep="")
    ylab<-paste("PCOA2(",as.numeric(sprintf("%.3f",result[2]))*100,"%)",sep="")
    outframe=data.frame(x = pc$Axis.1,y=pc$Axis.2)
    rownames(outframe)<-colnames(inputframe)
    colnames(outframe)<-c(xlab,ylab)
    outframe=cbind(outframe,group)
    outlist=c(list(outframe),list(PCOA))
    names(outlist)=c("outframe","PCOA")
    return(outlist)
  }else if(format==1){
    inputframe<- inputframe[which(rowSums(inputframe) > 0),]
    data.pca<- prcomp(t(inputframe), scale. = TRUE)
    summary_pca<-summary(data.pca)
    importance.pca=   summary_pca$importance
    rotation.pca<- summary_pca$rotation
    xy.pca<-summary_pca$x %>% as.data.frame()
    result.pca<-as.data.frame(rbind(importance.pca,rotation.pca,xy.pca))
    outframe<-cbind(xy.pca[,1:2],group)
    outlist<-c(list(outframe), list(data.pca))
    names(outlist)=c("outframe","data.pca")
    return(outlist)
  }else if(format==3){
    NMDSstat<-metaMDS(t(inputframe))
    outframe<-cbind(data.frame(MDS1 = NMDSstat$points[,1], MDS2 = NMDSstat$points[,2]),group)
    outlist<-c(list(outframe), list(NMDSstat))
    names(outlist)=c("outframe","NMDSstat")
    return(outlist)
  }
  else{warning("Please choose rigth format!!")}
}
