#' Calculate Network Natural Connectivity
#'
#' @param adj_matrix Adjacency data frame or matrix. Can be calculated from \code{\link{network_analysis}}
#'
#' @return Numeric value of natural connectivity
#' @export
#'
#' @examples
#' {
#'   ### Data preparation ###
#'   data("Two_group")
#'
#'   ### One input network analysis ###
#'   network_results <- network_analysis(
#'     taxobj = Two_group,
#'     taxlevel = "Base",
#'     reads = FALSE,
#'     n = 10,
#'     threshold = 0.6
#'   )
#'
#'   # Convert network results to a data frame for the adjacency matrix
#'   network_matrix <- as.data.frame(network_results[[3]])  # Complete adjacency matrix
#'
#'   # Check initial natural connectivity
#'   nc_initial <- nc(network_matrix)
#'   print(nc_initial)  # Print the initial natural connectivity
#' }
nc <- function(adj_matrix) {
  lambda0 <- as.matrix(adj_matrix)%>% abs() %>%
    eigen(only.values = TRUE)
  lambda<-lambda0$values %>% sort(decreasing = TRUE)
  lambda_sum <- 0
  for (i in 1:length(lambda)){
    lambda_sum = lambda_sum + exp(lambda[i])
  }
  log(lambda_sum/length(lambda)) %>% return()
}
