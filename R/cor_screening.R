#' @title cor_screening
#'
#' @description Computes a path of sparse adjacency matrix estimates using correlation thresholding
#'
#' @param X a n by p data matrix containing iid observations. Each row is a new sample.
#' @param nlambda the number of tuning parameters. Default \code{NULL}.
#' @param lambda_min_ratio the largest sparsity level of the sparse adjacency matrix estimates. Default value \code{0.05}.
#'
#' @return A list containing array of sparse adjacency matrix estimates, and corresponding tuning parameter values.
#'
#' library(huge)
#'
#' n <- 100
#' p <- 50
#'
#' Data <- huge.generator(n = n, d = p, graph = "hub")
#'
#' true_nmb_hubs <- ifelse(p > 40, ceiling(p/20), 2)
#'
#' Sigma <- Data$sigma
#'
#' A_true <- Data$theta
#'
#' A_true <- as.matrix(A_true)
#'
#' true_hubs <- colSums(A_true)
#'
#' true_hubs <- order(-true_hubs)[1:true_nmb_hubs]
#'
#' A_true[A_true != 0] <- 1
#'
#' diag(A_true) <- 0
#'
#' G_true <- igraph::graph_from_adjacency_matrix(A_true, mode = "undirected", diag = F)
#'
#' node_names <- as.character(1:p)
#'
#' vertex_attr(G_true) <- list(name = node_names)
#'
#' d_true <- igraph::degree(G_true)
#'
#' degree_true <- colSums(A_true)
#'
#' names(degree_true) <- 1:p
#'
#' X = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma)
#'
#' X = scale(X)
#'
#' nlambda = 70
#'
#' cor_path = cor_screening(X = X, nlambda = nlambda)
#'
#' @export
#' @importFrom stats cor
#' @importFrom huge huge.generator
#' @importFrom MASS mvrnorm
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph degree

cor_screening <- function(X, nlambda = NULL, lambda_min_ratio = NULL){

  if(is.null(nlambda)) stop("nlambda is missing with no default")

  X <- scale(X)

  S <- cor(X)

  p <- ncol(S)

  if(is.null(lambda_min_ratio)) lambda_min_ratio <- 0.05

  cor_path <- list(wi = array(0, dim = c(p, p, nlambda)), rholist = rep(0, nlambda))

  lambda_max <- max(max(S - diag(p)), -min(S - diag(p)))

  lambda_min <- lambda_min_ratio*lambda_max

  lambda <- exp(seq(log(lambda_max), log(lambda_min), length = nlambda))

  S <- abs(S)

  for(i in 1:nlambda){

    cor_path$wi[,,i][S > lambda[i]] <- 1

  }

  cor_path$rholist <- lambda

  return(cor_path)

}
