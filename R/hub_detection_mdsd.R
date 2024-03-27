#' @title hub_detection_mdsd
#'
#' @description Computes the MDSD distance from the solution path information.
#'
#' @param sol_path a list containing an array \code{wi} of spare matrices and the corresponding tuning parameter values \code{rho_list}.
#' @param node_names node names. If not provided, nodes are numbered from 1 to \code{ncol(sol_path$wi)}. Default \code{NULL}.
#' @param gamma the multiplier of the 'outlier' cut-off value. Default value \code{3}.
#' @param skew_thr the skewness threshold. Only a subset of models with estimated degree distribution skewness greater than \code{skew_thr} are inspected. Default value \code{0.5}.
#'
#' @return A list of potential hub nodes, MDSD distance values, and node degrees.
#'
#' @examples
#'
#' library(huge)
#' library(igraph)
#' library(MASS)
#'
#' n <- 100
#' p <- 50
#'
#' Data <- huge::huge.generator(n = n, d = p, graph = "hub")
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
#' G_true <- igraph::graph_from_adjacency_matrix(A_true, mode = "undirected", diag = FALSE)
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
#' X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma)
#'
#' X <- scale(X)
#'
#' nlambda <- 70
#'
#' cor_path <- cor_screening(X = X, nlambda = nlambda)
#'
#' res <- hub_detection_mdsd(sol_path = cor_path)
#'
#' @export
#' @importFrom huge huge.generator
#' @importFrom igraph graph_from_adjacency_matrix degree
#' @importFrom MASS mvrnorm
#' @importFrom e1071 skewness
#' @importFrom rlang is_empty

hub_detection_mdsd <- function(sol_path = NULL, node_names = NULL, gamma = NULL, skew_thr = NULL){

  if(is.null(gamma)) gamma <- 3

  if(is.null(skew_thr)) skew_thr <- 0.5

  p <- ncol(sol_path$wi)

  if(is.null(node_names)) node_names <- 1:p

  lambda <- sol_path$rholist

  nlambda <- length(lambda)

  Degree_lambda <- matrix(0, nrow = p, ncol = nlambda)

  colnames(Degree_lambda) <- lambda
  rownames(Degree_lambda) <- node_names

  for(j in 1:nlambda){

    A_temp <- sol_path$wi[,, j]

    A_temp[A_temp != 0] <- 1

    diag(A_temp) <- 0

    d_temp <- colSums(A_temp)

    names(d_temp) <- node_names

    Degree_lambda[node_names, j] <- d_temp

  }

  rm(A_temp)

  MDSD <- rep(0, p)

  for(j in 1:p){

    MDSD[j] <- mean((Degree_lambda[rep(j, p - 1), ] - Degree_lambda[-j, ])^2)

  }

  names(MDSD) <- rownames(Degree_lambda)

  d_skewness = apply(Degree_lambda, 2, e1071::skewness)

  d_skewness[is.na(d_skewness)] = 0

  burn = which(d_skewness <= skew_thr)

  if(!rlang::is_empty(burn)){

    MDSD_burn = rep(0, p)

    for(j in 1:p){

      MDSD_burn[j] = mean((Degree_lambda[rep(j, p - 1), -burn] - Degree_lambda[-j, -burn])^2)

    }

  }else{

    MDSD_burn = MDSD

  }

  names(MDSD_burn) <- rownames(Degree_lambda)

  hub_nodes_MDSD <- node_names[MDSD > gamma*mean(MDSD)]

  hub_nodes_MDSD_burn <- node_names[MDSD_burn > gamma*mean(MDSD_burn)]

  return(list(hub_nodes_MDSD = hub_nodes_MDSD,
              hub_nodes_MDSD_burn = hub_nodes_MDSD_burn,
              MDSD = MDSD,
              MDSD_burn = MDSD_burn,
              Degree = Degree_lambda,
              lambda = lambda,
              node_names = node_names))

}
