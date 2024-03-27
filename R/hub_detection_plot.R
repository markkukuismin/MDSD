
#' @title hub_detection_plot
#'
#' @description Makes plots of the MDSD distance and node degrees.
#'
#' @param data a \code{hub_detection_mdsd} object.
#' @param plot plot the MDSD plot and the node degrees plotted against \code{lambda}. Default \code{TRUE}.
#'
#' @return A list of \code{lambda}, node degree data, MDSD values, degree plot, and MDSD plot.
#'
#' @examples
#'
#' library(huge)
#' library(igraph)
#' library(MASS)
#' library(ggplot2)
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
#' plot_res <- hub_detection_plot(data = res, plot = FALSE)
#' @export
#' @importFrom huge huge.generator
#' @importFrom igraph graph_from_adjacency_matrix degree
#' @importFrom MASS mvrnorm
#' @importFrom ggplot2 ggplot

hub_detection_plot <- function(data = NULL, plot = TRUE){

  if(is.null(data)) stop("hub_detection_mdsd object is missing with no default.")

  lambda <- NULL
  Degree <- NULL
  node <- NULL
  MDSD <- NULL
  MDSD_burn <- NULL

  p <- nrow(data$Degree)

  degrees <- c(t(data$Degree))

  lambdas <- rep(data$lambda, times = p)

  node_names <- data$node_names

  all_nodes <- rep(node_names, each = ncol(data$Degree))

  degree_data <- data.frame(node = all_nodes,
                        Degree = degrees,
                        lambda = lambdas)

  mdsd_data <- data.frame(node = node_names,
                         MDSD = data$MDSD,
                         MDSD_burn = data$MDSD_burn)

  p_degree <- ggplot2::ggplot(degree_data,
                              ggplot2::aes(x = lambda, y = Degree, group = node)
  ) +
    ggplot2::geom_line() +
    ggplot2::ylab("Node degree")

  p_mdsd <- ggplot(mdsd_data, ggplot2::aes(x = node, y = MDSD)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = node, xend = node, y = 0, yend = MDSD)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ggplot2::ylab("MDSD value")

  p_mdsd_burn <- ggplot(mdsd_data, ggplot2::aes(x = node, y = MDSD_burn)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = node, xend = node, y = 0, yend = MDSD_burn)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ggplot2::ylab("MDSD value when skewness is considered")

  if(!plot) p_degree <- NULL
  if(!plot) p_mdsd <- NULL
  if(!plot) p_mdsd_burn <- NULL

  results <- list(lambda = data$lambda,
                 degree_data = degree_data,
                 mdsd_data = mdsd_data,
                 degree_plot = p_degree,
                 MDSD_plot = p_mdsd,
                 MDSD_burn_plot = p_mdsd_burn
                 )

  return(results)

}
