
hub_detection_plot <- function(data = NULL, plot = TRUE){

  if(is.null(data)) stop("hub_detection_mdsd object is missing with no default.")

  p <- nrow(data$Degree_lambda)

  degrees <- c(t(data$Degree_lambda))

  lambdas <- rep(data$lambda, times = p)

  node_names <- data$node_names

  all_nodes <- rep(node_names, each = ncol(data$Degree))

  degree_data <- data.frame(node = all_nodes,
                        Degree = degrees,
                        lambda = lambdas)

  mdsd_data <- data.frame(node = node_names,
                         MDSD = data$MDSD,
                         lambda = data$lambda)

  p_degree <- ggplot(degree_data,
              aes(x = lambda, y = Degree, group = node)
  ) +
    geom_line() +
    ylab("Node degree")

  p_mdsd <- ggplot(mdsd_data, aes(x = node, y = MDSD)) +
    geom_segment(
      aes(x = node, xend = node, y = 0, yend = MDSD)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ylab("MDSD value")

  if(!plot) p_degree <- NULL
  if(!plot) p_mdsd <- NULL

  results <- list(lambda = lambda,
                 degree_data = degree_data,
                 mdsd_data = mdsd_data,
                 degree_plot = p_degree,
                 MDSD_plot = p_mdsd
                 )

  return(results)

}
