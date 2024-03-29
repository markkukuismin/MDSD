% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hub_detection_plot.R
\name{hub_detection_plot}
\alias{hub_detection_plot}
\title{hub_detection_plot}
\usage{
hub_detection_plot(data = NULL, plot = TRUE)
}
\arguments{
\item{data}{a \code{hub_detection_mdsd} object.}

\item{plot}{plot the MDSD plot and the node degrees plotted against \code{lambda}. Default \code{TRUE}.}
}
\value{
A list of \code{lambda}, node degree data, MDSD values, degree plot, and MDSD plot.
}
\description{
Makes plots of the MDSD distance and node degrees.
}
\examples{

library(huge)
library(igraph)
library(MASS)
library(ggplot2)

n <- 100
p <- 50

Data <- huge::huge.generator(n = n, d = p, graph = "hub")

true_nmb_hubs <- ifelse(p > 40, ceiling(p/20), 2)

Sigma <- Data$sigma

A_true <- Data$theta

A_true <- as.matrix(A_true)

true_hubs <- colSums(A_true)

true_hubs <- order(-true_hubs)[1:true_nmb_hubs]

A_true[A_true != 0] <- 1

diag(A_true) <- 0

G_true <- igraph::graph_from_adjacency_matrix(A_true, mode = "undirected", diag = FALSE)

node_names <- as.character(1:p)

vertex_attr(G_true) <- list(name = node_names)

d_true <- igraph::degree(G_true)

degree_true <- colSums(A_true)

names(degree_true) <- 1:p

X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma)

X <- scale(X)

nlambda <- 70

cor_path <- cor_screening(X = X, nlambda = nlambda)

res <- hub_detection_mdsd(sol_path = cor_path)

plot_res <- hub_detection_plot(data = res, plot = FALSE)
}
