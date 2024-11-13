
library(dplyr)
library(tidyverse)
library(igraph)
library(ggfortify)
library(FactoMineR) 
library(factoextra)  
library(centiserve) # for Katz centrality


# COMPUTE CENTRALITY SNAILS

## read bipartite matrix

bipartite_matrix <- readRDS(here::here("output/temporal_output/bipartite_matrix.rds"))

############## COMPUTE CENTRALITY

bipartite_matrix <- as.matrix(bipartite_matrix)

# Define types manually based on the bipartite structure of the matrix
num_snails <- nrow(bipartite_matrix)
num_otus <- ncol(bipartite_matrix)

# Create the bipartite graph
g <- graph_from_biadjacency_matrix(bipartite_matrix, directed = FALSE)

# Set types: FALSE for snail nodes (rows), TRUE for OTU nodes (columns)
V(g)$type <- c(rep(FALSE, num_snails), rep(TRUE, num_otus))

# Filter snail nodes to verify correct assignment
snail_nodes <- V(g)[V(g)$type == FALSE]

degree_centrality <- degree(g, v = snail_nodes)
closeness_centrality <- closeness(g, v = snail_nodes)
katz_centrality <- katzcent(g, alpha = 0.05, vids = snail_nodes)

# Combine the centralities into a data frame and standardize them (z-score)
centralities_df <- data.frame(
  species = names(degree_centrality),
  degree_z = scale(degree_centrality),
  closeness_z = scale(closeness_centrality),
  katz_z = scale(katz_centrality)
)

# Perform PCA on the standardized centrality measures
pca_result <- PCA(centralities_df[, -1], scale.unit = TRUE, graph = FALSE)
pca_centrality <- as.data.frame(pca_result$ind$coord)[, 1]  # First principal component

# Add PCA centrality scores to the centralities data frame
centralities_df <- centralities_df %>%
  mutate(PCA_centrality = pca_centrality)

print(centralities_df)

#saveRDS(centralities_df, here::here("output/tables/centralities_df.rds"))




############## PLOT

pca_model <- prcomp(centralities_df[,-1], center = TRUE, scale. = TRUE)

# Summary to get % of variance explained by each component
pca_summary <- summary(pca_model)
explained_variance <- pca_summary$importance[2, ] * 100  # % variance for each PC

# Plot PCA
autoplot(pca_model, data = centralities_df, 
         loadings = TRUE, loadings.label = TRUE,
         loadings.label.size = 3, color = 'black') +
  labs(x = paste("PC1 (", round(explained_variance[1], 1), "%)", sep=""),
       y = paste("PC2 (", round(explained_variance[2], 1), "%)", sep="")) +
  theme_minimal()
