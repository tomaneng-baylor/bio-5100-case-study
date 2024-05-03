# Load required libraries
library(Rtsne)
library(ggplot2)
library(dplyr)

# Load and preprocess the data
file_names <- c("sgCntrl collagen I_normalized.csv",
                "sgCntrl fibronectin_normalized.csv",
                "sgCntrl plastic_normalized.csv",
                "sgRAI14(3)-6 collagen I_normalized.csv",
                "sgRAI14(3)-6 fibronectin_normalized.csv",
                "sgRAI14(3)-6 plastic_normalized.csv",
                "sgRAI14(3) collagen I_normalized.csv",
                "sgRAI14(3) fibronectin_normalized.csv",
                "sgRAI14(3) plastic_normalized.csv")

# Read all files and combine them into one data frame
data <- lapply(file_names, function(file) {
  temp <- read.csv(file, header = TRUE)
  temp$group <- gsub(".csv", "", file)
  temp
})

# Combine all data into one dataframe
data <- do.call(rbind, data)

# Perform t-SNE
tsne_result <- Rtsne(data[, -ncol(data)], dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

# Create a dataframe with t-SNE results
tsne_df <- as.data.frame(tsne_result$Y)
tsne_df$group <- data$group

# Abbreviate group names
tsne_df$group <- gsub("sgCntrl collagen I_normalized", "sgCntrl COL1", tsne_df$group)
tsne_df$group <- gsub("sgCntrl fibronectin_normalized", "sgCntrl FBN", tsne_df$group)
tsne_df$group <- gsub("sgCntrl plastic_normalized", "sgCntrl PL", tsne_df$group)
tsne_df$group <- gsub("sgRAI14(3)-6 collagen I_normalized", "sgRAI14 COL1", tsne_df$group)
tsne_df$group <- gsub("sgRAI14(3)-6 fibronectin_normalized", "sgRAI14 FBN", tsne_df$group)
tsne_df$group <- gsub("sgRAI14(3)-6 plastic_normalized", "sgRAI14 PL", tsne_df$group)
tsne_df$group <- gsub("sgRAI14(3) collagen I_normalized", "sgRAI14 COL1", tsne_df$group)
tsne_df$group <- gsub("sgRAI14(3) fibronectin_normalized", "sgRAI14 FBN", tsne_df$group)
tsne_df$group <- gsub("sgRAI14(3) plastic_normalized", "sgRAI14 PL", tsne_df$group)

# Plot t-SNE results
tsne_plot <- ggplot(tsne_df, aes(x = V1, y = V2, color = group)) +
  geom_point() +
  theme_minimal() +
  labs(color = "Group") +
  ggtitle("t-SNE Plot") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

# Save t-SNE plot
ggsave("tsne_plot.png", tsne_plot, width = 6, height = 6)

# Perform k-means clustering
set.seed(123)
kmeans_result <- kmeans(data[, -ncol(data)], centers = 5)

# Add cluster labels to t-SNE dataframe
tsne_df$cluster <- as.factor(kmeans_result$cluster)

# Compute proportion of sgCntrl cells in each cluster
cluster_props <- tsne_df %>%
  group_by(cluster, group) %>%
  summarise(n = n()) %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n))

# Rank clusters by sgCntrl content
cluster_rank <- cluster_props %>%
  filter(grepl("sgCntrl", group)) %>%
  arrange(desc(prop)) %>%
  select(cluster)

# Plot t-SNE with cluster labels
tsne_cluster_plot <- ggplot(tsne_df, aes(x = V1, y = V2, color = cluster)) +
  geom_point() +
  theme_minimal() +
  labs(color = "Cluster") +
  ggtitle("t-SNE Plot with Clusters") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

# Save t-SNE plot with clusters
ggsave("tsne_cluster_plot.png", tsne_cluster_plot, width = 6, height = 6)

# Print cluster ranking
cat("Cluster ranking based on sgCntrl content:\n")
print(cluster_rank)

