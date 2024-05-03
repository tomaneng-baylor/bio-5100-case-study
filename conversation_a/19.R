# Load required libraries
library(tidyverse)
library(Rtsne)

# Function to read and preprocess data
read_and_preprocess_data <- function(file_path) {
  # Read the data
  data <- read.csv(file_path, header = TRUE)
  
  # Extract cell IDs
  cell_ids <- data[["cell.id"]]
  
  # Remove cell IDs from the data
  data <- data[, -1]
  
  # Rename columns with non-alphanumeric characters
  colnames(data) <- make.names(colnames(data))
  
  # Run t-SNE
  set.seed(123)
  tsne_model <- Rtsne(data, dims = 2, perplexity = 30, verbose = TRUE)
  
  # Convert the result to a data frame
  tsne_data <- as.data.frame(tsne_model$Y)
  tsne_data$cell_id <- cell_ids
  
  # Extract group ID from file name
  group_id <- gsub(".+_([A-Za-z0-9]+)_normalized.csv", "\\1", file_path)
  
  # Abbreviate group IDs
  group_id <- gsub("collagen I", "COLI", group_id)
  group_id <- gsub("fibronectin", "FN1", group_id)
  
  # Add group ID to data
  tsne_data$group_id <- group_id
  
  return(tsne_data)
}

# List of file names
file_names <- c("sgCntrl collagen I_normalized.csv",
                "sgCntrl fibronectin_normalized.csv",
                "sgCntrl plastic_normalized.csv",
                "sgRAI14(3)-6 collagen I_normalized.csv",
                "sgRAI14(3)-6 fibronectin_normalized.csv",
                "sgRAI14(3)-6 plastic_normalized.csv",
                "sgRAI14(3) collagen I_normalized.csv",
                "sgRAI14(3) fibronectin_normalized.csv",
                "sgRAI14(3) plastic_normalized.csv")

# Directory containing data files
data_dir <- "data/"

# Read and preprocess each data file
tsne_data_list <- lapply(file.path(data_dir, file_names), read_and_preprocess_data)

# Combine data into a single data frame
tsne_data <- do.call(rbind, tsne_data_list)

# Perform k-means clustering on t-SNE-transformed data
set.seed(123)
k <- 3 # Number of clusters
tsne_data$cluster <- kmeans(tsne_data[,1:2], centers = k)$cluster

# Compute proportion of sgCntrl cells within each cluster
cluster_summary <- tsne_data %>%
  group_by(cluster, group_id) %>%
  summarise(count = n()) %>%
  spread(group_id, count, fill = 0) %>%
  ungroup() %>%
  mutate(total = rowSums(select(., -cluster)),
         sgCntrl_prop = ifelse(total == 0, 0, `sgCntrl` / total)) %>%
  arrange(desc(sgCntrl_prop))

# Plot t-SNE visualization with cluster colors
p <- ggplot(data = tsne_data, aes(x = V1, y = V2, color = as.factor(cluster))) +
  geom_point() +
  ggtitle("t-SNE Visualization of Pancreatic Cancer Cells") +
  theme_minimal() +
  scale_color_discrete(name = "Cluster")

# View the cluster summary
print(cluster_summary)

# Plot t-SNE visualization with cluster colors
p

