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
  group_id <- gsub(".+?([A-Za-z0-9]+) collagen I_normalized.csv", "\\1-COLI", file_path)
  group_id <- gsub(".+?([A-Za-z0-9]+) fibronectin_normalized.csv", "\\1-FN1", group_id)
  group_id <- gsub(".+?([A-Za-z0-9]+) plastic_normalized.csv", "\\1-PLSTC", group_id)
  
  return(list(data = tsne_data, group_id = group_id))
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
tsne_data <- do.call(rbind, lapply(tsne_data_list, function(x) {x$data}))

# Combine group IDs into a single vector
group_ids <- unlist(lapply(tsne_data_list, function(x) {x$group_id}))

# Generate labels vector
labels <- unique(group_ids)
labels <- gsub("^sgCntrl", "CTRL", labels)
labels <- gsub("sgRAI14\\(3\\)-6", "RAI14(3)-6", labels)
labels <- gsub("sgRAI14\\(3\\)", "RAI14(3)", labels)
labels <- gsub("collagen I", "COLI", labels)
labels <- gsub("fibronectin", "FN1", labels)
labels <- gsub("plastic", "PLSTC", labels)

# Plot t-SNE visualization
ggplot(data = tsne_data, aes(x = V1, y = V2, color = factor(group_id, levels = unique(group_ids)))) +
  geom_point() +
  ggtitle("t-SNE Visualization of Pancreatic Cancer Cells") +
  theme_minimal() +
  scale_color_discrete(labels = labels)

