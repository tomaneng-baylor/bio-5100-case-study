# Load required libraries
library(tidyverse)
library(Rtsne)

# Read the data
data <- read.delim("your_data_file.txt", header = TRUE, row.names = 1)

# Run t-SNE
set.seed(123)
tsne_model <- Rtsne(data, dims = 2, perplexity = 30, verbose = TRUE)

# Convert the result to a data frame
tsne_data <- as.data.frame(tsne_model$Y)
tsne_data$cell_id <- rownames(data)

# Plot t-SNE visualization
ggplot(tsne_data, aes(x = V1, y = V2)) +
  geom_point() +
  ggtitle("t-SNE Visualization of Pancreatic Cancer Cells") +
  theme_minimal()

