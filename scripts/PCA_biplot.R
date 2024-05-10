# Install and load the required packages
install.packages("ggplot2")
library(ggplot2)

# Read the count matrix from output.txt
count_matrix <- read.table("output.txt", sep = ",", header = TRUE, row.names = 1)

# Remove genes with zero counts across all samples
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

# Perform PCA on the count matrix without scaling
pca_result <- prcomp(t(count_matrix), scale. = FALSE)

# Extract the principal components (PCs)
pcs <- data.frame(pca_result$x)

# Add a group column to the data frame based on the sample names
pcs$group <- ifelse(grepl("CTL", rownames(pcs)), "CTL", "KO")

# Create a data frame for the legend
legend_df <- data.frame(
  group = c("CTL", "KO"),
  color = c("blue", "red")
)

# Create the PCA plot using ggplot2
pca_plot <- ggplot(pcs, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red")) +
  labs(title = "PCA Biplot",
       x = paste0("PC1: ", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "% variance"),
       y = paste0("PC2: ", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "% variance"),
       color = "Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

# Display the plot
print(pca_plot)
