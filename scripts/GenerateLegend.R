library(ComplexHeatmap)
library(grid)

# Define colors for groups
group_colors <- c("BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
                  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
                  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1],
                  "Sig17-" = wes_palette("AsteroidCity3")[2])

# Create a standalone legend
phenotype_legend <- Legend(
  title = "Phenotype",
  at = names(group_colors),
  labels = names(group_colors),
  legend_gp = gpar(fill = group_colors),
  direction = "horizontal",  # Ensure left-to-right layout
  nrow = 1,  # Force single-row layout
  labels_gp = gpar(fontsize = 10),  # Adjust font size
  legend_width = unit(10, "cm")  # Increase width for spacing
)

# Open a new PDF or plot window
pdf("/Users/jao/Desktop/Phenotype_Legend.pdf", width = 8, height = 2)
grid.draw(phenotype_legend)  # Draw only the legend
dev.off()



# Define colors for groups
group_colors <- c("BarrettsSig17+" = wes_palette("IsleofDogs1")[1],
                  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
                  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
                  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1],
                  "Sig17-" = wes_palette("AsteroidCity3")[2])

# Create a standalone legend
phenotype_legend <- Legend(
  title = "Phenotype",
  at = names(group_colors),
  labels = names(group_colors),
  legend_gp = gpar(fill = group_colors),
  direction = "horizontal",  # Ensure left-to-right layout
  nrow = 1,  # Force single-row layout
  labels_gp = gpar(fontsize = 10),  # Adjust font size
  legend_width = unit(10, "cm")  # Increase width for spacing
)

# Open a new PDF or plot window
pdf("/Users/jao/Desktop/Phenotype_LegendwithBE.pdf", width = 8, height = 2)
grid.draw(phenotype_legend)  # Draw only the legend
dev.off()
