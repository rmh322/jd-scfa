# INSTALLING AND LOADING PACKAGES 
packages <- c("vegan", "ggplot2", "dplyr", "colorRamps", "pheatmap", "openxlsx", "readxl", "tidyr",
              "writexl", "tibble", "ggrepel")
install.packages(packages[!packages %in% installed.packages()[,"Package"]])
lapply(packages, library, character.only = TRUE)

# SET WORKING DIRECTORY 
# starting with previously generated excel sheet
setwd("~/Documents/MSF Data/Jess MSc Project/Github Upload")
data_clean <- read_excel("jd-scfa_raw data.xlsx", sheet = 2)

# PREPARE DATA FOR PCA + BIPLOT
pca_data <- data_clean %>% 
  column_to_rownames("Protein") %>%
  t() %>% 
  as.data.frame()
# ADD GROUP INFORMATION
pca_groups <- data.frame(
  Sample = rownames(pca_data),
  Group = case_when(
    grepl("^A", rownames(pca_data)) ~ "LFD Sedentary",
    grepl("^B", rownames(pca_data)) ~ "LFD Exercise",
    grepl("^C", rownames(pca_data)) ~ "HFD Sedentary",
    grepl("^D", rownames(pca_data)) ~ "HFD Exercise"
  )
)
# MAKE SURE THAT GROUP NAMES ARE CORRECT 
print(unique(pca_groups$Group))

# PERFORM PCA
pca_result <- prcomp(pca_data, scale. = TRUE)
# PREPARE DATA FOR PLOTTING
pca_plot_data <- as.data.frame(pca_result$x) %>%
  cbind(pca_groups)
# CALCULATE VARIANCE EXPLAINED
var_explained <- summary(pca_result)$importance[2,] * 100
# CREATE PCA
proteomic_pca <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 6) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.1, size = 1) +
  scale_color_manual(values = c("LFD Sedentary" = "#FF7C7C", 
                                "LFD Exercise" = "#47C278", 
                                "HFD Sedentary" = "#ff8aff", 
                                "HFD Exercise" = "#8788ff")) +
  scale_fill_manual(values = c("LFD Sedentary" = "#FF7C7C", 
                               "LFD Exercise" = "#47C278", 
                               "HFD Sedentary" = "#ff8aff", 
                               "HFD Exercise" = "#8788ff")) +
  coord_fixed(ratio = 1.25) +
  theme_bw() +
  theme(aspect.ratio = 1.25,
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1.5),
        text = element_text(family = "Arial", size = 13),
        plot.title = element_text(family = "Arial", size = 13, face = "bold"),
        axis.title = element_text(family = "Arial", size = 30),
        axis.text = element_text(family = "Arial", size = 30, colour = "black")) +
  labs(title = "PCA of Total Proteome",
       x = paste0("PC1\n (", round(var_explained[1], 2), "% variance)"),
       y = paste0("PC2\n (", round(var_explained[2], 2), "% variance)"))
# DISPLAY THE PLOT AND SAVE 
print(proteomic_pca)
ggsave("Jess_Total Proteome PCA Plot.png", plot = proteomic_pca, width = 10, height = 8, units = "in", dpi = 300)

# RE-CREATE PCA TO INCLUDE SAMPLE IDs
# get unique group names
group_levels <- unique(pca_plot_data$Group)
# create a named vector for colours
group_colors <- setNames(c("#ff224b", "#74a8e3", "#ff8c00", "#228B22"), group_levels)
# CREATE PCA WITH SAMPLE IDs
pca_with_ids <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 0) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.1) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  coord_fixed(ratio = 1.25) +
  theme_bw() +
  theme(aspect.ratio = 1.5,
        panel.grid = element_blank(),
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(family = "Arial", size = 12, face = "bold"),
        axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 10)) +
  labs(title = "PCA of Proteome + IDs",
       x = paste0("PC1\n (", round(var_explained[1], 2), "% variance)"),
       y = paste0("PC2\n (", round(var_explained[2], 2), "% variance)"))
# DISPLAY THE PLOT AND SAVE
print(pca_with_ids)
ggsave("Jess_PCA Plot+ID.png", plot = pca_with_ids, width = 10, height = 8, units = "in", dpi = 300)


# GENERATE A BIPLOT 
# Remove the Protein column 
clean_prot <- data_clean
if("Protein" %in% colnames(clean_prot)) {
  clean_prot <- clean_prot[, -which(names(clean_prot) == "Protein")]
}
# Convert to numeric if not already
clean_prot <- as.data.frame(lapply(clean_prot, as.numeric), stringsAsFactors = FALSE)
# Perform PCA using vegan's rda function (can't figure out how to get around this)
pca <- rda(decostand(clean_prot, method = "hellinger"), scale = TRUE)
# Extract scores and loadings
scores_sites <- scores(pca, display = "sites", choices = c(1,2))
scores_species <- scores(pca, display = "species", choices = c(1,2))
# Create a data frame for sample scores
scores_sites_df <- as.data.frame(scores_sites)
scores_sites_df$Sample <- rownames(scores_sites_df)
# Create a data frame for variable loadings
scores_species_df <- as.data.frame(scores_species)
scores_species_df$Variable <- rownames(scores_species_df)

# CREATE THE BIPLOT
biplot <- ggplot() +
  geom_point(data = scores_sites_df, aes(x = PC1, y = PC2), size = 3) +
  geom_text_repel(data = scores_sites_df, aes(x = PC1, y = PC2, label = Sample), 
                  size = 3, max.overlaps = 10) +
  geom_segment(data = scores_species_df, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "#ff224b") +
  geom_text_repel(data = scores_species_df, aes(x = PC1, y = PC2, label = Variable), 
                  color = "#ff224b", size = 3, max.overlaps = 20) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(family = "Arial", size = 12, face = "bold"),
        axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 10)) +
  labs(title = "PCA Biplot",
       x = "PC1",
       y = "PC2")
# DISPLAY THE PLOT AND SAVE
print(biplot)
ggsave("Jess_PCA Biplot.png", plot = biplot, width = 10, height = 8, units = "in", dpi = 300)























