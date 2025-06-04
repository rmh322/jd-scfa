# INSTALLING AND LOADING PACKAGES 
packages <- c("vegan", "ggplot2", "dplyr", "colorRamps", "pheatmap", "openxlsx", "readxl", "tidyr",
              "writexl", "tibble", "ggrepel")
install.packages(packages[!packages %in% installed.packages()[,"Package"]])
lapply(packages, library, character.only = TRUE)

# SET WORKING DIRECTORY 
# starting with previously generated excel sheet
setwd("~/Documents/MSF Data/Jess MSc Project/Github Upload")
data_clean <- read_excel("jd-scfa_raw data.xlsx", sheet = 2)

# PREPARE DATA FOR VOLC PLOT AND ADD GROUP INFORMATION
long_data <- data_clean %>%
  pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Expression") %>%
  mutate(
    Group = case_when(
      grepl("^A", Sample) ~ "LFD_Sedentary",
      grepl("^B", Sample) ~ "LFD_Exercise",
      grepl("^C", Sample) ~ "HFD_Sedentary",
      grepl("^D", Sample) ~ "HFD_Exercise"
    )
  )
# Print the first few rows of the reshaped data
print(head(long_data))

# CALCULATE MEAN EXPRESSION AND FOLD CHANGES
mean_expression <- long_data %>%
  group_by(Protein, Group) %>%
  summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
fold_changes <- mean_expression %>%
  pivot_wider(names_from = Group, values_from = MeanExpression) %>%
  mutate(
    FC_LFD_Exercise = LFD_Exercise / LFD_Sedentary,
    FC_HFD_Sedentary = HFD_Sedentary / LFD_Sedentary,
    FC_HFD_Exercise = HFD_Exercise / LFD_Sedentary
  )
# Print the first few rows of fold changes
print(head(fold_changes))

# PERFORM STATISTICAL ANALYSIS
perform_t_test <- function(group_data, control_data) {
  tryCatch({
    t.test(group_data, control_data)$p.value
  }, error = function(e) NA)
}

p_values <- long_data %>%
  group_by(Protein) %>%
  summarise(
    p_LFD_Exercise = perform_t_test(Expression[Group == "LFD_Exercise"], Expression[Group == "LFD_Sedentary"]),
    p_HFD_Sedentary = perform_t_test(Expression[Group == "HFD_Sedentary"], Expression[Group == "LFD_Sedentary"]),
    p_HFD_Exercise = perform_t_test(Expression[Group == "HFD_Exercise"], Expression[Group == "LFD_Sedentary"])
  )

# COMBINE FOLD-CHANGSE + P-VALUES
results <- fold_changes %>%
  left_join(p_values, by = "Protein") %>%
  pivot_longer(
    cols = c(starts_with("FC_"), starts_with("p_")),
    names_to = c(".value", "Group"),
    names_pattern = "(FC|p)_(.*)"
  ) %>%
  mutate(
    log2FoldChange = log2(FC),
    neg_log10_pvalue = -log10(p),
    significant = p < 0.05 & !is.na(p)
  )
# Print the first few rows of results
print(head(results))

# MODIFY RESULTS DATA FRAME TO INCLUDE COLOURS
results <- results %>%
  mutate(
    color_category = case_when(
      abs(log2FoldChange) <= 1 ~ "Not Significant",
      significant ~ Group,
      TRUE ~ "Not Significant"
    )
  )

# CREATE VOLCANO PLOT (NO LABELS)
volcplot <- ggplot(results, aes(x = log2FoldChange, y = neg_log10_pvalue, color = Group)) +
  geom_point(size = 6, aes(alpha = significant)) +
  scale_color_manual(values = c("LFD_Exercise" = "#47C278", 
                                "HFD_Sedentary" = "#ff8aff", 
                                "HFD_Exercise" = "#8788ff")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 1, color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 1, color = "grey60") +
  coord_fixed(ratio = 1.25) +
  labs(
    title = "Volcano Plot v. LFD Sedentary",
    x = "log2\n (Fold Change)",
    y = "-log10\n (p-value)"
  ) +
  theme_bw() +
  theme(aspect.ratio = 1.25,
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1.5),
        text = element_text(family = "Arial", size = 13),
        plot.title = element_text(family = "Arial", size = 13, face = "bold"),
        axis.title = element_text(family = "Arial", size = 30),
        axis.text = element_text(family = "Arial", size = 30, colour = "black"),
        legend.position = "none"
  ) +
  geom_text_repel(
    data = subset(results, significant & abs(log2FoldChange) > 1),
    aes(label = Protein),
    size = 10,
    box.padding = 0.5,
    max.overlaps = 0
  ) +
  theme(legend.position = "none")
  
# DISPLAY THE PLOT AND SAVE
print(volcplot)
ggsave("Jess_All Proteome Volc Plot.png", plot = volcplot, width = 10, height = 8, units = "in", dpi = 300)

# CREATE A VECTOR OF PROTEINS YOU WANT TO HIGHLIGHT
proteins_of_interest <- c("THIKA", "THIM", "E9PSQ0", "ACD11", "ACADL", "ACADM", "ACADS", "ACADV", "THIL",
                         "ACSL1", "ACSL6", "CPT1B", "CPT2", "CACP", "ESHM", "ECI1", "ECI2", "ETFA", "ETFB", 
                         "ETFD", "HCDH", "ECHA", "ECHB", "HCD2", "PCCA", "MCAT", "SCOT1")
# Define color palette
group_colors <- c("LFD_Exercise" = "#74a8e3", "HFD_Sedentary" = "#228B22", "HFD_Exercise" = "#ff224b")
# Modify the volcano plot code
volcplot <- ggplot(results, aes(x = log2FoldChange, y = neg_log10_pvalue, color = Group)) +
  geom_point(aes(alpha = significant)) +
  scale_color_manual(values = c("LFD_Exercise" = "#74a8e3", "HFD_Sedentary" = "#228B22", "HFD_Exercise" = "#ff224b")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  labs(
    title = "Volcano Plot v. LFD Sedentary",
    x = "log2 (Fold Change)",
    y = "-log10 (P-value)"
  ) +
  theme_minimal() +
  theme(
    aspect.ratio = 1.5,
    panel.grid = element_blank(),
    text = element_text(family = "Arial", size = 25),
    plot.title = element_text(family = "Arial", size = 25, face = "bold"),
    axis.title = element_text(family = "Arial", size = 25),
    axis.text = element_text(family = "Arial", size = 20),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  # Add labels for proteins of interest
  geom_text_repel(
    data = subset(results, Protein %in% proteins_of_interest),
    aes(label = Protein, color = "black"),
    size = 4,
    box.padding = 0.7,
    point.padding = 0.5,
    force = 10,
    max.overlaps = 20  # Increased to allow more labels
  )

# DISPLAY THE PLOT AND SAVE
print(volcplot)
ggsave("Jess_All Pr Volc Plot_Highlighted.png", plot = volcplot, width = 10, height = 8, units = "in", dpi = 300)






# IDENTIFY TOP DIFFERENTIALLY EXPRESSED PROTEINS FOR EACH GROUP
top_proteins <- results %>%
  filter(significant) %>%
  group_by(Group) %>%
  top_n(3, abs(log2FoldChange)) %>%
  ungroup()

# CREATE VOLCANO PLOT (WITH LABELS)
volcplot_label <- ggplot(results, aes(x = log2FoldChange, y = neg_log10_pvalue, color = Group)) +
  geom_point(aes(alpha = significant)) +
  scale_color_manual(values = c("LFD_Exercise" = "#74a8e3", "HFD_Sedentary" = "#ff8c00", "HFD_Exercise" = "#228B22")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot v. LFD Sedentary",
    x = "log2 (Fold Change)",
    y = "-log10 (P-value)"
  ) +
  theme_minimal() +
  theme(aspect.ratio = 1.5,
        panel.grid = element_blank(),
        text = element_text(family = "Arial", size = 12),
        plot.title = element_text(family = "Arial", size = 12, face = "bold"),
        axis.title = element_text(family = "Arial", size = 12),
        axis.text = element_text(family = "Arial", size = 10)) +
  theme(legend.position = "none") +
  geom_text_repel(
    data = top_proteins,
    aes(label = Protein),
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf,
    force = 10
  )
# DISPLAY THE PLOT AND SAVE
print(volcplot_label)
ggsave("Jess_All Pr Volc Plot + Label.png", plot = volcplot, width = 10, height = 8, units = "in", dpi = 300)

# Print and save summary of top differentially expressed proteins
print(top_proteins)
write_xlsx(top_proteins, "Jess_All Pr (Top Differentially Expressed).xlsx")

# Print summary of significant proteins
significant_proteins <- results %>%
  filter(significant & abs(log2FoldChange) > 1) %>%
  arrange(Group, desc(abs(log2FoldChange)))

print(significant_proteins)
write_xlsx(significant_proteins, "Jess_All Pr (Significant Proteins).xlsx")