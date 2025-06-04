# INSTALLING AND LOADING PACKAGES 
packages <- c("vegan", "ggplot2", "dplyr", "colorRamps", "pheatmap", "openxlsx", "readxl", "tidyr",
              "writexl", "tibble", "ggrepel")
install.packages(packages[!packages %in% installed.packages()[,"Package"]])
lapply(packages, library, character.only = TRUE)

# SET WORKING DIRECTORY 
setwd("~/Documents/MSF Data/Jess MSc Project/Github Upload")
jdata <- read_excel("jd-scfa_raw data.xlsx", sheet = 2)

# CLEAN DATA SET TO INCLUDE ONLY PROTEINS EXPRESSED IN ALL TISSUES
data_clean <- na.omit(jdata) # clean to 454
str(data_clean)
print(paste("# of Pr after cleaning:", nrow(data_clean)))
# save "data_clean" as new excel sheet to use for next code
write_xlsx(data_clean, "Jess_Clean Protein List.xlsx") # copy this list into original excel sheet & delete ... 

# RESHAPE DATA FOR ANALYSIS
long_data <- data_clean %>%
  pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Expression") %>%
  mutate(
    Diet = ifelse(grepl("^[AB]", Sample), "LFD", "HFD"),
    Exercise = ifelse(grepl("^[AC]", Sample), "Sedentary", "Exercise"),
    Group = case_when(
      grepl("^A", Sample) ~ "LFD_Sedentary",
      grepl("^B", Sample) ~ "LFD_Exercise",
      grepl("^C", Sample) ~ "HFD_Sedentary",
      grepl("^D", Sample) ~ "HFD_Exercise"
    )
  )
table(long_data$Group, long_data$Sample)

# PERFORM APPROPRIATE STATISTICAL ANALYSIS (two-way ANOVA here)
analyze_protein <- function(df) {
  model <- aov(Expression ~ Diet * Exercise, data = df)
  summary_stats <- summary(model)[[1]]
  
  # CALCULATE FOLD CHANGES
  fc_diet <- mean(df$Expression[df$Diet == "HFD"]) / mean(df$Expression[df$Diet == "LFD"])
  fc_exercise <- mean(df$Expression[df$Exercise == "Exercise"]) / mean(df$Expression[df$Exercise == "Sedentary"])
  
  data.frame(
    Protein = df$Protein[1],
    p_Diet = summary_stats["Diet", "Pr(>F)"],
    p_Exercise = summary_stats["Exercise", "Pr(>F)"],
    p_Interaction = summary_stats["Diet:Exercise", "Pr(>F)"],
    fc_Diet = fc_diet,
    fc_Exercise = fc_exercise
  )
}
# APPLY THE ANALYSIS FOR EACH PROTEIN 
results <- long_data %>%
  group_by(Protein) %>%
  do(analyze_protein(.))
# ADJUST p-values FOR MULTIPLE COMPARISONS
results <- results %>%
  mutate(
    across(starts_with("p_"), ~p.adjust(., method = "fdr"), .names = "{.col}_adj")
  )

# CREATE A DATA FRAME WITH SIGNIFICANT PROTEINS ONLY
significant_proteins <- results %>%
  filter(p_Diet_adj < 0.05 | p_Exercise_adj < 0.05 | p_Interaction_adj < 0.05)


# GENERATE DATA FOR THE VOLCANO PLOT
volcano_data <- results %>%
  mutate(
    fcvolc_Diet = log2(fc_Diet),
    fcvolc_Exercise = log2(fc_Exercise),
    pvolc_Diet = -log10(p_Diet_adj),
    pvolc_Exercise = -log10(p_Exercise_adj),
    pvolc_Interaction = -log10(p_Interaction_adj)
  )
write_xlsx(volcano_data, "Jess_Volcano Plot Data.xlsx")

# GENERATE AND SAVE NEW EXCEL FILES WITH INFORMATION GENERATED
full_list <- left_join(data_clean, results, by = "Protein") %>%
  mutate(
    fcvolc_Diet = log2(fc_Diet),
    fcvolc_Exercise = log2(fc_Exercise),
    pvolc_Diet = -log10(p_Diet_adj),
    pvolc_Exercise = -log10(p_Exercise_adj),
    pvolc_Interaction = -log10(p_Interaction_adj)
  )
write_xlsx(full_list, "Jess_Clean Protein List + Stats.xlsx")
write_xlsx(significant_proteins, "Jess_Significant Pr.xlsx")

# DONE HERE - NOW NEED TO VISUALIZE DATA W DIFFERENT SCRIPT


