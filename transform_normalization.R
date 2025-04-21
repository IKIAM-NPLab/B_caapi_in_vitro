library(ggplot2)

# ------------------------------
# Log Transformation of Raw Data
# ------------------------------

# Extract the feature table from the imputed dataset
feat_table <- exprs(imputed)

# Apply log1p transformation to improve feature distribution
feat_table_log <- log1p(feat_table)

# ------------------------------
# Internal Standard Normalization
# ------------------------------

# Extract feature data from clustered object (prior to clustering compression)
feat_table_is <- exprs(clustered)

# Retrieve internal standard (e.g., caffeine) as a reference vector
caffeine_vector <- feat_table_is["Rtx5_EI_194_052603269692a25_099312", ]

# Normalize log-transformed data by internal standard
is_norm_log <- sweep(feat_table_log, 2, caffeine_vector, "/")

# ------------------------------
# Homoscedasticity & Normality Test (IS normalization)
# ------------------------------

# Create a metaboset object
modes_is_norm_log <- construct_metabosets(
  exprs = is_norm_log,
  pheno_data = imputed@phenoData@data,
  feature_data = imputed@featureData@data,
  group_col = "Group"
)

# Extract specific mode
mode_is_norm_log <- modes_is_norm_log$Rtx5_EI

# Perform tests for homoscedasticity
homoscedasticity_test_is_log <- perform_homoscedasticity_tests(
  mode_is_norm_log, formula_char = "Feature ~ Group", all_features = TRUE
)

# Summary results
total <- nrow(homoscedasticity_test_is_log)
passed <- sum(homoscedasticity_test_is_log$Bartlett_P_FDR > 0.05)
failed <- sum(homoscedasticity_test_is_log$Bartlett_P_FDR <= 0.05)

cat("Features that PASSED homoscedasticity:", passed, "\n")
cat("Features that FAILED homoscedasticity:", failed, "\n")
cat("Total features evaluated:", total, "\n")

# Shapiro-Wilk test for normality
exprs_is_norm_log <- exprs(mode_is_norm_log)
shapiro_pvalues <- apply(exprs_is_norm_log, 1, function(x) shapiro.test(x)$p.value)

normality_is_log <- data.frame(
  Feature_ID = rownames(exprs_is_norm_log),
  Shapiro_p = shapiro_pvalues,
  Shapiro_FDR = p.adjust(shapiro_pvalues, method = "fdr")
)

# Summary of normality test
total <- nrow(normality_is_log)
passed <- sum(normality_is_log$Shapiro_FDR > 0.05)
failed <- sum(normality_is_log$Shapiro_FDR <= 0.05)

cat("Features that PASSED normality:", passed, "\n")
cat("Features that FAILED normality:", failed, "\n")
cat("Total features evaluated:", total, "\n")


#plotting the results

# Create summary data for both tests
comparison_df <- data.frame(
  Test = rep(c("Homoscedasticity", "Normality"), each = 2),
  Result = rep(c("Passed", "Failed"), 2),
  Count = c(
    sum(homoscedasticity_test_is_log$Bartlett_P_FDR > 0.05),
    sum(homoscedasticity_test_is_log$Bartlett_P_FDR <= 0.05),
    sum(normality_is_log$Shapiro_FDR > 0.05),
    sum(normality_is_log$Shapiro_FDR <= 0.05)
  )
)

# Plot
ggplot(comparison_df, aes(x = Test, y = Count, fill = Result)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Comparison of Homoscedasticity and Normality Tests",
    x = "Statistical Test",
    y = "Number of Features"
  ) +
  scale_fill_manual(values = c("Passed" = "#4CAF50", "Failed" = "#F44336")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )


# ------------------------------
# Mass Normalization
# ------------------------------

# Use previously transformed log1p feature table
# Extract sample mass values
masses <- imputed@phenoData@data$Mass

# Normalize by sample mass
mass_norm_log <- sweep(feat_table_log, 2, masses, "/")

# Create metaboset for mass-normalized data
modes_mass_norm_log <- construct_metabosets(
  exprs = mass_norm_log,
  pheno_data = imputed@phenoData@data,
  feature_data = imputed@featureData@data,
  group_col = "Group"
)

mode_mass_norm_log <- modes_mass_norm_log$Rtx5_EI

# Homoscedasticity test
homoscedasticity_test_mass_norm_log <- perform_homoscedasticity_tests(
  mode_mass_norm_log, formula_char = "Feature ~ Group", all_features = TRUE
)

# Summary
total <- nrow(homoscedasticity_test_mass_norm_log)
passed <- sum(homoscedasticity_test_mass_norm_log$Bartlett_P_FDR > 0.05)
failed <- sum(homoscedasticity_test_mass_norm_log$Bartlett_P_FDR <= 0.05)

cat("Features that PASSED homoscedasticity:", passed, "\n")
cat("Features that FAILED homoscedasticity:", failed, "\n")
cat("Total features evaluated:", total, "\n")

# Normality test
exprs_mass_norm_log <- exprs(mode_mass_norm_log)
shapiro_pvalues <- apply(exprs_mass_norm_log, 1, function(x) shapiro.test(x)$p.value)

normality_mass_norm_log <- data.frame(
  Feature_ID = rownames(exprs_mass_norm_log),
  Shapiro_p = shapiro_pvalues,
  Shapiro_FDR = p.adjust(shapiro_pvalues, method = "fdr")
)

# Summary of normality test
total <- nrow(normality_mass_norm_log)
passed <- sum(normality_mass_norm_log$Shapiro_FDR > 0.05)
failed <- sum(normality_mass_norm_log$Shapiro_FDR <= 0.05)

cat("Features that PASSED normality:", passed, "\n")
cat("Features that FAILED normality:", failed, "\n")
cat("Total features evaluated:", total, "\n")

#plotting the results
# Create summary data for both tests
comparison_df <- data.frame(
  Test = rep(c("Homoscedasticity", "Normality"), each = 2),
  Result = rep(c("Passed", "Failed"), 2),
  Count = c(
    sum(homoscedasticity_test_mass_norm_log$Bartlett_P_FDR > 0.05),
    sum(homoscedasticity_test_mass_norm_log$Bartlett_P_FDR <= 0.05),
    sum(normality_mass_norm_log$Shapiro_FDR > 0.05),
    sum(normality_mass_norm_log$Shapiro_FDR <= 0.05)
  )
)

# Plot
ggplot(comparison_df, aes(x = Test, y = Count, fill = Result)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Comparison of Homoscedasticity and Normality Tests",
    x = "Statistical Test",
    y = "Number of Features"
  ) +
  scale_fill_manual(values = c("Passed" = "#4CAF50", "Failed" = "#F44336")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

# ------------------------------
# Combined Normalization (IS + Mass)
# ------------------------------

# Normalize internal standard-normalized data by mass
is_mass_norm_log <- sweep(is_norm_log, 2, masses, "/")

# Construct metaboset
modes_is_mass_norm_log <- construct_metabosets(
  exprs = is_mass_norm_log,
  pheno_data = imputed@phenoData@data,
  feature_data = imputed@featureData@data,
  group_col = "Group"
)

mode_is_mass_norm_log <- modes_is_mass_norm_log$Rtx5_EI

# Homoscedasticity test
homoscedasticity_test_is_mass_log <- perform_homoscedasticity_tests(
  mode_is_mass_norm_log, formula_char = "Feature ~ Group", all_features = TRUE
)

# Summary
total <- nrow(homoscedasticity_test_is_mass_log)
passed <- sum(homoscedasticity_test_is_mass_log$Bartlett_P_FDR > 0.05)
failed <- sum(homoscedasticity_test_is_mass_log$Bartlett_P_FDR <= 0.05)

cat("Features that PASSED homoscedasticity:", passed, "\n")
cat("Features that FAILED homoscedasticity:", failed, "\n")
cat("Total features evaluated:", total, "\n")

# Normality test
exprs_is_mass_norm_log <- exprs(mode_is_mass_norm_log)
shapiro_pvalues <- apply(exprs_is_mass_norm_log, 1, function(x) shapiro.test(x)$p.value)

normality_is_mass_norm <- data.frame(
  Feature_ID = rownames(exprs_is_mass_norm_log),
  Shapiro_p = shapiro_pvalues,
  Shapiro_FDR = p.adjust(shapiro_pvalues, method = "fdr")
)

# Summary of normality test
total <- nrow(normality_is_mass_norm)
passed <- sum(normality_is_mass_norm$Shapiro_FDR > 0.05)
failed <- sum(normality_is_mass_norm$Shapiro_FDR <= 0.05)

cat("Features that PASSED normality:", passed, "\n")
cat("Features that FAILED normality:", failed, "\n")
cat("Total features evaluated:", total, "\n")

#plotting the results
# Create summary data for both tests
comparison_df <- data.frame(
  Test = rep(c("Homoscedasticity", "Normality"), each = 2),
  Result = rep(c("Passed", "Failed"), 2),
  Count = c(
    sum(homoscedasticity_test_is_mass_log$Bartlett_P_FDR > 0.05),
    sum(homoscedasticity_test_is_mass_log$Bartlett_P_FDR <= 0.05),
    sum(normality_is_mass_norm$Shapiro_FDR > 0.05),
    sum(normality_is_mass_norm$Shapiro_FDR <= 0.05)
  )
)

# Plot
ggplot(comparison_df, aes(x = Test, y = Count, fill = Result)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Comparison of Homoscedasticity and Normality Tests",
    x = "Statistical Test",
    y = "Number of Features"
  ) +
  scale_fill_manual(values = c("Passed" = "#4CAF50", "Failed" = "#F44336")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )
