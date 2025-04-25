# Build metaboset from internal standard-normalized data
modes_is_norm <- construct_metabosets(
  exprs = is_norm,
  pheno_data = imputed@phenoData@data,
  feature_data = imputed@featureData@data,
  group_col = "Group"
)

# Select specific mode (e.g., GC-MS method)
mode_is_norm <- modes_is_norm$Rtx5_EI

# "SummarizedExperiment" package installation
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

# Convert feature height table to SummarizedExperiment class
data <- SummarizedExperiment(assays = exprs(mode_is_norm),
                             colData = mode_is_norm@phenoData@data)


# Package for generalized logarithmic transform
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("pmp")
library(pmp)
# Generalised logarithmic transform
glog_exprs <- glog_transformation(df = data@assays@data@listData[[1]],
                                  classes = data$Group,
                                  qc_label = "QC")

# Adding glog transformation to notame MetaboSet
glog_set <- noflag
exprs(glog_set) <- glog_exprs

# Boxplot
bp <- plot_sample_boxplots(glog_set,
                                 order_by = "QC",
                                 fill_by = "QC")
# PCA
pca <- plot_pca(glog_set,
                      center = TRUE,
                      shape = "Group",
                      color = "Group",)
# Plot
pca + bp


# Perform tests for homoscedasticity
homoscedasticity_test_is <- perform_homoscedasticity_tests(
  glog_set,
  formula_char = "Feature ~ Group",
  all_features = TRUE
)

# Summary of results
total <- nrow(homoscedasticity_test_is)
passed <- sum(homoscedasticity_test_is$Fligner_P_FDR > 0.05)
failed <- sum(homoscedasticity_test_is$Fligner_P_FDR <= 0.05)

# Print results
cat("Features that PASSED homoscedasticity:", passed, "\n")
cat("Features that FAILED homoscedasticity:", failed, "\n")
cat("Total features evaluated:", total, "\n")


# Extract expression matrix from the metaboset object
exprs_is_norm_trans <- exprs(glog_set)

# Apply Shapiro-Wilk test to each feature (row)
shapiro_results_is_norm_trans <- apply(exprs_is_norm_trans, 1, function(x) shapiro.test(x)$p.value)

# Create a data frame with raw and adjusted p-values
normality_is <- data.frame(
  Feature_ID = rownames(exprs_is_norm_trans),
  Shapiro_p = shapiro_results_is_norm_trans,
  Shapiro_FDR = p.adjust(shapiro_results_is_norm_trans, method = "fdr")
)

# Summary of normality test results
total <- nrow(normality_is)
passed <- sum(normality_is$Shapiro_FDR > 0.05)
failed <- sum(normality_is$Shapiro_FDR <= 0.05)

# Print results
cat("Features that PASSED normality:", passed, "\n")
cat("Features that FAILED normality:", failed, "\n")
cat("Total features evaluated:", total, "\n")



#######################################

# Build metaboset from mass-normalized data
modes_mass_norm <- construct_metabosets(
  exprs = mass_norm,
  pheno_data = imputed@phenoData@data,
  feature_data = imputed@featureData@data,
  group_col = "Group"
)

# Select the relevant mode (e.g., Rtx5_EI method)
mode_mass_norm <- modes_mass_norm$Rtx5_EI


# Convert feature height table to SummarizedExperiment class
data_mass <- SummarizedExperiment(assays = exprs(mode_mass_norm),
                             colData = mode_mass_norm@phenoData@data)


# Generalised logarithmic transform
glog_exprs_mass <- glog_transformation(df = data_mass@assays@data@listData[[1]],
                                  classes = data_mass$Group,
                                  qc_label = "QC")

# Adding glog transformation to notame MetaboSet
glog_set_mass <- mode_mass_norm
exprs(glog_set_mass) <- glog_exprs_mass

# Boxplot
bp <- plot_sample_boxplots(glog_set_mass,,
                           order_by = "QC",
                           fill_by = "QC")
# PCA
pca <- plot_pca(glog_set_mass,
                center = TRUE,
                shape = "Group",
                color = "Group",)
# Plot
pca + bp


# Perform tests for homoscedasticity
homoscedasticity_test_mass <- perform_homoscedasticity_tests(
  glog_set_mass,
  formula_char = "Feature ~ Group",
  all_features = TRUE
)

# Summary of results
total <- nrow(homoscedasticity_test_mass)
passed <- sum(homoscedasticity_test_mass$Fligner_P_FDR > 0.05)
failed <- sum(homoscedasticity_test_mass$Fligner_P_FDR <= 0.05)

# Print results
cat("Features that PASSED homoscedasticity:", passed, "\n")
cat("Features that FAILED homoscedasticity:", failed, "\n")
cat("Total features evaluated:", total, "\n")


# Extract expression matrix from the metaboset object
exprs_mass_norm_trans <- exprs(glog_set_mass)

# Apply Shapiro-Wilk test to each feature (row)
shapiro_results_mass_norm_trans <- apply(exprs_mass_norm_trans, 1, function(x) shapiro.test(x)$p.value)

# Create a data frame with raw and adjusted p-values
normality_mass <- data.frame(
  Feature_ID = rownames(exprs_mass_norm_trans),
  Shapiro_p = shapiro_results_mass_norm_trans,
  Shapiro_FDR = p.adjust(shapiro_results_mass_norm_trans, method = "fdr")
)

# Summary of normality test results
total <- nrow(normality_is)
passed <- sum(normality_is$Shapiro_FDR > 0.05)
failed <- sum(normality_is$Shapiro_FDR <= 0.05)

# Print results
cat("Features that PASSED normality:", passed, "\n")
cat("Features that FAILED normality:", failed, "\n")
cat("Total features evaluated:", total, "\n")

##########################################33

# Build metaboset using IS + Mass normalized data
modes_is_mass_norm <- construct_metabosets(
  exprs = is_mass_norm,
  pheno_data = imputed@phenoData@data,
  feature_data = imputed@featureData@data,
  group_col = "Group"
)

# Select the analytical mode (e.g., Rtx5_EI)
mode_is_mass_norm <- modes_is_mass_norm$Rtx5_EI

# Convert feature height table to SummarizedExperiment class
data_is_mass <- SummarizedExperiment(assays = exprs(mode_is_mass_norm),
                                  colData = mode_is_mass_norm@phenoData@data)


# Generalised logarithmic transform
glog_exprs_is_mass <- glog_transformation(df = data_is_mass@assays@data@listData[[1]],
                                       classes = data_is_mass$Group,
                                       qc_label = "QC")

# Adding glog transformation to notame MetaboSet
glog_set_is_mass <- mode_mass_norm
exprs(glog_set_is_mass) <- glog_exprs_is_mass

# Boxplot
bp <- plot_sample_boxplots(glog_set_is_mass,,
                           order_by = "QC",
                           fill_by = "QC")
# PCA
pca <- plot_pca(glog_set_is_mass,
                center = TRUE,
                shape = "Group",
                color = "Group",)
# Plot
pca + bp


# Perform Levene's test for homoscedasticity
homoscedasticity_test_is_mass_t <- perform_homoscedasticity_tests(
  glog_set_is_mass,
  formula_char = "Feature ~ Group",
  all_features = TRUE
)

# Summarize results
total <- nrow(homoscedasticity_test_is_mass_t)
passed <- sum(homoscedasticity_test_is_mass_t$Fligner_P_FDR > 0.05)
failed <- sum(homoscedasticity_test_is_mass_t$Fligner_P_FDR <= 0.05)

# Print output
cat("Features that PASSED homoscedasticity:", passed, "\n")
cat("Features that FAILED homoscedasticity:", failed, "\n")
cat("Total features evaluated:", total, "\n")


# Extract expression matrix from the IS + Mass normalized metaboset
exprs_is_mass_norm <- exprs(glog_set_is_mass)

# Apply Shapiro-Wilk test to each feature (row-wise)
shapiro_results_is_mass_norm <- apply(exprs_is_mass_norm, 1, function(x) shapiro.test(x)$p.value)

# Compile results into a data frame with adjusted p-values (FDR)
normality_is_mass_norm <- data.frame(
  Feature_ID = rownames(exprs_is_mass_norm),
  Shapiro_p = shapiro_results_is_mass_norm,
  Shapiro_FDR = p.adjust(shapiro_results_is_mass_norm, method = "fdr")
)

# Summarize results
total <- nrow(normality_is_mass_norm)
passed <- sum(normality_is_mass_norm$Shapiro_FDR > 0.05)
failed <- sum(normality_is_mass_norm$Shapiro_FDR <= 0.05)

# Display summary
cat("Features that PASSED normality:", passed, "\n")
cat("Features that FAILED normality:", failed, "\n")
cat("Total features evaluated:", total, "\n")
