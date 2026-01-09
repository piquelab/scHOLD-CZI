##################################################### DEGs LPS #####################################
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(stringr)

rm(list=ls())

# Define the significance threshold and DEG count threshold
ALPHA <- 0.1
MIN_DEGS_THRESHOLD <- 50

# Set directories
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/LPS_enrich/DEGs/")
outDir_save <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/LPS_enrich/DEGs/contingency"

# Define variables and directions
ps_vars <- c("PSS_all_mean", "ISEL_Mean")
ps_vars_short <- c("PSS", "ISEL")
directions <- list(
    c("Up", "Up"),   # Var Up vs LPS Up
    c("Down", "Up"), # Var Down vs LPS Up
    c("Up", "Down"), # Var Up vs LPS Down
    c("Down", "Down")# Var Down vs LPS Down
)

# --- 0.1 LOAD AND PREPARE RNA-SEQ DATA (DEGS) ---

# Load Psychosocial DEGS (res)
fname <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/allres_allvars_noCombat_DESeq.txt"
res <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")
# Assuming 'gene' is the column containing gene identifiers in the file
res <- res %>% filter(!cluster=="C5")

# Load LPS DEGS (treats)

# load LPS degs
fname=paste0("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/treatment_noCOMBAT/deseqres/ALL.0.1.13.deseqres.RNA-LPS_vs_RNA-CTRL.treatment_noCOMBAT.txt")
treats <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")
# Assuming 'gene' is the column containing gene identifiers in the file
treats <- treats %>% filter(!cluster=="C5")

# --- APPLY FILTER C5 (R5 DC) as per previous script logic ---
res <- res %>% filter(!cluster=="C5")
treats <- treats %>% filter(!cluster=="C5")


# Define RNA-specific cluster and cell type map
rna_cell_type_map <- c(
    "C0" = "R0 T CD4+",
    "C1" = "R1 T CD8+",
    "C2" = "R2 NK",
    "C3" = "R3 Monocyte",
    "C4" = "R4 B" # C5 removed, so C4 is the last one
)

# Define factor levels and color palette for plotting (using a subset of colors)
rna_cell_type_levels <- unname(rna_cell_type_map)
# Removed one color as C5 is filtered out
rna_colors_array <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#AA4B56") 
rna_colors <- setNames(rna_colors_array, rna_cell_type_levels)

# Define the clusters to analyze (based on LPS data)
clusters_to_analyze <- unique(treats$cluster)
clusters_to_analyze <- sort(clusters_to_analyze)


# --- 1.5 Determine Variable/Cluster Pairs with > MIN_DEGS_THRESHOLD Significant DEGs (padj < ALPHA) ---

# Calculate the count of significant DEGs for each variable-cluster pair
deg_counts_filter <- res %>% 
    # Must use padj < ALPHA to define "significant" for the count threshold
    filter(padj < ALPHA, var %in% ps_vars) %>%
    group_by(var, cluster) %>%
    # CRITICAL: Ensure we only count unique gene identifiers within each group
    distinct(identifier, .keep_all = TRUE) %>% 
    summarise(Sig_DEG_Count = n(), .groups = 'drop') %>%
    filter(Sig_DEG_Count >= MIN_DEGS_THRESHOLD)

# --- 1. Core Fisher Test Function for RNA-Seq (DEGs) ---

#' Runs a two-sided Fisher's Exact Test for enrichment on cluster-specific RNA-Seq data.
run_fisher_test_dynamic_rna <- function(ps_data, ps_dir_sign, lps_dir_sign, test_name, treats_data, cluster_name, ps_variable) {

  # --- A. Create Binary Status for PS Variable ---
  ps_status <- ps_data %>%
    # Use 'identifier' (gene) as the unique key
    select(identifier, padj, logFC) %>%
    rename(padj_Var = padj, logFC_Var = logFC) %>%
    mutate(
      Var_Status = ifelse(
        # Check for significant and correct direction (this is the filter)
        padj_Var < ALPHA & (if (ps_dir_sign == "Up") logFC_Var > 0 else logFC_Var < 0),
        TRUE,
        FALSE
      )
    )

  # --- B. Create Binary Status for LPS Treatment ---
  lps_status <- treats_data %>%
    # Use 'identifier' (gene) as the unique key
    select(identifier, padj, logFC) %>%
    rename(padj_LPS = padj, logFC_LPS = logFC) %>%
    mutate(
      LPS_Status = ifelse(
        # Check for significant and correct direction (this is the filter)
        padj_LPS < ALPHA & (if (lps_dir_sign == "Up") logFC_LPS > 0 else logFC_LPS < 0),
        TRUE,
        FALSE
      )
    )

  # --- C. Merge Datasets (Common Universe of Genes within this Cluster) ---
  merged_status <- inner_join(
    ps_status,
    lps_status,
    by = "identifier" # Gene ID is the common key
  )

  # --- D. Calculate the Four Cell Counts (A, B, C, D) ---
  A <- merged_status %>% filter(Var_Status == TRUE & LPS_Status == TRUE) %>% nrow()
  B <- merged_status %>% filter(Var_Status == TRUE & LPS_Status == FALSE) %>% nrow()
  C <- merged_status %>% filter(Var_Status == FALSE & LPS_Status == TRUE) %>% nrow()
  D <- merged_status %>% filter(Var_Status == FALSE & LPS_Status == FALSE) %>% nrow()

  # --- E. Construct the Contingency Table (with descriptive labels) ---
  var_target_label = paste0(ps_variable, "_", ps_dir_sign, "_DEGs")
  var_not_target_label = paste0("NOT_", ps_variable, "_", ps_dir_sign, "_DEGs")
  lps_target_label = paste0("LPS_", lps_dir_sign, "_DEGs")
  lps_not_target_label = paste0("NOT_LPS_", lps_dir_sign, "_DEGs")

  contingency_matrix <- matrix(
    c(A, C, B, D), # Fill column-wise
    nrow = 2,
    dimnames = list(
      Var_Status = c(var_target_label, var_not_target_label),
      LPS_Status = c(lps_target_label, lps_not_target_label)
    )
  )

  # --- F. Run Two-Sided Fisher's Exact Test ---
  fisher_result <- tryCatch({
    fisher.test(contingency_matrix, alternative = "two.sided")
  }, error = function(e) {
    list(p.value = NA, estimate = NA, conf.int = c(NA, NA))
  })


  # --- G. Format Results ---
  results_df <- data.frame( 
    Test_Name = test_name,
    Cluster = cluster_name,
    Var_Variable = ps_variable,
    Var_Target_Dir = ps_dir_sign,
    LPS_Target_Dir = lps_dir_sign,
    N_Enrichment = A,
    N_Universe = A + B + C + D, # This is the large number (all genes tested)
    P_Value = fisher_result$p.value,
    Odds_Ratio = fisher_result$estimate,
    CI_Lower = fisher_result$conf.int[1],
    CI_Upper = fisher_result$conf.int[2],
    stringsAsFactors = FALSE
  )

  return(list(stats = results_df, matrix = contingency_matrix))
}


# Initialize lists to store results
all_stats_list <- list()
all_matrices_list <- list()
counter <- 1

# --- 2. Loop Through Clusters, Variables, and Directions (Filtered Tests) ---

# Loop over each cluster
for (current_cluster_name in clusters_to_analyze) { # RENAMED: Use current_cluster_name

    # 2.1 Filter LPS data for the current cluster once
    lps_subset_raw <- treats %>%
        filter(cluster == current_cluster_name) %>% # Use the new distinct variable
        select(identifier, padj, logFC) %>%
        distinct(identifier, .keep_all = TRUE) 

    # Loop over each psychosocial variable (PSS and ISEL)
    for (ps_var_idx in seq_along(ps_vars)) {
        ps_var <- ps_vars[ps_var_idx]
        ps_var_short <- ps_vars_short[ps_var_idx]

        # CRITICAL FIX: Check if this Variable-Cluster pair meets the MIN_DEGS_THRESHOLD
        current_count_row <- deg_counts_filter %>%
            filter(cluster == current_cluster_name, var == ps_var) # Use the new distinct variable
        
        if (nrow(current_count_row) == 0) {
            # This pair either has 0 significant DEGs or is below the threshold
            cat(paste("Skipping:", current_cluster_name, "|", ps_var_short, " - Below DEG threshold (<", MIN_DEGS_THRESHOLD, ").\n"))
            next
        }

        # 2.2 Filter Var data for the current variable AND cluster
        ps_subset_raw <- res %>%
            filter(var == ps_var, cluster == current_cluster_name) %>% # Use the new distinct variable
            select(identifier, padj, logFC) %>%
            distinct(identifier, .keep_all = TRUE) 

        # Loop over the 4 directional combinations
        for (dir_pair in directions) {
            ps_dir <- dir_pair[1]
            lps_dir <- dir_pair[2]

            # Test Name construction
            test_name <- paste0("DEG_", current_cluster_name, "_", ps_var_short, ps_dir, "_LPS", lps_dir)

            cat(paste("Running Test", counter, ":", test_name, "\n"))
            
            # DIAGNOSTIC: Check Raw Significant Count in the PS data before merging.
            # This value should equal A + B in the final contingency table.
            ps_raw_sig_count <- ps_subset_raw %>%
                filter(
                    padj < ALPHA, 
                    if (ps_dir == "Up") logFC > 0 else logFC < 0
                ) %>%
                nrow()
            cat(paste("    -> Raw PSS/ISEL Significant Genes (Expected A+B):", ps_raw_sig_count, "\n"))


            # Run the dynamic test function
            test_output <- run_fisher_test_dynamic_rna(
                ps_data = ps_subset_raw,
                ps_dir_sign = ps_dir,
                lps_dir_sign = lps_dir,
                test_name = test_name,
                treats_data = lps_subset_raw,
                cluster_name = current_cluster_name,
                ps_variable = ps_var_short
            )

            # Store the results
            all_stats_list[[counter]] <- test_output$stats

            # Store the matrix separately
            all_matrices_list[[test_name]] <- test_output$matrix
            counter <- counter + 1
        }
    }
}

# --- 3. Final Output and Saving ---

# Combine all statistics into a single dataframe
final_results_deg <- bind_rows(all_stats_list)

# 3.1 Save Combined DEG Results
print("--- SUMMARY OF ALL CELL-TYPE-SPECIFIC DEG FISHER'S EXACT TESTS (RNA-SEQ) ---")
print(final_results_deg)
fname_stats_deg <- paste0(outDir_save, "25.3.summary_DEG_fisher_tests_by_cluster.txt")
write.table(final_results_deg, file = fname_stats_deg, sep = "\t", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)
cat(paste0("\nSaved combined DEG summary statistics to: ", fname_stats_deg, "\n"))


# 3.2 Save all contingency matrices individually
for (test_name in names(all_matrices_list)) {
    mat <- all_matrices_list[[test_name]]
    # Output path confirmed to be one level up as requested previously: outDir_save / matrix_...
    fname_matrix <- paste0(outDir_save, "/", "matrix_", test_name, ".txt")

    df_to_save <- as.data.frame(mat)
    # The row labels (Variable_Status) are converted to a column here
    df_to_save <- cbind(Variable_Status = rownames(mat), df_to_save)
    # Ensure nested directory is created if necessary (though usually handled by R)
    dir.create(dirname(fname_matrix), showWarnings = FALSE, recursive = TRUE)
    write.table(df_to_save, file = fname_matrix, sep = "\t", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)
}
cat(paste0("Saved ", length(all_matrices_list), " individual contingency matrices to: ", outDir_save, "/matrix_DEG_*.txt\n"))


# --- 4. Generate Integrated Forest Plot (Log Odds Ratio) ---

# Function to generate and save the integrated forest plot for DEGs
plot_integrated_forest_plot_rna <- function(data_df, outDir_path, cell_type_lookup, levels, colors) {

    # Define the factor levels and labeller for the psychosocial variable
    var_facet_levels <- c("PSS", "ISEL")
    var_facet_labels <- c("PSS" = "Psychological Stress", "ISEL" = "Social Support")
    
    plot_data <- data_df %>%
        mutate(
            Cell_Type = cell_type_lookup[Cluster],
            Log_Odds = log(Odds_Ratio),
            Log_CI_Lower = log(CI_Lower),
            Log_CI_Upper = log(CI_Upper),
            
            # Use explicit direction labels for the y-axis
            Plot_Label = paste0(Cell_Type, " | ", Var_Variable, " ", Var_Target_Dir, " vs LPS ", LPS_Target_Dir),
            
            # Set Cell_Type as factor using the provided levels
            Cell_Type_Factor = factor(Cell_Type, levels = levels[levels %in% Cell_Type]),
            
            # Set Var_Variable as a factor for controlled facet ordering (PSS top, ISEL bottom)
            Var_Variable_Factor = factor(Var_Variable, levels = var_facet_levels)
        ) %>%
        # Order the final Y-axis by Cell Type Factor, then Variable, then Direction
        arrange(Var_Variable_Factor, Cell_Type_Factor, Test_Name) %>%
        # Reverse the final y-axis for plotting order
        mutate(Test_Label = factor(Plot_Label, levels = rev(unique(Plot_Label)))) 

    plot_title <- "Integrated Fisher's Exact Test Enrichment for DEGs (Cell Type x Variable)"
    plot_x_label <- "Log Odds Ratio"

    p_forest <- ggplot(plot_data, aes(x = Log_Odds, y = Test_Label, color = Cell_Type_Factor)) +
        # Add the vertical line at Log Odds = 0 (Odds Ratio = 1, the null hypothesis)
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +

        # Add the confidence intervals as error bars (horizontal)
        geom_errorbarh(aes(xmin = Log_CI_Lower, xmax = Log_CI_Upper), height = 0.2) +

        # Add the point estimate (Log Odds Ratio)
        geom_point(size = 3) + 

        # Facet by the new factored variable and apply custom labels
        facet_grid(
            rows = vars(Var_Variable_Factor), # Use the new factored variable
            scales = "free_y", 
            space = "free_y", 
            switch = "y",
            labeller = labeller(Var_Variable_Factor = var_facet_labels) # Apply the descriptive labels
        ) + 

        # Apply custom colors based on Cell Type Factor
        scale_color_manual(
            name = "Cell Type", 
            values = colors
        ) +

        # Customize labels and theme
        labs(
            x = plot_x_label,
            y = "Cell Type | Comparison",
            title = plot_title,
            subtitle = paste0("Significance threshold (padj) = ", ALPHA, " | Minimum DEGs per cluster = ", MIN_DEGS_THRESHOLD),
            color = "Cell Type"
        ) +
        theme_bw() + # Clean background
        theme(
            # Apply increased font sizes
            plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
            plot.subtitle = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 10), 
            legend.position = "bottom",
            # Increase facet label size (strip text)
            strip.text.y = element_text(size = 14, face = "bold")
        )

    # Save the plot
    fname_plot <- paste0(outDir_path, "25.4.forest_plot_log_odds_integrated_rna_v3.png")
    # Save a larger plot to accommodate all test results
    ggsave(fname_plot, plot = p_forest, width = 12, height = 18, units = "in", dpi = 300)

    cat(paste0("\nSaved integrated forest plot to: ", fname_plot, "\n"))
}

# 4.1 Generate Integrated Plot
plot_integrated_forest_plot_rna(
    final_results_deg, 
    outDir_save, 
    rna_cell_type_map, 
    rna_cell_type_levels, 
    rna_colors
)







rm(list=ls())

##################################################### DAMS #####################################
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(stringr)

# Define the significance threshold
ALPHA <- 0.1

# --- 0. Initial Setup and Data Loading ---

# Set directories
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/LPS_enrich/DAMs/")
# Ensure the output directories are correctly set based on the user's snippet
outDir_save <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/LPS_enrich/DAMs/contingency"


# --- 0.1 LOAD AND PREPARE ATAC-SEQ DATA (DAMS) ---

# Load Psychosocial DAMS (vardams) - uses the user's updated path
fname_vardams <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
vardams <- read.table(fname_vardams, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

# Add identifier (gene/motif name)
vardams$identifier <- vardams$gene

# Load LPS DAMS (treats_dam) - uses the user's updated path
#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_treat/2_th20_plotData.comb.txt.gz"
fname_treats_dam = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_treat/2_th20_plotData.comb.txt.gz"
treats_dam <- read.table(fname_treats_dam, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
treats_dam$identifier <- treats_dam$gene


# Filter out A8 DC (Cluster C8) from LPS data as per user's original code
#treats_dam <- treats_dam %>% filter(Cluster != "C10")
treats_dam <- treats_dam %>% filter(Cluster %in% c("C0", "C1", "C2", "C4", "C6"))


# Define all clusters, variables, and directional combinations
# Note: These are the clusters remaining after filtering C8 from treats_dam, as implied by the user's code
clusters_to_analyze <- unique(vardams$Cluster)
clusters_to_analyze <- clusters_to_analyze[clusters_to_analyze %in% c("C0", "C1", "C2", "C4", "C6")] # Also ensure C8 is out of PS
clusters_to_analyze <- sort(clusters_to_analyze)

ps_vars <- c("PSS_all_mean", "ISEL_Mean")
# Short names used internally in data processing
ps_vars_short <- c("PSS", "ISEL")

# Define cell type labels for clusters
cell_type_map <- c(
    "C0" = "A0 T CD4+",
    "C5" = "A5 T CD4+",
    "C6" = "A6 B",
    "C7" = "A7 T CD4+",
    "C1" = "A1 T CD8+",
    "C2" = "A2 NK",
    "C3" = "A3 T CD4+",
    "C4" = "A4 Monocyte",
    "C8" = "A8 T CD4+", 
    "C9" = "A9 T CD4+"#, 
    #"C10" = "A8 T CD4+"

)


# List of directional pairs: Var_DIR vs LPS_DIR
directions <- list(
    c("Up", "Up"),  # Up vs Up
    c("Down", "Up"), # Down vs Up
    c("Up", "Down"), # Up vs Down
    c("Down", "Down")# Down vs Down
)

# --- NEW: Define factor levels and color palette for plotting ---
cell_type_levels <- c("A0 T CD4+", "A1 T CD8+", "A2 NK", "A4 Monocyte", "A6 B", "A5 T CD4+", "A3 T CD4+", "A7 T CD4+", "A10 DC")
col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#AA4B56","#FF7F00","#FF7F00","#FF7F00", "#D4B9DA")
names(col2) <- cell_type_levels

# --- NEW: Define Psychosocial Variable Labels and Plotting Order ---
ps_labels_map <- c(
    "PSS" = "Psychological Stress",
    "ISEL" = "Social Support"
)

# Define the desired facet order for plotting (Stress on top)
ps_plot_levels <- c("Psychological Stress", "Social Support")

# --- 1.5 Determine Variable/Cluster Pairs with at least one Significant DAM (padj_t < ALPHA) ---

# Filter the raw psychosocial DAMs to find pairs that have at least one significant hit (padj_t < 0.1)
damssig_filter <- vardams %>% 
    filter(padj_t < ALPHA) %>%
    # Select the columns that define the pair (Variable and Cluster)
    select(psycho_variable, Cluster) %>%
    unique() %>%
    # Exclude C8 clusters from the filter list
    filter(Cluster != "C8") %>%
    filter(Cluster != "C3")


# --- 1. Core Fisher Test Function (Modified to use 'identifier' as key) ---

#' Runs a two-sided Fisher's Exact Test for enrichment on cluster-specific data.
#' NOTE: Assumes input dataframes (ps_data and treats_data) already have
#' 'padj' and 'logFC' and are already filtered to a single cluster.
run_fisher_test_dynamic <- function(ps_data, ps_dir_sign, lps_dir_sign, test_name, treats_data, cluster_name, ps_variable) {

  # --- A. Create Binary Status for PS Variable ---
  ps_status <- ps_data %>%
    select(identifier, padj, logFC) %>%
    rename(padj_Var = padj, logFC_Var = logFC) %>%
    mutate(
      Var_Status = ifelse(
        # Check for significant and correct direction
        padj_Var < ALPHA & (if (ps_dir_sign == "Up") logFC_Var > 0 else logFC_Var < 0),
        TRUE,
        FALSE
      )
    )

  # --- B. Create Binary Status for LPS Treatment ---
  lps_status <- treats_data %>%
    select(identifier, padj, logFC) %>%
    rename(padj_LPS = padj, logFC_LPS = logFC) %>%
    mutate(
      LPS_Status = ifelse(
        # Check for significant and correct direction
        padj_LPS < ALPHA & (if (lps_dir_sign == "Up") logFC_LPS > 0 else logFC_LPS < 0),
        TRUE,
        FALSE
      )
    )

  # --- C. Merge Datasets (Common Universe of Motifs within this Cluster) ---
  merged_status <- inner_join(
    ps_status,
    lps_status,
    by = "identifier"
  )

  # --- D. Calculate the Four Cell Counts (A, B, C, D) ---
  A <- merged_status %>% filter(Var_Status == TRUE & LPS_Status == TRUE) %>% nrow()
  B <- merged_status %>% filter(Var_Status == TRUE & LPS_Status == FALSE) %>% nrow()
  C <- merged_status %>% filter(Var_Status == FALSE & LPS_Status == TRUE) %>% nrow()
  D <- merged_status %>% filter(Var_Status == FALSE & LPS_Status == FALSE) %>% nrow()

  # --- E. Construct the Contingency Table ---
  # Use descriptive labels for the matrix rows and columns
  var_target_label = paste0(ps_variable, "_", ps_dir_sign, "_DAMs") # e.g., "PSS_Up_DAMs"
  var_not_target_label = paste0("NOT_", ps_variable, "_", ps_dir_sign, "_DAMs") # e.g., "NOT_PSS_Up_DAMs"
  lps_target_label = paste0("LPS_", lps_dir_sign, "_DAMs") # e.g., "LPS_Up_DAMs"
  lps_not_target_label = paste0("NOT_LPS_", lps_dir_sign, "_DAMs") # e.g., "NOT_LPS_Up_DAMs"

  contingency_matrix <- matrix(
    c(A, C, B, D), # Fill column-wise
    nrow = 2,
    dimnames = list(
      # Rows: Variable Status
      Var_Status = c(var_target_label, var_not_target_label),
      # Columns: LPS Status
      LPS_Status = c(lps_target_label, lps_not_target_label)
    )
  )

  # --- F. Run Two-Sided Fisher's Exact Test ---
  fisher_result <- tryCatch({
    fisher.test(contingency_matrix, alternative = "two.sided")
  }, error = function(e) {
    # Return dummy result if test fails (e.g., zero count)
    list(p.value = NA, estimate = NA, conf.int = c(NA, NA))
  })


  # --- G. Format Results ---
  results_df <- data.frame(
    Test_Name = test_name,
    Cluster = cluster_name,
    Var_Variable = ps_variable,
    Var_Target_Dir = ps_dir_sign,
    LPS_Target_Dir = lps_dir_sign,
    N_Enrichment = A,
    N_Universe = A + B + C + D,
    P_Value = fisher_result$p.value,
    Odds_Ratio = fisher_result$estimate,
    CI_Lower = fisher_result$conf.int[1],
    CI_Upper = fisher_result$conf.int[2],
    stringsAsFactors = FALSE
  )

  return(list(stats = results_df, matrix = contingency_matrix))
}


# Initialize lists to store results
all_stats_list <- list()
all_matrices_list <- list()
counter <- 1

# --- 2. Loop Through Clusters, Variables, and Directions (Filtered Tests) ---

# Loop over each cluster
for (cluster in clusters_to_analyze) {

    # 2.1 Filter LPS data for the current cluster once
    lps_subset_raw <- treats_dam %>%
        filter(Cluster == cluster) %>%
        # Rename LPS DAMs columns (beta -> logFC)
        rename(logFC = beta) %>%
        # Select relevant columns for the function
        select(identifier, padj, logFC, Cluster)

    # Loop over each psychosocial variable (PSS and ISEL)
    for (i in seq_along(ps_vars)) {
        ps_var <- ps_vars[i]
        ps_var_short <- ps_vars_short[i]

        # NEW: Skip if this Variable-Cluster pair has no significant DAMs (padj_t < ALPHA)
        current_pair_check <- damssig_filter %>%
            filter(Cluster == cluster, psycho_variable == ps_var) %>%
            nrow()

        if (current_pair_check == 0) {
            # Skip this iteration of the variable loop (and its 4 directional tests)
            cat(paste("Skipping:", cluster, "|", ps_var_short, " - No significant DAMs (padj <", ALPHA, ") found.\n"))
            next
        }

        # 2.2 Filter Var data for the current variable AND cluster
        ps_subset_raw <- vardams %>%
            filter(psycho_variable == ps_var, Cluster == cluster) %>%
            # Rename Var DAMs columns (padj_t -> padj, estimate -> logFC)
            rename(padj = padj_t, logFC = estimate) %>%
            # Select relevant columns for the function
            select(identifier, padj, logFC, Cluster)

        # Loop over the 4 directional combinations
        for (dir_pair in directions) {
            ps_dir <- dir_pair[1]
            lps_dir <- dir_pair[2]

            # Test Name construction
            test_name <- paste0("DAM_", cluster, "_", ps_var_short, ps_dir, "_LPS", lps_dir)

            cat(paste("Running Test", counter, ":", test_name, "\n"))

            # Run the dynamic test function
            test_output <- run_fisher_test_dynamic(
                ps_data = ps_subset_raw,
                ps_dir_sign = ps_dir,
                lps_dir_sign = lps_dir,
                test_name = test_name,
                treats_data = lps_subset_raw,
                cluster_name = cluster,
                ps_variable = ps_var_short
            )

            # Store the results
            all_stats_list[[counter]] <- test_output$stats

            # Store the matrix separately
            all_matrices_list[[test_name]] <- test_output$matrix
            counter <- counter + 1
        }
    }
}

# --- 3. Final Output and Saving ---

# Combine all statistics into a single dataframe
final_results_dam <- bind_rows(all_stats_list)

# 3.1 Save Combined DAM Results
print("--- SUMMARY OF ALL CELL-TYPE-SPECIFIC DAM FISHER'S EXACT TESTS (ATAC-SEQ) ---")
print(final_results_dam)
fname_stats_dam <- paste0(outDir_save, "25.1.summary_DAM_fisher_tests_by_cluster.txt")
write.table(final_results_dam, file = fname_stats_dam, sep = "\t", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)
cat(paste0("\nSaved combined DAM summary statistics to: ", fname_stats_dam, "\n"))


# 3.2 Save all contingency matrices individually
for (test_name in names(all_matrices_list)) {
    mat <- all_matrices_list[[test_name]]
    # Ensure the directory structure from the user's snippet is used
    fname_matrix <- paste0(outDir_save, "/", "matrix_", test_name, ".txt")

    df_to_save <- as.data.frame(mat)
    df_to_save <- cbind(Variable_Status = rownames(mat), df_to_save)
    write.table(df_to_save, file = fname_matrix, sep = "\t", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE)
}
cat(paste0("Saved ", length(all_matrices_list), " individual contingency matrices to: ", outDir_save, "matrix_DAM_*.txt\n"))


# --- 4. Generate Integrated Forest Plot (Log Odds Ratio) ---

# Function to generate and save the integrated forest plot
plot_integrated_forest_plot <- function(data_df, outDir_path, cell_type_lookup, levels, colors, ps_label_map, ps_plot_order) {

    plot_data <- data_df %>%
        mutate(
            # NEW: Translate short variable names to long, descriptive labels
            Var_Variable_Long = ps_label_map[Var_Variable],
            
            # NEW: Create a factor for the variable order to ensure Stress is on top
            Var_Factor = factor(Var_Variable_Long, levels = ps_plot_order),

            Cell_Type = cell_type_lookup[Cluster],
            Log_Odds = log(Odds_Ratio),
            Log_CI_Lower = log(CI_Lower),
            Log_CI_Upper = log(CI_Upper),
            
            # UPDATED: Use long label in Plot_Label
            Plot_Label = paste0(Cell_Type, " | ", Var_Variable_Long, " ", Var_Target_Dir, " vs LPS ", LPS_Target_Dir),
            
            # NEW: Set Cell_Type as factor using the provided levels (only keep existing types)
            Cell_Type_Factor = factor(Cell_Type, levels = levels[levels %in% Cell_Type])
        ) %>%
        # Order the final Y-axis by the new variable factor (Stress first), then Cell Type Factor, then Test Name
        arrange(Var_Factor, Cell_Type_Factor, Test_Name) %>%
        # Reverse the final y-axis for plotting order
        mutate(Test_Label = factor(Plot_Label, levels = rev(unique(Plot_Label)))) 

    plot_title <- "Integrated Fisher's Exact Test Enrichment for DAMs (Cell Type x Variable)"
    plot_x_label <- "Log Odds Ratio"

    p_forest <- ggplot(plot_data, aes(x = Log_Odds, y = Test_Label, color = Cell_Type_Factor)) +
        # Add the vertical line at Log Odds = 0 (Odds Ratio = 1, the null hypothesis)
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +

        # Add the confidence intervals as error bars (horizontal)
        geom_errorbarh(aes(xmin = Log_CI_Lower, xmax = Log_CI_Upper), height = 0.2) +

        # Add the point estimate (Log Odds Ratio)
        geom_point(size = 3) + # Slightly increased point size

        # UPDATED: Facet by the new Var_Factor to control order and display full label
        facet_grid(rows = vars(Var_Factor), scales = "free_y", space = "free_y", switch = "y") + 

        # NEW: Apply custom colors based on Cell Type Factor
        scale_color_manual(
            name = "Cell Type", 
            values = colors
        ) +

        # Customize labels and theme
        labs(
            x = plot_x_label,
            y = "Cell Type | Comparison",
            title = plot_title,
            subtitle = paste0("Significance threshold (padj) = ", ALPHA),
            color = "Cell Type"
        ) +
        theme_bw() + # Clean background
        theme(
            # UPDATED: Increase Font Sizes
            plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
            plot.subtitle = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 10), 
            legend.position = "bottom",
            # UPDATED: Increase facet label size (strip text)
            strip.text.y = element_text(size = 14, face = "bold")
        )

    # Save the plot
    fname_plot <- paste0(outDir_path, "25.2.forest_plot_log_odds_integrated_by_celltype_v3_final.png")
    # Save a larger plot to accommodate all test results
    ggsave(fname_plot, plot = p_forest, width = 12, height = 18, units = "in", dpi = 300)

    cat(paste0("\nSaved integrated forest plot to: ", fname_plot, "\n"))
}

# 4.1 Generate Integrated Plot with new labels, levels, and colors
plot_integrated_forest_plot(
    final_results_dam, 
    outDir_save, 
    cell_type_map, 
    cell_type_levels, 
    col2,
    ps_labels_map,
    ps_plot_levels # Passes the new order (Stress on top)
)


