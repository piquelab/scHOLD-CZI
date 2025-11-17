library(scales)
library(pheatmap)
require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(dplyr)
library(tidyr)
library(sva)

###############################################################################
##################################### add the band. add cel type band too. 

rm(list=ls())

outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/heatmap/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

# Variables to loop over
variables <- c("cytocomp", "PSS_all_mean", "ISEL_Mean") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
variable_names <- c("Cytokines", "Psychological Stress", "Social Support") #) #"socioeconomic status", "diastolic blood pressure", 

treatment <- "CTRL"

cluster_number <- c("0","1","2", "3", "4")#, "5")#, "6")#, "7")
clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R5 DC")

 #col2 <- c("#E69F00","#984EA3","#D55E00", "#CC79A7", "#56B4E9")
col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#D97986")#, "#D4B9DA")

names(col2) <- celltypes

#cor_results <- list()

    # Create a list to store z-scores for each variable
    zscore_list <- list()
    gene_list <- list()
    reslist <- list()
    rescat_list <- list()

    expanded_var <- c()
    expanded_var_name <- c()

        for(i in 1:length(variables)) {
        #tryCatch({

        myvar <- variables[i]
        var_name <- variable_names[i]
            cat("########################################################\n")
            cat("##Processing: ", var_name, "\n")

            for(j in 1:length(cluster_number)){
                #tryCatch({
                cell <- cluster_number[j]
                celltype <- celltypes[j]
               cat("################## Processing cluster: ", cell, "##############","\n")
                        #plotsDir <- paste(platePrefix,"/cluster_",cell, sep='')
                        #system(paste("mkdir -p", plotsDir))
                        ##
                               
                                # load DESeq results - Cindy
                    fname=paste(dataDir, platePrefix,"C",cell, ".deseqres_", myvar, plateSufix, ".txt", sep="")
                    resC <- read.table(fname,  sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
                    resC <- resC %>% dplyr::select(-any_of("baseMean"))
                    #fname=paste(v2dir, statsDirV2, '/', platePrefixV2,"_",myvar, "_", cell,  "_significant_DESeq_results", ".txt", sep="")
                    ressigC <- resC %>% filter(padj < 0.1)
                        cat("# of genes tested in Cindy for - ", myvar, " in cell type ", cell," :", nrow(resC)," \n")
                        cat("# of DEGs in Cindy for - ", myvar, " in cell type ", cell," :", nrow(ressigC)," \n")

                # assign z scores
                resC$zscore <- resC$logFC / resC$SE

                # merge results
                #resC$symb <- paste0(resC$cluster, "_", resC$identifier)
                resC$symb <- resC$identifier

                #resm <- merge(resV3, resC, by="symb", suffixes = c("_A", "_C"))
                #ressigm <- resm %>% filter(padj_A < 0.1 | padj_C < 0.1)
                reslist[[j]] <- resC
     
            rescat <- do.call(rbind, reslist)
             #rescat_list[[i]] <- rescat

            # Append to the zscore_list and gene_list
            zscore_list[[paste0(cell, "_", myvar)]] <- resC$zscore
            gene_list[[paste0(cell, "_", myvar)]] <- resC$symb
            
            # Add expanded name to the list
            expanded_var <- c(expanded_var, paste0(cell, "_", myvar))
            #expanded_var_name <- c(expanded_var_name, paste0(cell, "_", var_name))
            expanded_var_name <- c(expanded_var_name, paste0(var_name, " - ", celltype))

           }

            # Store z-scores and gene names
            #zscore_list[[i]] <- rescat$zscore  # Store z-scores 
            #gene_list[[i]] <- rescat$symb         # Store gene symbols
   }

             #results_all <- do.call(rbind, rescat_list)

    # Initialize matrices for correlations, p-values, and color mapping
    zscore_matrix <- matrix(NA, nrow = length(expanded_var), ncol = length(expanded_var))
    p_matrix <- matrix(NA, nrow = length(expanded_var), ncol = length(expanded_var))
    color_matrix <- matrix(NA, nrow = length(expanded_var), ncol = length(expanded_var))  # New matrix for color coding
    rownames(zscore_matrix) <- expanded_var_name
    colnames(zscore_matrix) <- expanded_var_name

    text_colors <- matrix("black", nrow = length(expanded_var), ncol = length(expanded_var))

    # Pairwise correlation loop (include only common genes for each pair of variables)
    for (v1 in 1:length(expanded_var)) {
        for (v2 in 1:length(expanded_var)) {
            # Get the common genes between the two variables
            common_genes <- intersect(gene_list[[v1]], gene_list[[v2]])
            
            # Filter z-scores for the common genes
            zscores_v1 <- zscore_list[[v1]][gene_list[[v1]] %in% common_genes]
            zscores_v2 <- zscore_list[[v2]][gene_list[[v2]] %in% common_genes]
            
            # Calculate Pearson correlation and p-value
            cor_test <- cor.test(zscores_v1, zscores_v2, method = "pearson")
            zscore_matrix[v1, v2] <- cor_test$estimate
            p_matrix[v1, v2] <- cor_test$p.value  # Store the p-value

            # Set the text color to white if correlation is > 0.60 or < -0.60
            if (cor_test$estimate > 0.60 | cor_test$estimate < -0.60) {
                text_colors[v1, v2] <- "white"
            }

            # Set color to white for non-significant correlations (p > 0.05)
            if (cor_test$p.value > 0.05) {
                color_matrix[v1, v2] <- NA  # Assign NA for non-significant correlations
            } else {
                color_matrix[v1, v2] <- zscore_matrix[v1, v2]  # Retain correlation value for significant correlations
            }

            # Mask non-significant correlations, replace them with zero. 
            if (cor_test$p.value > 0.05) {
                zscore_matrix[v1, v2] <- 0  # Assign NA for non-significant correlations
            } else {
                zscore_matrix[v1, v2] <- zscore_matrix[v1, v2]  # Retain correlation value for significant correlations
            }



        }
    }

colnames(zscore_matrix)

# Save the combined dataframe as a single text file
fname <- paste0(getwd(), "/02.6.zscore-correlations_cytokines.txt")
write.table(zscore_matrix, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = T)

#celltype <- c("R0 CD4+ T cell", "R1 CD4+ T cell", "R2 CD8+ T cell", "R3 NK cell", "R4 Monocyte", "R5 B cell")#, "R6 Dendritic cell")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R5 DC")

 #col2 <- c("#E69F00","#984EA3","#D55E00", "#CC79A7", "#56B4E9")
col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#D97986")#, "#D4B9DA")

names(col2) <- celltype

expanded_annotations <- rep(variable_names, each = length(celltypes))
#expanded_var_name <- paste0(uniondegs$Variable, " - ", uniondegs$celltype)
#expanded_var_name <- unique(expanded_var_name)

# Create row and column annotations based on the expanded variable names
row_annotation <- data.frame(Variable = expanded_annotations)

rownames(row_annotation) <- expanded_var_name
row_annotation$CellType <- paste0(gsub(".* - ", "", rownames(row_annotation)))
col_annotation <- row_annotation

#var_color <- c("#E34234", "#29465B")#, "#FFB6C1")#, "#E34234", "#ADD8E6",) 
var_color <- c("#FFA500", "#E34234", "#3A84D8")#, "#FFB6C1")#, "#E34234", "#ADD8E6",) "#CC5500", "#C1440E", "#A23E00",
#var_color <- c("#FFA500", "#FFB733" "#E34234", "#3A84D8")#, "#FFB6C1")#, "#E34234", "#ADD8E6",) "#CC5500", "#C1440E", "#A23E00",

 #cell_colors <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#D97986")
cell_colors <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#D97986")#, "#D4B9DA")

names(var_color) <- variable_names
names(cell_colors) <- celltypes

annotation_colors <- list(Variable  = var_color, CellType = cell_colors)


   # Define a custom color palette that includes green for non-significant correlations
    color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

    # Plot the clustered heatmap

    heatmap_filename_clustered <- paste(outFolder, "02.1.zscore-correlation_allgenes_heatmap_3var_5cell_annot.png", sep = "")

    png(heatmap_filename_clustered, width = 2300, height = 1500, res = 240)
    #png(heatmap_filename_clustered, width = 5000, height = 4000, res = 240)

    # Set the non-significant correlations to be displayed in green
    p2 <- pheatmap(zscore_matrix,
                   main = paste0("z-score correlation"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("blue", "white", "red"))(100),  # Red for negative, blue for positive
                   breaks = seq(-1, 1, length.out = 101),
                   display_numbers = TRUE,
                   #number_color = text_colors,
                   fontsize = 10, fontsize_number = 7,               
                   #show_colnames = FALSE,  # Removes column labels
                   na_col = "white", # Non-significant correlations are white
                   annotation_row = row_annotation,  # Add row annotations
                   #annotation_col = col_annotation,  # Add column annotations
                   annotation_colors = annotation_colors # Add colors for annotations
                   )  

    # Save the clustered heatmap
    print(p2)
    dev.off()

#figfn <- paste0(outFolder, "01.4.heatmap_logFC_unionDEGs_2Var_5celltypes_annotations_all_genes_scalerowoff_celltypeorder.pdf")
figfn <- paste0(outFolder, "02.1.zscore-correlation_allgenes_heatmap_3var_5cell_annot.pdf")
pdf(figfn, width=14, height=7)
print(p2)
dev.off()


#simplified
        heatmap_filename_clustered <- paste(outFolder, "02.2.zscore-correlation_allgenes_heatmap_3var_5cell-simplified.png", sep = "")
        #heatmap_filename_clustered <- paste("05.3.zscore-correlation_allgenes_heatmap_3var-6celltypes-simplified_proinf.png", sep = "")

    png(heatmap_filename_clustered, width = 2300, height = 1000, res = 240)
    #png(heatmap_filename_clustered, width = 3500, height = 2500, res = 240)
    
    # Set the non-significant correlations to be displayed in green
    p2 <- pheatmap(zscore_matrix,
                   main = paste0("z-score correlation"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("blue", "white", "red"))(100),  # Red for negative, blue for positive
                   breaks = seq(-1, 1, length.out = 101),
                   display_numbers = FALSE,
                   #number_color = text_colors,
                   fontsize = 10, fontsize_number = 5.5,    
                   show_colnames = FALSE,  # Removes column labels
                   show_rownames = TRUE,   # Ensures row labels are shown           
                   #show_colnames = FALSE,  # Removes column labels
                   na_col = "white", # Non-significant correlations are white
                   annotation_row = row_annotation,  # Add row annotations
                   #annotation_col = col_annotation,  # Add column annotations
                   annotation_colors = annotation_colors # Add colors for annotations
                   )  

    # Save the clustered heatmap
    print(p2)
    dev.off()


#figfn <- paste0(outFolder, "01.4.heatmap_logFC_unionDEGs_2Var_5celltypes_annotations_all_genes_scalerowoff_celltypeorder.pdf")
figfn <- paste0(outFolder, "02.2.zscore-correlation_allgenes_heatmap_3var_5cell-simplified.pdf")
pdf(figfn, width=14, height=7)
print(p2)
dev.off()
